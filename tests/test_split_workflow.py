import errno
import gzip
import io
import json
import os
import re
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

from flask import Flask
from scripts.split_web import register_split_routes
from scripts.i18n import register_i18n, message
from scripts.split_workflow import InputError, parse_parameters, DEFAULTS, describe_failure, FAILURES
from scripts.utils import create_blast_db
from src.common import paired_scan
from src.common.job_errors import ToolFailedError
from src.common.paired_scan import ScanError, ScanTimeout, scan_pairs
from src.common.paired_specificity import analyze_probe_sequences, load_reference
from scripts.reference_library import ReferenceLibrary
from werkzeug.datastructures import FileStorage


SEQ = 'ACGT' * 20


def fake_design(**kwargs):
    """Only replace the expensive engine; exercise actual HTTP/storage adapters."""
    assessed = bool(kwargs['ref_genome'])
    probe = {'id': 1, 'keep': not assessed,
             'reason_code': 'cross_species_exact_pair' if assessed else 'selected',
             'evaluation': {'quality_score': 5, 'specificity_score': 100 if assessed else None},
             'specificity': {'status': 'completed' if assessed else 'not_evaluated',
                             'on_target_exact': assessed, 'counts': {'cross_species': [1, 1, 1, 1]} if assessed else {},
                             'paired_matches': []}}
    for side, start in [('left_probe', 1), ('right_probe', 16)]:
        probe[side] = {'sequence': SEQ[:15], 'position': {'start': start, 'end': start + 14},
                       'tm': 50, 'gc_content': 50, 'blast_results': {'private_raw_hsp': 'retained in ZIP'}}
    output = Path(kwargs['output_dir'])
    (output / (kwargs['name'] + '_filtered.json')).write_text(json.dumps(
        {'probe_sets': [probe], 'target_policy': kwargs['target_policy'],
         'search': {'status': 'completed', 'seconds': 1.5, 'reference': 'secret-path', 'executable': 'bin'} if assessed else None}))
    (output / 'dual_probe.csv').write_text('Probe ID,Sequence\n')


class GuidedWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.reference = self.root / 'bol_bju_cdna.fasta'
        self.reference.write_text('>BjuTarget\n' + SEQ + '\n>BoOther\n' + SEQ + '\n')
        repository = Path(__file__).resolve().parents[1]
        self.app = Flask(__name__, template_folder=str(repository / 'templates'), static_folder=str(repository / 'static'))
        # Jobs run inside submit() so each request finishes before the test continues.
        self.app.config.update(TESTING=True, UPLOAD_FOLDER=str(self.root), SPLIT_JOBS_INLINE=True)
        for endpoint in ['index', 'triplet', 'triplet_batch', 'bridge', 'bridge_list', 'bridge_query',
                         'bridge_sequence', 'bridge_analyze', 'bridge_validate', 'formatting', 'subgenomes',
                         'probe_recommend', 'filter', 'help', 'orders', 'references']:
            self.app.add_url_rule('/' + endpoint, endpoint, lambda: 'placeholder')
        register_split_routes(self.app, {'BP_0079': 'A' * 19}, {'BjuTarget': 'T' * 80}, {}, {},
                              lambda: [{'filename': self.reference.name, 'display_name': 'Bo/Bju'}])
        self.client = self.app.test_client()
        self.input = {'input_mode': 'single', 'name': 'ABCG13', 'geneID': 'bjutarget', 'bp_id': 'BP_0079',
                      'sequenceType': 'reference', 'reference_mode': 'existing',
                      'existing_ref_genome': self.reference.name}

    def tearDown(self):
        self.tmp.cleanup()

    def run_job(self, prepared, parameters):
        """Queue a job and return (run response, finished result JSON or None)."""
        response = self.client.post(prepared['run_url'], json=parameters)
        if response.status_code not in (200, 202):
            return response, None
        result = self.client.get(response.json['result_url'], headers={'Accept': 'application/json'})
        return response, result.json if result.status_code == 200 else None

    def prepare(self, **overrides):
        response = self.client.post('/split/prepare', data=dict(self.input, **overrides))
        self.assertEqual(response.status_code, 200, response.get_json())
        return response.get_json()

    def test_single_and_batch_pages_share_guided_defaults(self):
        for url in ['/split', '/split/batch']:
            page = self.client.get(url)
            self.assertEqual(page.status_code, 200)
            self.assertIn('序列与参考', page.text)
            self.assertIn('结果复核', page.text)
            self.assertIn('name="min_gap" value="2"', page.text)
        template = self.client.get('/split/batch-template.csv')
        self.assertIn('task_name,BP_ID,gene_id,max_selected,allowed_ids', template.data.decode('utf-8-sig'))
        template.close()

    def test_bojbu_target_comes_from_selected_reference_case_insensitively(self):
        prepared = self.prepare()
        task = prepared['tasks'][0]
        self.assertEqual(task['gene_id'], 'BjuTarget')
        self.assertNotIn('sequence', task)
        manifest = json.loads((self.root / 'split_results' / prepared['id'] / 'request.json').read_text())
        self.assertEqual(manifest['tasks'][0]['sequence'], SEQ)  # not the conflicting internal cDNA
        self.assertEqual(prepared['reference']['species_counts'], {'Bju': 1, 'Bo': 1})

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_multiple_chromosome_sets_use_saved_library_annotation_and_snapshot(self, design):
        library = ReferenceLibrary(self.root, self.root / 'instance' / 'references.json')
        self.app.extensions['reference_library'] = library
        library.configure(self.reference.name, 'Bo/Bju', True,
                          FileStorage(stream=io.BytesIO(b'BjuTarget,A\nBoOther,C\n'), filename='groups.csv'))
        prepared = self.prepare(process_multiple_sets='on')
        self.assertTrue(prepared['process_multiple_sets'])
        directory = self.root / 'split_results' / prepared['id']
        manifest = json.loads((directory / 'request.json').read_text())
        self.assertTrue(manifest['process_multiple_sets'])
        self.assertEqual(prepared['tasks'][0]['subgenome'], 'A')
        # Subsequent edits in the shared library do not alter a prepared task.
        library.configure(self.reference.name, 'Bo/Bju', True,
                          FileStorage(stream=io.BytesIO(b'BjuTarget,B\nBoOther,C\n'), filename='updated.csv'))
        response, result = self.run_job(prepared, DEFAULTS)
        self.assertEqual(response.status_code, 202, response.json)
        policy = design.call_args.kwargs['target_policy']
        self.assertEqual(policy['mode'], 'subgenome')
        self.assertEqual(policy['target_subgenome'], 'A')
        self.assertEqual(Path(policy['subgenome_map']), (directory / 'subgenomes.tsv').resolve())
        self.assertIn('BjuTarget\tA', Path(policy['subgenome_map']).read_text())
        self.assertEqual(result['parameters']['target_mode'], 'subgenome')
        self.assertEqual(self.client.post(prepared['run_url'], json=DEFAULTS).status_code, 200)
        plain = self.prepare()
        self.assertFalse(plain['process_multiple_sets'])
        self.assertIsNone(plain['tasks'][0]['subgenome'])
        self.assertFalse((self.root / 'split_results' / plain['id'] / 'subgenomes.tsv').exists())
        response, _ = self.run_job(plain, DEFAULTS)
        self.assertEqual(response.status_code, 202)
        self.assertEqual(design.call_args.kwargs['target_policy']['mode'], 'species')
        page = self.client.get('/split').text
        self.assertIn('是否考虑多套染色体', page)
        self.assertNotIn('name="subgenome_map"', page)

    def test_multiple_sets_reject_missing_ploidy_annotation_and_target_membership(self):
        library = ReferenceLibrary(self.root, self.root / 'instance' / 'references.json')
        self.app.extensions['reference_library'] = library
        for ploidy, annotation in [(False, None), (True, None),
                                  (True, FileStorage(stream=io.BytesIO(b'BoOther,C\n'), filename='groups.csv'))]:
            library.configure(self.reference.name, 'Bo/Bju', ploidy, annotation)
            response = self.client.post('/split/prepare', data=dict(self.input, process_multiple_sets='on'))
            self.assertEqual(response.status_code, 400, response.json)
            self.assertEqual(response.json['field'], 'process_multiple_sets')
        self.assertEqual(len(list((self.root / 'split_results').iterdir())), 0)

    def test_pasted_sequence_keeps_target_id_and_is_not_overwritten(self):
        pasted = 'G' * 40
        prepared = self.prepare(sequenceType='custom', seq='>user\n' + pasted.lower())
        manifest = json.loads((self.root / 'split_results' / prepared['id'] / 'request.json').read_text())
        self.assertEqual(manifest['tasks'][0]['sequence'], pasted)
        invalid = self.client.post('/split/prepare', data=dict(self.input, sequenceType='custom', seq=pasted, geneID=''))
        self.assertEqual(invalid.status_code, 400)
        self.assertEqual(invalid.json['field'], 'geneID')

    def test_batch_headers_count_overrides_and_no_silent_skipping(self):
        data = 'task_name,BP_ID,gene_id,max_selected\none,BP_0079,BjuTarget,3\ntwo,BP_0079,BoOther,\n'
        prepared = self.prepare(input_mode='batch', batch_text=data)
        self.assertEqual([t['max_selected'] for t in prepared['tasks']], [3, None])
        invalid = self.client.post('/split/prepare', data=dict(self.input, input_mode='batch',
            batch_text='task_name,BP_ID,gene_id\none,BP_1001,BjuTarget\ntwo,BP_0079,Missing\n'))
        self.assertEqual(invalid.status_code, 400)
        self.assertIn('第 2 行', invalid.json['error'])
        self.assertIn('第 3 行', invalid.json['error'])
        self.assertEqual(len(list((self.root / 'split_results').iterdir())), 1)

    def test_allowed_transcripts_are_validated_against_the_reference(self):
        prepared = self.prepare(allowed_ids='boother; BOOTHER')
        self.assertEqual(prepared['tasks'][0]['allowed_ids'], ['BoOther'])
        prepared = self.prepare(input_mode='batch', batch_text='one,BP_0079,BjuTarget,3,BoOther\ntwo,BP_0079,BoOther,,\n')
        self.assertEqual([t['allowed_ids'] for t in prepared['tasks']], [['BoOther'], []])
        for allowed in ['Missing', 'BjuTarget']:
            response = self.client.post('/split/prepare', data=dict(self.input, allowed_ids=allowed))
            self.assertEqual(response.status_code, 400, allowed)
        self.assertEqual(self.client.post('/split/prepare', data=dict(self.input, input_mode='batch',
            batch_text='one,BP_0079,BjuTarget,3,BoOther,extra\n')).status_code, 400)

    def subgenome_upload(self, name, records):
        return (io.BytesIO(''.join('>%s\n%s\n' % (key, SEQ) for key in records).encode()), name)

    def test_subgenome_files_are_merged_into_one_reference_and_mapping(self):
        prepared = self.prepare(reference_mode='subgenome_files', existing_ref_genome='', geneID='BjuA1',
                                subgenome_files=[self.subgenome_upload('A.fasta', ['BjuA1', 'BjuA2']),
                                                 self.subgenome_upload('B.fa', ['BjuB1'])])
        self.assertEqual(prepared['reference']['sequence_count'], 3)
        self.assertEqual(prepared['reference']['subgenome'],
                         {'source': '按文件划分', 'annotated': 3, 'unannotated': 0, 'outside_reference': 0,
                          'subgenome_counts': {'A': 2, 'B': 1}})
        self.assertEqual(prepared['tasks'][0]['subgenome'], 'A')
        directory = self.root / 'split_results' / prepared['id']
        self.assertIn('BjuB1\tB', (directory / 'subgenomes.tsv').read_text())
        self.assertEqual((directory / 'reference.fasta').read_text().count('>'), 3)
        self.assertEqual(self.reference.read_text().count('>'), 2)  # the shared library is untouched

    def test_subgenome_file_uploads_reject_duplicate_ids_names_and_single_files(self):
        for files, field in [
                ([self.subgenome_upload('A.fasta', ['BjuA1']), self.subgenome_upload('B.fa', ['BjuA1'])], 'subgenome_files'),
                ([self.subgenome_upload('A.fasta', ['BjuA1']), self.subgenome_upload('A.fa', ['BjuB1'])], 'subgenome_files'),
                ([self.subgenome_upload('A.fasta', ['BjuA1', 'BjuA2'])], 'subgenome_files'),
                ([self.subgenome_upload('A.txt', ['BjuA1']), self.subgenome_upload('B.fa', ['BjuB1'])], 'subgenome_files')]:
            response = self.client.post('/split/prepare', data=dict(
                self.input, reference_mode='subgenome_files', existing_ref_genome='', geneID='BjuA1',
                subgenome_files=files))
            self.assertEqual(response.status_code, 400, files)
            self.assertEqual(response.json['field'], field)
        self.assertEqual(len(list((self.root / 'split_results').iterdir())), 0)
        mixed = self.client.post('/split/prepare', data=dict(
            self.input, reference_mode='subgenome_files',
            subgenome_files=[self.subgenome_upload('A.fasta', ['BjuA1']), self.subgenome_upload('B.fa', ['BjuB1'])]))
        self.assertEqual(mixed.status_code, 400)

    def test_subgenome_mapping_file_is_validated_against_the_reference(self):
        mapping = b'transcript_id,subgenome\nBjuTarget,A\nBoOther,B\nAbsent,C\n'
        prepared = self.prepare(subgenome_map=(io.BytesIO(mapping), 'subgenomes.csv'))
        self.assertEqual(prepared['reference']['subgenome']['subgenome_counts'], {'A': 1, 'B': 1})
        self.assertEqual(prepared['reference']['subgenome']['outside_reference'], 1)
        self.assertEqual(prepared['tasks'][0]['subgenome'], 'A')
        for upload, field in [((io.BytesIO(b'Absent,C\n'), 'subgenomes.csv'), 'subgenome_map'),
                              ((io.BytesIO(b'BjuTarget\n'), 'subgenomes.csv'), 'subgenome_map'),
                              ((io.BytesIO(mapping), 'subgenomes.pdf'), 'subgenome_map')]:
            response = self.client.post('/split/prepare', data=dict(self.input, subgenome_map=upload))
            self.assertEqual(response.status_code, 400)
            self.assertEqual(response.json['field'], field)

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_subgenome_mode_needs_annotation_and_an_annotated_target(self, design):
        plain = self.prepare()
        response = self.client.post(plain['run_url'], json=dict(DEFAULTS, target_mode='subgenome'))
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json['field'], 'target_mode')
        self.assertIn('第一步', response.json['error'])
        annotated = self.prepare(input_mode='batch', batch_text='one,BP_0079,BjuTarget\ntwo,BP_0079,BoOther\n',
                                 subgenome_map=(io.BytesIO(b'BjuTarget\tA\n'), 'subgenomes.tsv'))
        response = self.client.post(annotated['run_url'], json=dict(DEFAULTS, target_mode='subgenome'))
        self.assertEqual(response.status_code, 400)
        self.assertIn('two', response.json['error'])
        full = self.prepare(subgenome_map=(io.BytesIO(b'BjuTarget\tA\nBoOther\tB\n'), 'subgenomes.tsv'))
        response, result = self.run_job(full, dict(DEFAULTS, target_mode='subgenome', allowed_ids=''))
        self.assertEqual(response.status_code, 202, response.json)
        policy = design.call_args.kwargs['target_policy']
        self.assertEqual(policy['mode'], 'subgenome')
        self.assertEqual(policy['target_subgenome'], 'A')
        self.assertTrue(policy['subgenome_map'].endswith('subgenomes.tsv'))
        self.assertEqual(result['tasks'][0]['target_policy'], policy)

    def test_target_mode_parameter_is_validated(self):
        self.assertEqual(parse_parameters({})['target_mode'], 'species')
        self.assertEqual(parse_parameters({'target_mode': 'subgenome'})['target_mode'], 'subgenome')
        self.assertEqual(parse_parameters({'target_mode': 'transcript', 'exclusion_max_mismatches': 2})['exclusion_max_mismatches'], 2)
        for overrides in [{'target_mode': 'copy'}, {'exclusion_max_mismatches': 4}]:
            with self.assertRaises(InputError):
                parse_parameters(overrides)

    def test_ambiguous_target_and_duplicate_task_names_are_rejected(self):
        self.reference.write_text('>BjuTarget.1\n' + SEQ + '\n>BjuTarget.2\n' + SEQ + '\n')
        response = self.client.post('/split/prepare', data=self.input)
        self.assertEqual(response.status_code, 400)
        self.assertIn('多个转录本', response.json['error'])
        response = self.client.post('/split/prepare', data=dict(self.input, input_mode='batch',
            batch_text='same,BP_0079,BjuTarget.1\nSAME,BP_0079,BjuTarget.2'))
        self.assertIn('重复', response.json['error'])

    def test_uploaded_gzip_is_staged_once_without_overwriting_existing_reference(self):
        original = self.reference.read_bytes()
        compressed = gzip.compress(('>BjuUpload\n' + SEQ).encode())
        response = self.client.post('/split/prepare', data=dict(self.input, reference_mode='upload',
            existing_ref_genome='', geneID='BjuUpload', ref_genome=(io.BytesIO(compressed), 'bol_bju_cdna.fasta.gz')))
        self.assertEqual(response.status_code, 200, response.json)
        self.assertEqual(self.reference.read_bytes(), original)
        staged = self.root / 'split_results' / response.json['id'] / 'reference.fasta'
        self.assertIn('BjuUpload', staged.read_text())

    def test_bad_reference_choice_and_unsafe_task_paths_are_rejected(self):
        for overrides in [{'existing_ref_genome': '../outside.fasta'}, {'name': '../overwrite'},
                          {'name': 'a,b'}, {'reference_mode': 'none'}, {'bp_id': 'BP_1001'}]:
            response = self.client.post('/split/prepare', data=dict(self.input, **overrides))
            self.assertEqual(response.status_code, 400, overrides)
        response = self.client.post('/split/prepare', data=dict(self.input,
            ref_genome=(io.BytesIO(b'>BjuTarget\n' + SEQ.encode()), 'ref.fa')))
        self.assertEqual(response.status_code, 400)

    def test_parameter_bounds_and_order_are_checked(self):
        self.assertEqual(parse_parameters({})['min_gap'], 2)
        for overrides in [{'min_tm': 'nan'}, {'max_selected': 0}, {'max_selected': '1.5'},
                          {'min_gap': -1}, {'min_gc': 80, 'max_gc': 60}, {'min_length': 20, 'max_length': 15}]:
            with self.assertRaises(InputError):
                parse_parameters(overrides)

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_full_enumeration_mode_reaches_engine_and_saved_result(self, design):
        self.assertEqual(parse_parameters({})['candidate_mode'], 'first_valid')
        with self.assertRaises(InputError): parse_parameters({'candidate_mode': 'unknown'})
        self.assertIn('name="candidate_mode"', self.client.get('/split/batch').text)
        prepared = self.prepare()
        response, result = self.run_job(prepared, dict(DEFAULTS, candidate_mode='exhaustive'))
        self.assertEqual(response.status_code, 202)
        self.assertEqual(design.call_args.kwargs['candidate_mode'], 'exhaustive')
        self.assertEqual(result['parameters']['candidate_mode'], 'exhaustive')
        summary = self.client.get(result['result_url'] + '?summary=1', headers={'Accept': 'application/json'}).json
        self.assertTrue(summary['lazy_tasks'])
        self.assertNotIn('probe_sets', summary['tasks'][0])
        task = self.client.get(summary['task_url_prefix'] + '0')
        self.assertEqual(task.status_code, 200)
        self.assertEqual(task.json, result['tasks'][0])
        self.assertEqual(self.client.get(summary['task_url_prefix'] + '13').status_code, 404)
        archive = self.client.get(result['download_url'])
        with ZipFile(io.BytesIO(archive.data)) as z:
            self.assertEqual(json.loads(z.read('run_summary.json'))['parameters']['candidate_mode'], 'exhaustive')
        archive.close()

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_saved_results_download_and_duplicate_run_are_consistent(self, design):
        prepared = self.prepare()
        response, result = self.run_job(prepared, DEFAULTS)
        self.assertEqual(response.status_code, 202, response.json)
        self.assertEqual(response.json['status'], 'completed')
        self.assertEqual(result['tasks'][0]['cross_exact'], 1)
        self.assertEqual(result['tasks'][0]['selected'], 0)
        self.assertEqual(result['tasks'][0]['excluded'], 1)
        self.assertTrue(result['search_exhaustive'])
        self.assertEqual(result['tasks'][0]['search'], {'status': 'completed', 'seconds': 1.5})
        self.assertEqual(result['tasks'][0]['target_policy'], {'mode': 'species', 'allowed_ids': [], 'exclusion_max_mismatches': 0})
        self.assertEqual(design.call_args.kwargs['target_policy'], {'mode': 'species', 'allowed_ids': [], 'exclusion_max_mismatches': 0})
        self.assertNotIn('blast_results', result['tasks'][0]['probe_sets'][0]['left_probe'])
        self.assertEqual(self.client.get(result['result_url']).status_code, 200)
        archive = self.client.get(result['download_url'])
        with ZipFile(io.BytesIO(archive.data)) as z:
            self.assertIn('run_summary.json', z.namelist())
            full = json.loads(z.read('ABCG13_filtered.json'))
            self.assertIn('blast_results', full['probe_sets'][0]['left_probe'])
        archive.close()
        self.assertEqual(self.client.post(prepared['run_url'], json=DEFAULTS).status_code, 200)  # already done
        # Completed jobs made before candidate_mode existed keep their original behavior.
        result_path = self.root / 'split_results' / prepared['id'] / 'result.json'
        old_result = json.loads(result_path.read_text())
        old_result['parameters'].pop('candidate_mode')
        result_path.write_text(json.dumps(old_result))
        self.assertEqual(self.client.post(prepared['run_url'], json=DEFAULTS).status_code, 200)
        self.assertEqual(design.call_count, 1)
        self.assertEqual(self.client.post(prepared['run_url'], json=dict(DEFAULTS, min_gap=3)).status_code, 409)

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_no_reference_stays_unassessed_and_legacy_post_still_downloads(self, design):
        response = self.client.post('/split', data=dict(self.input, reference_mode='none',
            existing_ref_genome='', sequenceType='custom', seq=SEQ, geneID=''))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.mimetype, 'application/zip')
        with ZipFile(io.BytesIO(response.data)) as z:
            report = json.loads(z.read('run_summary.json'))
            self.assertEqual(report['tasks'][0]['unassessed'], 1)
            self.assertIsNone(report['reference'])
            self.assertFalse(report['search_exhaustive'])
        response.close()
        batch = self.client.post('/split/batch', data=dict(self.input, target_mode='transcript', exclusion_max_mismatches='1',
            batch_file=(io.BytesIO(b'one,BP_0079,BjuTarget,2,BoOther\n'), 'batch.csv')))
        self.assertEqual(batch.mimetype, 'application/zip')
        self.assertEqual(design.call_args.kwargs['max_selected'], 2)
        self.assertEqual(design.call_args.kwargs['target_policy'],
                         {'mode': 'allowlist', 'allowed_ids': ['BoOther'], 'exclusion_max_mismatches': 1})
        batch.close()

    def test_changed_reference_and_running_task_fail_without_restarting(self):
        prepared = self.prepare()
        directory = self.root / 'split_results' / prepared['id']
        (directory / 'running.lock').touch()
        self.assertEqual(self.client.post(prepared['run_url'], json={}).status_code, 409)
        (directory / 'running.lock').unlink()
        self.reference.write_text('>BjuTarget\n' + SEQ + 'A\n')
        self.assertEqual(self.client.post(prepared['run_url'], json={}).status_code, 400)
        self.assertFalse((directory / 'running.lock').exists())

    def test_engine_failure_is_not_reported_as_empty_success(self):
        prepared = self.prepare()
        with patch('scripts.split_workflow.design_split', side_effect=RuntimeError('engine failed')):
            with self.assertLogs(self.app.logger, level='ERROR'):
                response = self.client.post(prepared['run_url'], json={})
        self.assertEqual(response.status_code, 202)
        status = self.client.get(response.json['status_url']).json
        self.assertEqual(status['status'], 'failed')
        self.assertIn('RuntimeError', status['error'])
        self.assertEqual(self.client.get('/split/results/' + prepared['id'],
                                         headers={'Accept': 'application/json'}).status_code, 409)
        self.assertIn('id="progress-panel"', self.client.get('/split/results/' + prepared['id']).text)
        with patch('scripts.split_workflow.design_split', side_effect=fake_design):
            response, result = self.run_job(prepared, {})
        self.assertEqual(response.status_code, 202)
        self.assertEqual(result['status'], 'completed')

    def failing_design(self, failing_task, error, stage='specificity'):
        def design(**kwargs):
            for name in ('reference', 'candidates', 'specificity', 'selection'):
                kwargs['progress'](name)
                if kwargs['name'] == failing_task and name == stage:
                    raise error() if callable(error) else error
            fake_design(**kwargs)
        return design

    def test_failed_batch_task_names_task_stage_and_cause_without_paths(self):
        data = 'task_name,BP_ID,gene_id\none,BP_0079,BjuTarget\ntwo,BP_0079,BoOther\n'
        prepared = self.prepare(input_mode='batch', batch_text=data)
        secret = str(self.root)
        error = ScanError(f'paired_scan failed (exit 1): cannot read {secret}/bol_bju_cdna.fasta', 1,
                          f'open {secret}/x failed')
        with patch('scripts.split_workflow.design_split', side_effect=self.failing_design('two', error)):
            with self.assertLogs(self.app.logger, level='ERROR') as logs:
                response = self.client.post(prepared['run_url'], json={})
        self.assertIn(secret, '\n'.join(logs.output))  # the full exception still reaches the log
        status = self.client.get(response.json['status_url'])
        self.assertEqual(status.json['status'], 'failed')
        self.assertEqual((status.json['error_code'], status.json['failed_task'], status.json['failed_stage']),
                         ('scan_failed', 'two', 'specificity'))
        self.assertEqual(status.json['failed_stage_label'], '双臂配对扫描与单臂比对')
        self.assertEqual(status.json['error'], FAILURES['scan_failed'])
        self.assertNotIn(secret, status.text)
        self.assertNotIn('Traceback', status.text)
        manifest = json.loads((self.root / 'split_results' / prepared['id'] / 'request.json').read_text())
        self.assertEqual(manifest['error_code'], 'scan_failed')
        english = self.client.get(response.json['status_url'] + '?lang=en').json
        self.assertTrue(english['error'].startswith('The paired-arm scan did not finish normally'))
        listing = self.client.get('/jobs', headers={'Accept': 'application/json'}).json['jobs']
        self.assertEqual(listing[0]['failed_task'], 'two')
        self.assertNotIn(secret, self.client.get('/jobs').text)
        self.assertNotIn(secret, self.client.get('/split/results/' + prepared['id']).text)
        # Rerunning clears the recorded failure.
        with patch('scripts.split_workflow.design_split', side_effect=fake_design):
            response, result = self.run_job(prepared, {})
        status = self.client.get(response.json['status_url']).json
        self.assertEqual((status['status'], status['error_code'], status['failed_task']), ('completed', None, None))

    def test_reference_removed_during_run_is_reported_as_missing(self):
        prepared = self.prepare()
        def vanish():
            self.reference.unlink()
            return FileNotFoundError(errno.ENOENT, 'No such file', str(self.reference))
        with patch('scripts.split_workflow.design_split', side_effect=self.failing_design('ABCG13', vanish, 'reference')):
            with self.assertLogs(self.app.logger, level='ERROR'):
                response = self.client.post(prepared['run_url'], json={})
        status = self.client.get(response.json['status_url']).json
        self.assertEqual((status['error_code'], status['failed_task'], status['failed_stage']),
                         ('reference_missing', 'ABCG13', 'reference'))
        self.assertEqual(status['failed_stage_label'], '准备参考库索引')


HAN = re.compile(r'[\u3400-\u9fff]')


class FailureReasonTests(unittest.TestCase):
    """Each failure class gets an actionable message without server paths or exception text."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.secret = os.path.join(self.tmp.name, 'private-server-folder')
        os.makedirs(self.secret)
        self.fasta = os.path.join(self.secret, 'ref.fasta')
        with open(self.fasta, 'w') as handle:
            handle.write('>BjuTarget\n' + SEQ + '\n')
        self.app = Flask(__name__)
        register_i18n(self.app)

    def raised(self, function, *args, **kwargs):
        try:
            function(*args, **kwargs)
        except Exception as error:  # noqa: BLE001 - the point is what the engine raises
            return error
        self.fail('%s did not raise' % function.__name__)

    def check(self, error, code, *expected, manifest=None):
        actual, text = describe_failure(error, manifest)
        self.assertEqual(actual, code, text)
        for part in expected:
            self.assertIn(part, text)
        for leak in (self.tmp.name, 'private-server-folder', 'Traceback', '/bin/', 'exit '):
            self.assertNotIn(leak, text)
        with self.app.test_request_context('/?lang=en'):
            english = message(text)
        self.assertIsNone(HAN.search(english), english)
        return text

    def write(self, body):
        with open(self.fasta, 'w') as handle:
            handle.write(body)

    def test_rejected_reference_fasta_names_sequence_and_reason(self):
        self.write('>prot_1 description\nMEEPQSDLACGT\n')
        text = self.check(self.raised(create_blast_db, self.fasta), 'reference_fasta', 'prot_1', '非核酸字符 ELPQ')
        self.assertTrue(text.startswith('参考 FASTA 未通过检查：序列 prot_1 含有非核酸字符'))
        self.write('>empty\n>next\nACGT\n')
        self.check(self.raised(create_blast_db, self.fasta), 'reference_fasta', '序列 empty 没有碱基')
        self.write('')
        self.check(self.raised(create_blast_db, self.fasta), 'reference_fasta', '没有任何序列')
        self.write('>dup\nACGT\n>dup\nACGT\n')
        self.check(self.raised(load_reference, self.fasta), 'reference_fasta', '序列编号 dup 重复')
        self.write('ACGT\n>late\nACGT\n')
        self.check(self.raised(load_reference, self.fasta), 'reference_fasta', '标题行之前')
        self.write('>a\nACGT\n> \nACGT\n')
        self.check(self.raised(load_reference, self.fasta), 'reference_fasta', '没有编号的标题行')
        self.write('>' + 'x' * 300 + '\nACGTJ\n')
        text = self.check(self.raised(create_blast_db, self.fasta), 'reference_fasta')
        self.assertLess(len(text), 200)  # long IDs are shortened

    def test_missing_or_changed_reference(self):
        missing = os.path.join(self.secret, 'gone.fasta')
        self.check(self.raised(create_blast_db, missing), 'reference_missing', '已被删除或清理')
        self.check(self.raised(load_reference, missing), 'reference_missing')
        self.check(self.raised(scan_pairs, [('ACGT' * 4, 'ACGT' * 4, 0)], missing), 'reference_missing')
        # Whatever broke afterwards, a vanished or edited reference is reported as the cause.
        manifest = {'reference_path': missing, 'reference': {'fingerprint': {'mtime_ns': 1, 'size': 1}}}
        self.check(KeyError('x'), 'reference_missing', manifest=manifest)
        manifest = dict(manifest, reference_path=self.fasta)
        self.check(ScanError('boom'), 'reference_changed', '被替换或修改', manifest=manifest)
        stat = os.stat(self.fasta)
        manifest['reference']['fingerprint'] = {'mtime_ns': stat.st_mtime_ns, 'size': stat.st_size}
        self.check(ScanError('boom'), 'scan_failed', manifest=manifest)

    def test_blast_missing_or_failing(self):
        with patch.dict(os.environ, {'HCR_BLAST_BIN': os.path.join(self.secret, 'nothing')}):
            self.check(self.raised(analyze_probe_sequences, ['ACGTACGTAC'], self.fasta),
                       'blast_unavailable', 'blastn', 'HCR_BLAST_BIN')
            self.write('>BjuTarget\nACGTACGT\n')
            self.check(self.raised(create_blast_db, self.fasta), 'blast_unavailable', 'makeblastdb')
        failure = subprocess.CalledProcessError(2, ['blastn'], stderr='BLAST Database error: ' + self.secret)
        with patch('src.common.paired_specificity.blast_executable', return_value='blastn'), \
                patch('src.common.paired_specificity.subprocess.run', side_effect=failure):
            error = self.raised(analyze_probe_sequences, ['ACGTACGTAC'], self.fasta)
        self.assertIsInstance(error, subprocess.CalledProcessError)  # existing handlers still match
        self.check(error, 'blast_failed', 'blastn')
        def makeblastdb(**kwargs):
            def run():
                raise subprocess.CalledProcessError(1, ['makeblastdb'], stderr='Error: ' + self.secret)
            return run
        with patch('scripts.utils.NcbimakeblastdbCommandline', side_effect=makeblastdb), \
                patch('scripts.utils.blast_executable', return_value='makeblastdb'):
            self.check(self.raised(create_blast_db, self.fasta), 'blast_failed', 'makeblastdb')
        self.check(ToolFailedError('makeblastdb', 1, stderr='No space left on device'), 'disk_full')

    def test_pair_scanner_unavailable_failed_or_timed_out(self):
        pair = ('ACGTACGTACGTACGT', 'ACGTACGTACGTACGT', 0)
        with patch.dict(os.environ, {'HCR_SCAN_BIN': os.path.join(self.secret, 'absent')}):
            error = self.raised(scan_pairs, [pair], self.fasta)
        self.assertIsInstance(error, FileNotFoundError)
        self.check(error, 'scanner_unavailable', 'HCR_SCAN_BIN')
        environment = {k: v for k, v in os.environ.items() if k not in ('HCR_SCAN_BIN', 'CXX')}
        with patch.dict(os.environ, environment, clear=True), \
                patch.object(paired_scan, 'PROJECT_ROOT', Path(self.secret)), \
                patch('src.common.paired_scan.shutil.which', return_value=None):
            self.check(self.raised(scan_pairs, [pair], self.fasta), 'scanner_unavailable', 'C++')
        with patch.object(paired_scan, 'scanner_executable', return_value='/bin/sh'):
            self.check(self.raised(scan_pairs, [pair], self.fasta), 'scan_failed', '可直接重新运行')
        self.check(ScanTimeout('paired_scan exceeded 5 s on ref.fasta'), 'scan_timeout', '时间限制')
        self.check(ScanError('paired_scan failed (exit -9): ', -9, ''), 'out_of_memory', '内存不足')
        self.check(ScanError('failed', 1, "terminate called after throwing 'std::bad_alloc'"), 'out_of_memory')

    def test_disk_permission_unexpected_and_input_errors(self):
        self.check(OSError(errno.ENOSPC, 'No space left on device', self.secret), 'disk_full', '磁盘空间不足')
        wrapped = RuntimeError('write failed')
        wrapped.__cause__ = OSError(errno.EDQUOT, 'Disc quota exceeded', self.secret)
        self.check(wrapped, 'disk_full')
        self.check(PermissionError(errno.EACCES, 'Permission denied', self.secret), 'permission', '权限')
        self.check(OSError(errno.EROFS, 'Read-only file system', self.secret), 'permission')
        self.check(MemoryError(), 'out_of_memory')
        text = self.check(KeyError(self.secret), 'unexpected', 'KeyError')
        self.assertEqual(text, '计算未完成（KeyError）。请检查参考文件与比对环境后重新运行。')
        self.assertEqual(describe_failure(InputError('请选择参考库。')), ('input', '请选择参考库。'))


if __name__ == '__main__':
    unittest.main()
