import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

from flask import Flask
from scripts.split_web import register_split_routes
from scripts.split_workflow import InputError, parse_parameters, DEFAULTS
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


if __name__ == '__main__':
    unittest.main()
