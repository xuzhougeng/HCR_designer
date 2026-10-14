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
    (output / (kwargs['name'] + '_filtered.json')).write_text(json.dumps({'probe_sets': [probe]}))
    (output / 'dual_probe.csv').write_text('Probe ID,Sequence\n')


class GuidedWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.reference = self.root / 'bol_bju_cdna.fasta'
        self.reference.write_text('>BjuTarget\n' + SEQ + '\n>BoOther\n' + SEQ + '\n')
        repository = Path(__file__).resolve().parents[1]
        self.app = Flask(__name__, template_folder=str(repository / 'templates'), static_folder=str(repository / 'static'))
        self.app.config.update(TESTING=True, UPLOAD_FOLDER=str(self.root))
        for endpoint in ['index', 'triplet', 'bridge_list', 'help']:
            self.app.add_url_rule('/' + endpoint, endpoint, lambda: 'placeholder')
        register_split_routes(self.app, {'BP_0079': 'A' * 19}, {'BjuTarget': 'T' * 80}, {}, {},
                              lambda: [{'filename': self.reference.name, 'display_name': 'Bo/Bju'}])
        self.client = self.app.test_client()
        self.input = {'input_mode': 'single', 'name': 'ABCG13', 'geneID': 'bjutarget', 'bp_id': 'BP_0079',
                      'sequenceType': 'reference', 'reference_mode': 'existing',
                      'existing_ref_genome': self.reference.name}

    def tearDown(self):
        self.tmp.cleanup()

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
        self.assertIn('task_name,BP_ID,gene_id,max_selected', template.data.decode('utf-8-sig'))
        template.close()

    def test_bojbu_target_comes_from_selected_reference_case_insensitively(self):
        prepared = self.prepare()
        task = prepared['tasks'][0]
        self.assertEqual(task['gene_id'], 'BjuTarget')
        self.assertNotIn('sequence', task)
        manifest = json.loads((self.root / 'split_results' / prepared['id'] / 'request.json').read_text())
        self.assertEqual(manifest['tasks'][0]['sequence'], SEQ)  # not the conflicting internal cDNA
        self.assertEqual(prepared['reference']['species_counts'], {'Bju': 1, 'Bo': 1})

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
        response = self.client.post(prepared['run_url'], json=DEFAULTS)
        self.assertEqual(response.status_code, 200, response.json)
        result = response.json
        self.assertEqual(result['tasks'][0]['cross_exact'], 1)
        self.assertEqual(result['tasks'][0]['selected'], 0)
        self.assertFalse(result['search_exhaustive'])
        self.assertNotIn('blast_results', result['tasks'][0]['probe_sets'][0]['left_probe'])
        self.assertEqual(self.client.get(result['result_url']).status_code, 200)
        archive = self.client.get(result['download_url'])
        with ZipFile(io.BytesIO(archive.data)) as z:
            self.assertIn('run_summary.json', z.namelist())
            full = json.loads(z.read('ABCG13_filtered.json'))
            self.assertIn('blast_results', full['probe_sets'][0]['left_probe'])
        archive.close()
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
        response.close()
        batch = self.client.post('/split/batch', data=dict(self.input,
            batch_file=(io.BytesIO(b'one,BP_0079,BjuTarget,2\n'), 'batch.csv')))
        self.assertEqual(batch.mimetype, 'application/zip')
        self.assertEqual(design.call_args.kwargs['max_selected'], 2)
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
        self.assertEqual(response.status_code, 500)
        self.assertEqual(self.client.get('/split/results/' + prepared['id']).status_code, 409)
        with patch('scripts.split_workflow.design_split', side_effect=fake_design):
            self.assertEqual(self.client.post(prepared['run_url'], json={}).status_code, 200)


if __name__ == '__main__':
    unittest.main()
