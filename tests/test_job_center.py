import json
from datetime import timedelta
from pathlib import Path
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from flask import Flask
from scripts.split_web import register_split_routes
from scripts.split_workflow import DEFAULTS, SplitWorkflow, read_synthesis
from scripts.admin_web import register_admin
from scripts.storage_maintenance import utcnow


SEQ = 'ACGTTGCAAGGCTTACCGATGCATTAGCCGGATCCAATTGGCATCGATCGGATTACCGTAGCTAGGCTAACG'
REPO = Path(__file__).resolve().parents[1]
release = threading.Event()


def fake_design(**kwargs):
    """Report each stage, then write the files the workflow reads back."""
    for stage in ('candidates', 'specificity', 'selection'):
        kwargs['progress'](stage)
    if kwargs['name'] == 'slow':
        release.wait(5)
    output = Path(kwargs['output_dir'])
    probe = {'id': 1, 'keep': True, 'reason_code': 'selected', 'evaluation': {},
             'specificity': {'status': 'completed', 'on_target_exact': True, 'counts': {'cross_species': [0, 0, 0, 0]},
                             'paired_matches': []},
             'left_probe': {'sequence': 'GCAACGT', 'position': {'start': 1, 'end': 7}, 'tm': 50, 'gc_content': 50},
             'right_probe': {'sequence': 'AAGCCTT', 'position': {'start': 8, 'end': 14}, 'tm': 50, 'gc_content': 50}}
    (output / (kwargs['name'] + '_filtered.json')).write_text(json.dumps({'probe_sets': [probe]}))
    (output / 'dual_probe_dual.csv').write_text(
        'Probe ID,Sequence\n{0}-1-L,GCAACGTANCHOR\n{0}-1-R-{1},BRIDGEAAGCCTT\n'.format(kwargs['name'], kwargs['bridge_probe_id']))


class JobCenterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.app = Flask(__name__, template_folder=str(REPO / 'templates'), static_folder=str(REPO / 'static'))
        self.app.config.update(TESTING=True, UPLOAD_FOLDER=str(self.root))
        for endpoint in ['index', 'triplet', 'triplet_batch', 'bridge', 'bridge_list', 'bridge_query', 'bridge_sequence',
                         'bridge_analyze', 'bridge_validate', 'formatting', 'subgenomes', 'references', 'probe_recommend', 'filter', 'help', 'orders']:
            self.app.add_url_rule('/' + endpoint, endpoint, lambda: 'placeholder')
        self.service_args = ({'BP_0079': 'A' * 19}, {'BjuTarget': SEQ}, {}, {})
        register_split_routes(self.app, *self.service_args, lambda: [])
        self.client = self.app.test_client()
        release.clear()

    def tearDown(self):
        release.set()
        self.app.extensions['split_jobs'].executor.shutdown(wait=True)
        self.tmp.cleanup()

    def prepare(self, name='ABCG13'):
        response = self.client.post('/split/prepare', data={
            'input_mode': 'single', 'name': name, 'geneID': 'BjuTarget', 'bp_id': 'BP_0079',
            'sequenceType': 'cdna', 'reference_mode': 'none'})
        self.assertEqual(response.status_code, 200, response.json)
        return response.json

    def wait(self, status_url, until=('completed', 'failed'), stage=None):
        for _ in range(100):
            status = self.client.get(status_url).json
            if status['status'] in until and (stage is None or (status['progress'] or {}).get('stage') == stage):
                return status
            time.sleep(.05)
        self.fail('job did not finish: %s' % status)

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_jobs_run_in_background_one_at_a_time(self, design):
        slow, quick = self.prepare('slow'), self.prepare('quick')
        first = self.client.post(slow['run_url'], json=DEFAULTS)
        self.assertEqual(first.status_code, 202)
        second = self.client.post(quick['run_url'], json=DEFAULTS)
        self.assertEqual(second.status_code, 202)
        running = self.wait(first.json['status_url'], until=('running',), stage='selection')
        self.assertEqual(running['progress']['label'], '排序与筛选')
        self.assertEqual(running['progress']['task_name'], 'slow')
        waiting = self.client.get(second.json['status_url']).json
        self.assertEqual((waiting['status'], waiting['queue_position']), ('queued', 1))
        self.assertEqual(self.client.post(slow['run_url'], json=DEFAULTS).status_code, 409)  # no duplicate runs
        release.set()
        self.assertEqual(self.wait(first.json['status_url'])['status'], 'completed')
        done = self.wait(second.json['status_url'])
        self.assertEqual(done['summary']['selected'], 1)
        result = self.client.get(done['result_url'], headers={'Accept': 'application/json'}).json
        task = result['tasks'][0]
        self.assertEqual(task['sequence'], SEQ)
        self.assertEqual(task['probe_sets'][0]['synthesis']['R'],
                         {'id': 'quick-1-R-BP_0079', 'sequence': 'BRIDGEAAGCCTT'})
        page = self.client.get(done['result_url'])
        self.assertEqual(page.status_code, 200)
        self.assertIn('id="track"', page.text)

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_job_list_hides_unsubmitted_and_shows_newest_first(self, design):
        self.prepare('never_run')
        older, newer = self.prepare('older'), self.prepare('newer')
        for job in (older, newer):
            self.wait(self.client.post(job['run_url'], json=DEFAULTS).json['status_url'])
        listing = self.client.get('/jobs', headers={'Accept': 'application/json'}).json['jobs']
        self.assertEqual([job['tasks'][0]['name'] for job in listing], ['newer', 'older'])
        self.assertTrue(all(job['status'] == 'completed' for job in listing))
        self.assertIn('任务记录', self.client.get('/jobs').text)

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_restart_marks_unfinished_jobs_interrupted_and_they_can_rerun(self, design):
        prepared = self.prepare()
        service = SplitWorkflow(self.root, *self.service_args)
        service.queue(prepared['id'], DEFAULTS)  # queued, but the process "dies" before running it
        (self.root / 'split_results' / prepared['id'] / 'running.lock').touch()
        self.assertEqual(service.recover(), [prepared['id']])
        status = self.client.get('/split/jobs/%s/status' % prepared['id']).json
        self.assertEqual(status['status'], 'interrupted')
        self.assertEqual(status['parameters']['min_gap'], DEFAULTS['min_gap'])
        rerun = self.client.post(prepared['run_url'], json=status['parameters'])
        self.assertEqual(rerun.status_code, 202)
        self.assertEqual(self.wait(rerun.json['status_url'])['status'], 'completed')

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_results_saved_before_this_version_gain_sequences_and_synthesis(self, design):
        prepared = self.prepare()
        self.wait(self.client.post(prepared['run_url'], json=DEFAULTS).json['status_url'])
        path = self.root / 'split_results' / prepared['id'] / 'result.json'
        old = json.loads(path.read_text())
        for task in old['tasks']:
            del task['sequence']
            for probe in task['probe_sets']:
                del probe['synthesis']
        path.write_text(json.dumps(old))
        task = SplitWorkflow(self.root, *self.service_args).result(prepared['id'])['tasks'][0]
        self.assertEqual(task['sequence'], SEQ)
        self.assertEqual(task['probe_sets'][0]['synthesis']['L']['id'], 'ABCG13-1-L')

    def test_read_synthesis_ignores_other_tasks(self):
        path = self.root / 'dual.csv'
        path.write_text('Probe ID,Sequence\nA-1-L,AAA\nA-1-R-BP_1,CCC\nAB-2-L,GGG\n')
        self.assertEqual(read_synthesis(path, 'A'), {1: {'L': {'id': 'A-1-L', 'sequence': 'AAA'},
                                                         'R': {'id': 'A-1-R-BP_1', 'sequence': 'CCC'}}})

    @patch('scripts.split_workflow.design_split', side_effect=fake_design)
    def test_admin_cleanup_protects_running_jobs_and_references_then_removes_completed_results(self, design):
        with tempfile.TemporaryDirectory() as config:
            self.app.config['ADMIN_CREDENTIALS_FILE'] = str(Path(config) / 'admin.json')
            register_admin(self.app)
            initialized = self.app.test_cli_runner().invoke(
                args=['init-admin', '--username', 'manager'], input='integration-password\nintegration-password\n')
            self.assertEqual(initialized.exit_code, 0, initialized.output)
            self.client.get('/admin/login')
            with self.client.session_transaction() as session:
                csrf = session['admin_csrf']
            self.client.post('/admin/login', data=dict(username='manager', password='integration-password', csrf_token=csrf))
            self.client.get('/admin')
            with self.client.session_transaction() as session:
                csrf = session['admin_csrf']
            headers = {'X-CSRF-Token': csrf}
            reference = self.root / 'reference.fasta'
            reference.write_text('>BjuTarget\n' + SEQ + '\n')
            prepared = self.client.post('/split/prepare', data=dict(
                input_mode='single', name='slow', geneID='BjuTarget', bp_id='BP_0079',
                sequenceType='reference', reference_mode='existing', existing_ref_genome=reference.name)).json
            first = self.client.post(prepared['run_url'], json=DEFAULTS)
            self.wait(first.json['status_url'], until=('running',), stage='selection')
            key = 'split_results/' + prepared['id']
            for entry in (key, reference.name):
                self.assertEqual(self.client.post('/admin/delete', json={'keys': [entry]}, headers=headers).status_code, 409)
            report = self.app.extensions['storage_maintenance'].cleanup(now=utcnow() + timedelta(days=100))
            self.assertEqual(report['skipped'], [key])
            release.set()
            done = self.wait(first.json['status_url'])
            self.assertEqual(done['status'], 'completed')
            removed = self.client.post('/admin/delete', json={'keys': [key]}, headers=headers)
            self.assertEqual(removed.status_code, 200)
            self.assertTrue(reference.exists())
            for route in (done['status_url'], done['result_url'], done['download_url']):
                self.assertEqual(self.client.get(route).status_code, 404)
            self.assertEqual(self.client.get('/jobs', headers={'Accept': 'application/json'}).json['jobs'], [])


class OrderCheckTests(unittest.TestCase):
    def test_bridge_conflicts_and_unknown_ids(self):
        import app as web
        client = web.app.test_client()
        response = client.post('/orders/check', json={'bp_ids': ['BP_0079', 'bp_0745', 'BP_9999'], 'kmer': 6})
        self.assertEqual(response.status_code, 200)
        data = response.json
        self.assertEqual(data['checked'], ['BP_0079', 'BP_0745'])
        self.assertEqual(data['unknown'], ['BP_9999'])
        self.assertEqual([(c['a'], c['b']) for c in data['conflicts']], [('BP_0079', 'BP_0745')])
        self.assertEqual(client.post('/orders/check', json={'bp_ids': ['BP_0079'], 'kmer': 'x'}).status_code, 400)
        self.assertEqual(client.get('/orders').status_code, 200)


if __name__ == '__main__':
    unittest.main()
