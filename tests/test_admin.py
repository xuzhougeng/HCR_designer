import json
import os
from pathlib import Path
import tempfile
import threading
from datetime import datetime, timezone
import unittest
from unittest.mock import patch

from flask import Flask
from werkzeug.security import check_password_hash

from scripts.admin_web import register_admin
from scripts.storage_maintenance import CleanupScheduler, StorageMaintenance, three_months_ago


REPO = Path(__file__).resolve().parents[1]
NOW = datetime(2026, 10, 15, 12, tzinfo=timezone.utc)
OLD = datetime(2026, 7, 14, 12, tzinfo=timezone.utc)
PASSWORD = 'test-admin-password-123'


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.storage = StorageMaintenance(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def job(self, identifier='a' * 32, status='completed', when=OLD, **extra):
        directory = self.root / 'split_results' / identifier
        directory.mkdir(parents=True)
        manifest = dict(id=identifier, status=status, created_at=when.isoformat(), **extra)
        if status == 'completed':
            manifest['completed_at'] = when.isoformat()
        (directory / 'request.json').write_text(json.dumps(manifest))
        (directory / 'results.zip').write_text('probes')
        return directory

    def age(self, path):
        for p in sorted(path.rglob('*'), reverse=True) if path.is_dir() else []:
            os.utime(p, (OLD.timestamp(), OLD.timestamp()))
        os.utime(path, (OLD.timestamp(), OLD.timestamp()))

    def test_retention_calendar_months_and_exact_boundary(self):
        self.assertEqual(three_months_ago(datetime(2026, 5, 31, tzinfo=timezone.utc)).day, 28)
        self.assertEqual(three_months_ago(datetime(2024, 5, 31, tzinfo=timezone.utc)).day, 29)
        self.assertEqual(three_months_ago(datetime(2026, 1, 31, tzinfo=timezone.utc)),
                         datetime(2025, 10, 31, tzinfo=timezone.utc))
        expired = self.job()
        boundary = self.job('b' * 32, when=three_months_ago(NOW))
        recent = self.job('c' * 32, when=NOW)
        preview = self.storage.cleanup(dry_run=True, now=NOW)
        self.assertEqual(preview['deleted'], ['split_results/' + 'a' * 32])
        self.assertTrue(expired.exists())
        report = self.storage.cleanup(now=NOW)
        self.assertEqual(report['deleted'], preview['deleted'])
        self.assertFalse(expired.exists())
        self.assertTrue(boundary.exists() and recent.exists())
        self.assertEqual(self.storage.cleanup(now=NOW)['deleted'], [])

    def test_prepared_results_expire_but_active_and_damaged_jobs_are_protected(self):
        prepared = self.job(status='prepared')
        queued = self.job('b' * 32, status='queued')
        running = self.job('c' * 32, status='running')
        locked = self.job('d' * 32)
        (locked / 'running.lock').touch()
        damaged = self.job('e' * 32)
        (damaged / 'request.json').write_text('broken')
        self.age(damaged)
        report = self.storage.cleanup(now=NOW)
        self.assertFalse(prepared.exists())
        self.assertTrue(all(p.exists() for p in (queued, running, locked, damaged)))
        self.assertEqual(len(report['skipped']), 4)
        for p in (queued, running, locked, damaged):
            with self.assertRaises(ValueError):
                self.storage.delete([str(p.relative_to(self.root))])

    def test_triplet_folder_zip_and_old_recommendations_are_cleaned_together(self):
        folder = self.root / 'triplet_results' / ('a' * 32)
        folder.mkdir(parents=True)
        (folder / 'probes.csv').write_text('probe')
        archive = folder.parent / ('triplet_batch_results_' + folder.name + '.zip')
        archive.write_text('zip')
        orphan = folder.parent / 'triplet_results_old.zip'
        orphan.write_text('zip')
        recommendation = self.root / ('probe_recommend_' + 'f' * 32 + '.csv')
        recommendation.write_text('probes')
        for path in (folder, archive, orphan, recommendation):
            self.age(path)
        report = self.storage.cleanup(now=NOW)
        self.assertEqual(len(report['deleted']), 3)
        self.assertFalse(any(p.exists() for p in (folder, archive, orphan, recommendation)))

    def test_reference_and_indexes_grouped_manual_only(self):
        names = ['ref.fa', 'ref.fa.gz', 'ref.fa.nhr', 'ref.fa.00.nsq', 'ref.fa.njs']
        for name in names + ['ref.fa.notes', 'uploaded-probes.csv']:
            (self.root / name).write_text('data')
            self.age(self.root / name)
        self.assertEqual(self.storage.cleanup(now=NOW)['deleted'], [])
        reference = next(e for e in self.storage.inventory(NOW) if e['key'] == 'ref.fa')
        self.assertEqual(reference['file_count'], 5)
        self.assertEqual(reference['size'], 20)
        report = self.storage.delete(['ref.fa'])
        self.assertEqual(report['bytes'], 20)
        self.assertTrue(all(not (self.root / name).exists() for name in names))
        self.assertTrue((self.root / 'ref.fa.notes').exists())
        self.assertTrue((self.root / 'uploaded-probes.csv').exists())

    def test_uppercase_compressed_references_and_invalid_status_records_are_handled(self):
        (self.root / 'reference.FASTA.GZ').write_text('compressed')
        directory = self.job(status=['broken'])
        entry = next(e for e in self.storage.inventory(NOW) if e['key'] == 'reference.FASTA')
        self.assertEqual(entry['size'], 10)
        self.assertTrue(entry['blocked'])  # an unreadable job state cannot prove this reference is unused
        report = self.storage.cleanup(now=NOW)
        self.assertEqual(report['skipped'], ['split_results/' + directory.name])
        self.assertTrue(directory.exists())

    def test_active_shared_reference_and_legacy_requests_block_deletion(self):
        reference = self.root / 'ref.fasta'
        reference.write_text('reference')
        self.job(status='queued', reference_path=str(reference))
        with self.assertRaises(ValueError):
            self.storage.delete(['ref.fasta'])
        self.storage.legacy_active = 1
        (self.root / 'uploaded.csv').touch()
        with self.assertRaises(ValueError):
            self.storage.delete(['uploaded.csv'])

    def test_path_traversal_symlinks_and_invalid_keys_cannot_delete_outside_storage(self):
        outside = self.root.parent / (self.root.name + '-outside')
        outside.mkdir()
        try:
            (outside / 'keep').write_text('keep')
            (self.root / 'linked').symlink_to(outside, target_is_directory=True)
            directory = self.job()
            (directory / 'external').symlink_to(outside, target_is_directory=True)
            for keys in (['../' + outside.name], [str(outside)], ['linked'], [None], {}, []):
                with self.assertRaises(ValueError):
                    self.storage.delete(keys)
            self.storage.cleanup(now=NOW)
            self.assertEqual((outside / 'keep').read_text(), 'keep')
            self.assertTrue(directory.exists())
        finally:
            (outside / 'keep').unlink()
            outside.rmdir()

    def test_failed_deletion_reports_error_and_continues(self):
        first = self.job()
        second = self.job('b' * 32)
        import shutil
        original = shutil.rmtree

        def remove(path, *args, **kwargs):
            if path == first:
                raise PermissionError('denied')
            return original(path, *args, **kwargs)

        with patch('scripts.storage_maintenance.shutil.rmtree', side_effect=remove):
            report = self.storage.cleanup(now=NOW)
        self.assertEqual(report['errors'], ['split_results/' + first.name])
        self.assertEqual(report['deleted'], ['split_results/' + second.name])
        self.assertTrue(first.exists())
        self.assertFalse(second.exists())

    def test_scheduler_runs_repeatedly_and_recovers_after_errors(self):
        completed = threading.Event()
        calls = []

        def cleanup():
            calls.append(1)
            if len(calls) == 1:
                raise OSError('temporary failure')
            completed.set()
            return dict(deleted=[], errors=[])

        with patch.object(self.storage, 'cleanup', side_effect=cleanup):
            from unittest.mock import Mock
            scheduler = CleanupScheduler(self.storage, Mock(), interval=.01)
            scheduler.start()
            self.assertTrue(completed.wait(2))
            scheduler.stop.set()
            scheduler.thread.join(2)
        self.assertGreaterEqual(len(calls), 2)


class AdminTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.app = Flask(__name__, template_folder=str(REPO / 'templates'), static_folder=str(REPO / 'static'))
        self.app.config.update(TESTING=True, UPLOAD_FOLDER=str(self.root / 'uploads'),
                               ADMIN_CREDENTIALS_FILE=str(self.root / 'instance' / 'admin.json'))
        (self.root / 'uploads').mkdir()
        for name in ('index', 'split', 'split_batch', 'triplet', 'triplet_batch', 'jobs', 'orders', 'bridge_list',
                     'bridge_query', 'bridge', 'bridge_validate', 'formatting', 'subgenomes', 'references', 'probe_recommend', 'filter', 'help'):
            self.app.add_url_rule('/' + name, name, lambda: 'placeholder')
        register_admin(self.app)
        self.client = self.app.test_client()
        self.runner = self.app.test_cli_runner()

    def tearDown(self):
        self.tmp.cleanup()

    def initialize(self, reset=False):
        result = self.runner.invoke(args=['init-admin', '--username', 'manager'] + (['--reset'] if reset else []),
                                    input=PASSWORD + '\n' + PASSWORD + '\n')
        self.assertEqual(result.exit_code, 0, result.output)

    def token(self, page='/admin/login'):
        self.client.get(page)
        with self.client.session_transaction() as session:
            return session['admin_csrf']

    def login(self):
        response = self.client.post('/admin/login', data=dict(username='manager', password=PASSWORD,
                                                            csrf_token=self.token()))
        self.assertEqual(response.status_code, 302)
        return self.token('/admin')

    def test_initialization_hashes_password_restricts_permissions_and_refuses_overwrite(self):
        self.assertEqual(self.client.get('/admin/login').status_code, 503)
        self.initialize()
        path = Path(self.app.config['ADMIN_CREDENTIALS_FILE'])
        data = json.loads(path.read_text())
        self.assertNotIn(PASSWORD, path.read_text())
        self.assertTrue(check_password_hash(data['password_hash'], PASSWORD))
        self.assertEqual(path.stat().st_mode & 0o777, 0o600)
        result = self.runner.invoke(args=['init-admin', '--username', 'other'])
        self.assertNotEqual(result.exit_code, 0)
        self.assertEqual(json.loads(path.read_text())['username'], 'manager')

    def test_uninitialized_and_unauthenticated_mutations_are_rejected(self):
        self.assertEqual(self.client.get('/admin').status_code, 302)
        for route in ('/admin/delete', '/admin/cleanup', '/admin/logout'):
            self.assertEqual(self.client.post(route, json={}).status_code, 401)
        self.initialize()
        self.assertEqual(self.client.get('/admin', headers={'Accept': 'application/json'}).status_code, 401)

    def test_authentication_csrf_logout_and_session_reset(self):
        self.initialize()
        self.assertEqual(self.client.post('/admin/login', data={'username': 'manager', 'password': PASSWORD}).status_code, 400)
        token = self.login()
        response = self.client.get('/admin')
        self.assertIn('空间管理', response.text)
        self.assertEqual(response.headers['Cache-Control'], 'no-store')
        self.assertEqual(self.client.post('/admin/cleanup', json={}).status_code, 400)
        self.assertEqual(self.client.post('/admin/logout', data={'csrf_token': token}).status_code, 302)
        self.assertEqual(self.client.get('/admin').status_code, 302)
        self.login()
        self.initialize(reset=True)
        self.assertEqual(self.client.get('/admin').status_code, 302)

    def test_failed_login_is_generic_and_rate_limited(self):
        self.initialize()
        token = self.token()
        for index in range(5):
            response = self.client.post('/admin/login', data=dict(username='unknown', password='bad', csrf_token=token))
            self.assertEqual(response.status_code, 401)
            self.assertIn('账号或密码错误', response.text)
        self.assertEqual(self.client.post('/admin/login', data=dict(username='manager', password=PASSWORD,
                                                                  csrf_token=token)).status_code, 429)

    def test_listing_deletion_and_cleanup_only_touch_selected_or_expired_uploads(self):
        self.initialize()
        token = self.login()
        uploads = self.root / 'uploads'
        (uploads / 'probes.csv').write_text('probe')
        directory = uploads / 'split_results' / ('a' * 32)
        directory.mkdir(parents=True)
        (directory / 'request.json').write_text(json.dumps(dict(status='completed', completed_at=OLD.isoformat())))
        payload = self.client.get('/admin', headers={'Accept': 'application/json'}).json
        self.assertEqual(len(payload['entries']), 2)
        self.assertTrue(all('_paths' not in item for item in payload['entries']))
        response = self.client.post('/admin/delete', json={'keys': ['../../resources']}, headers={'X-CSRF-Token': token})
        self.assertEqual(response.status_code, 409)
        response = self.client.post('/admin/delete', json={'keys': ['probes.csv']}, headers={'X-CSRF-Token': token})
        self.assertEqual(response.json['deleted'], ['probes.csv'])
        self.assertFalse((uploads / 'probes.csv').exists())
        response = self.client.post('/admin/cleanup', json={}, headers={'X-CSRF-Token': token})
        self.assertEqual(response.json['deleted'], ['split_results/' + 'a' * 32])
        self.assertFalse(directory.exists())

    def test_cli_dry_run_preserves_files_without_needing_admin_login(self):
        directory = self.root / 'uploads' / 'triplet_results'
        directory.mkdir()
        archive = directory / 'old.zip'
        archive.write_text('old results')
        os.utime(archive, (OLD.timestamp(), OLD.timestamp()))
        preview = self.runner.invoke(args=['cleanup-results', '--dry-run'])
        self.assertEqual(preview.exit_code, 0, preview.output)
        self.assertEqual(json.loads(preview.output)['deleted'], ['triplet_results/old.zip'])
        self.assertTrue(archive.exists())
        actual = self.runner.invoke(args=['cleanup-results'])
        self.assertEqual(actual.exit_code, 0, actual.output)
        self.assertFalse(archive.exists())

    def test_enabled_scheduler_never_starts_from_cli_and_starts_only_once_on_serving(self):
        app = Flask('scheduler-test')
        app.config.update(UPLOAD_FOLDER=str(self.root / 'uploads'),
                          ADMIN_CREDENTIALS_FILE=str(self.root / 'scheduler-admin.json'))
        with patch.dict(os.environ, HCR_CLEANUP_ENABLED='1'):
            register_admin(app)
        scheduler = app.extensions['storage_cleanup_scheduler']
        with patch.object(scheduler, 'start') as start:
            self.assertEqual(app.test_cli_runner().invoke(args=['cleanup-results', '--dry-run']).exit_code, 0)
            start.assert_not_called()
            client = app.test_client()
            client.get('/nonexistent')
            start.assert_called_once()
            # A real start assigns ident; reflect that without starting a destructive worker.
            with patch.object(type(scheduler.thread), 'ident', property(lambda self: 1)):
                client.get('/nonexistent')
            start.assert_called_once()


if __name__ == '__main__':
    unittest.main()
