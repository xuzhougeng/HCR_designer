from datetime import datetime, timedelta, timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flask import Flask
from werkzeug.datastructures import FileStorage

from scripts.admin_web import register_admin
from scripts.reference_library import ReferenceLibrary
from scripts.reference_web import register_reference_routes
from scripts.split_workflow import SplitWorkflow
from scripts.storage_maintenance import REFERENCE_RETENTION_DAYS, StorageMaintenance


FASTA = b'>BjuVA01G1 description\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n'


def uploaded(name='reference.fa', data=FASTA):
    return FileStorage(stream=io.BytesIO(data), filename=name)


class ReferenceLibraryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.uploads = self.root / 'uploads'
        self.uploads.mkdir()
        self.library = ReferenceLibrary(self.uploads, self.root / 'instance' / 'references.json')

    def tearDown(self):
        self.tmp.cleanup()

    def test_upload_hash_configuration_and_restart(self):
        result = self.library.upload(uploaded(), 'Bju reference', True)
        entry = result['reference']
        self.assertFalse(result['duplicate'])
        self.assertEqual(entry['sha256'], hashlib.sha256(FASTA).hexdigest())
        self.assertEqual(entry['sequence_count'], 1)
        self.assertTrue(entry['polyploid'])
        restarted = ReferenceLibrary(self.uploads, self.library.catalog_path)
        self.assertEqual(restarted.list(), self.library.list())
        updated = restarted.configure(entry['filename'], 'New name', False)
        self.assertEqual(updated['display_name'], 'New name')
        self.assertFalse(updated['polyploid'])
        self.assertEqual((self.uploads / entry['filename']).read_bytes(), FASTA)

    def test_renamed_and_recompressed_duplicates_reuse_existing_file_and_settings(self):
        first = self.library.upload(uploaded(), 'Original', True)['reference']
        for name, data in [('renamed.fasta', FASTA), ('compressed.fa.gz', gzip.compress(FASTA, mtime=1)),
                           ('compressed-again.fna.GZ', gzip.compress(FASTA, mtime=99))]:
            duplicate = self.library.upload(uploaded(name, data), 'Duplicate name', False)
            self.assertTrue(duplicate['duplicate'])
            # A duplicate upload counts as a use, so only the expiry fields move.
            usage = ('last_used_at', 'expires_at', 'days_left', 'expiring_soon')
            self.assertEqual({k: v for k, v in duplicate['reference'].items() if k not in usage},
                             {k: v for k, v in first.items() if k not in usage})
            self.assertGreaterEqual(duplicate['reference']['last_used_at'], first['last_used_at'])
        self.assertEqual(len(list(self.uploads.iterdir())), 1)

    def test_same_filename_different_content_never_overwrites_existing_or_indexes(self):
        first = self.library.upload(uploaded(), 'First', False)['reference']
        index = self.uploads / (first['filename'] + '.nhr')
        index.write_bytes(b'original index')
        second = self.library.upload(uploaded(data=FASTA + b'ACGT\n'), 'Second', True)
        self.assertFalse(second['duplicate'])
        self.assertNotEqual(first['filename'], second['reference']['filename'])
        self.assertEqual((self.uploads / first['filename']).read_bytes(), FASTA)
        self.assertEqual(index.read_bytes(), b'original index')

    def test_legacy_files_are_listed_then_hashed_for_deduplication(self):
        old = self.uploads / 'legacy.fasta'
        old.write_bytes(FASTA)
        self.assertIsNone(self.library.list()[0]['sha256'])
        self.assertIsNone(self.library.list()[0]['polyploid'])
        duplicate = self.library.upload(uploaded('new.fa'), 'New upload', True)
        self.assertTrue(duplicate['duplicate'])
        self.assertEqual(duplicate['reference']['filename'], old.name)
        self.assertEqual(len(self.library.list()), 1)

    def test_unicode_filenames_remain_usable_and_external_paths_are_rejected(self):
        entry = self.library.upload(uploaded('中文参考.fa'), '中文参考', None)['reference']
        self.assertTrue(self.library.path(entry['filename']).is_file())
        self.assertEqual(len(self.library.list()), 1)
        for name in ('../outside.fa', str(self.root / 'outside.fa'), None):
            with self.assertRaises(ValueError):
                self.library.configure(name, 'Invalid', True)
        (self.uploads / 'symlink.fa').symlink_to(self.library.path(entry['filename']))
        self.assertEqual(len(self.library.list()), 1)

    def test_invalid_uploads_leave_no_partial_files(self):
        for name, data in [('bad.fa', b'ACGT'), ('empty.fa', b''), ('empty-seq.fa', b'>id\n'),
                           ('duplicates.fa', b'>id\nA\n>ID\nA\n'), ('protein.fa', b'>id\nLFEQ\n'),
                           ('bad.gz', b'bad'), ('broken.fa.gz', gzip.compress(FASTA)[:-5])]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.library.upload(uploaded(name, data), 'Bad', None)
            self.assertEqual(list(self.uploads.iterdir()), [])
            self.assertFalse(list(self.library.catalog_path.parent.glob('*.tmp')))

    def test_hash_cache_detects_same_size_rewrite_with_preserved_mtime(self):
        entry = self.library.upload(uploaded(), 'Original', None)['reference']
        path = self.library.path(entry['filename'])
        stat = path.stat()
        path.write_bytes(FASTA.replace(b'ACGT', b'TGCA'))
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        self.assertIsNone(self.library.list()[0]['sha256'])
        updated = self.library.index(path.name)
        self.assertNotEqual(updated['sha256'], entry['sha256'])
        self.assertEqual(updated['sequence_count'], 1)

    def test_deleted_references_disappear_and_can_be_reuploaded(self):
        entry = self.library.upload(uploaded(), 'First', True)['reference']
        self.library.path(entry['filename']).unlink()
        self.assertEqual(self.library.list(), [])
        self.assertFalse(self.library.upload(uploaded(), 'Again', False)['duplicate'])

    def test_annotation_upload_persistence_validation_and_reference_change(self):
        entry = self.library.upload(uploaded(), 'Bju', True,
                                    uploaded('groups.csv', b'transcript_id,subgenome\nBjuVA01G1,A\n'))['reference']
        self.assertTrue(entry['subgenome_ready'])
        self.assertEqual(entry['subgenome']['subgenome_counts'], {'A': 1})
        restarted = ReferenceLibrary(self.uploads, self.library.catalog_path)
        original, metadata = restarted.annotation(entry['filename'])
        self.assertIn(b'BjuVA01G1\tA', original)
        for data in (b'wrong,B\n', b'BjuVA01G1,A\nBjuVA01G1,B\n', b'BjuVA01G1\n'):
            with self.assertRaises(ValueError):
                restarted.configure(entry['filename'], 'Bad', False, uploaded('bad.csv', data))
            self.assertEqual(restarted.list()[0], entry)
            self.assertEqual(restarted.annotation(entry['filename'])[0], original)
        # Renamed identical FASTA must not silently replace an existing annotation.
        duplicate = restarted.upload(uploaded('renamed.fa'), 'Duplicate', False, uploaded('other.csv', b'BjuVA01G1,B\n'))
        self.assertTrue(duplicate['duplicate'])
        self.assertEqual(restarted.annotation(entry['filename'])[0], original)
        restarted.path(entry['filename']).write_bytes(FASTA + b'ACGT\n')
        self.assertFalse(restarted.list()[0]['subgenome_ready'])
        with self.assertRaises(ValueError):
            restarted.annotation(entry['filename'])

    def test_annotation_and_settings_roll_back_together_on_save_error(self):
        entry = self.library.upload(uploaded(), 'Bju', True, uploaded('groups.csv', b'BjuVA01G1,A\n'))['reference']
        original, _ = self.library.annotation(entry['filename'])
        with patch.object(self.library, '_save', side_effect=OSError('disk full')), self.assertRaises(OSError):
            self.library.configure(entry['filename'], 'Updated', False, uploaded('new.csv', b'BjuVA01G1,B\n'))
        self.assertEqual(self.library.list()[0], entry)
        self.assertEqual(self.library.annotation(entry['filename'])[0], original)

    def test_polyploid_configuration_is_saved_in_prepared_task(self):
        entry = self.library.upload(uploaded(), 'Bju configured', True)['reference']
        service = SplitWorkflow(self.uploads, {'BP_1': 'A' * 19}, {}, {}, {}, reference_library=self.library)
        preview = service.prepare(dict(input_mode='single', name='test', geneID='BjuVA01G1', bp_id='BP_1',
                                       sequenceType='reference', reference_mode='existing', existing_ref_genome=entry['filename']), {})
        self.assertTrue(preview['reference']['polyploid'])
        self.assertEqual(preview['reference']['display_name'], 'Bju configured')


class ReferenceHttpTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.app = Flask(__name__)
        self.app.config.update(TESTING=True, UPLOAD_FOLDER=str(self.root / 'uploads'),
                               ADMIN_CREDENTIALS_FILE=str(self.root / 'admin.json'),
                               REFERENCE_CATALOG_FILE=str(self.root / 'references.json'))
        register_admin(self.app)
        register_reference_routes(self.app)
        self.client = self.app.test_client()

    def tearDown(self):
        self.tmp.cleanup()

    def login(self):
        result = self.app.test_cli_runner().invoke(args=['init-admin', '--username', 'admin'], input='test-reference-password\ntest-reference-password\n')
        self.assertEqual(result.exit_code, 0)
        # Obtain the CSRF token without needing the full application's template shell.
        with self.client.session_transaction() as session:
            session['admin_csrf'] = 'test-token'
        self.assertEqual(self.client.post('/admin/login', data=dict(username='admin', password='test-reference-password',
                                                                 csrf_token='test-token')).status_code, 302)
        with self.client.session_transaction() as session:
            session['admin_csrf'] = 'logged-in-token'
        return {'X-CSRF-Token': 'logged-in-token'}

    def test_public_upload_configuration_and_admin_only_deletion(self):
        self.assertEqual(self.client.get('/references', headers={'Accept':'application/json'}).status_code, 200)
        for route in ('/references/upload', '/references/configure', '/references/hash'):
            self.assertEqual(self.client.post(route, json={}).status_code, 400)  # still requires CSRF
        with self.client.session_transaction() as session:
            session['admin_csrf'] = 'public-token'
        headers = {'X-CSRF-Token':'public-token'}
        self.assertEqual(self.client.post('/references/upload').status_code, 400)
        response = self.client.post('/references/upload', headers=headers, content_type='multipart/form-data', data=dict(
            reference=(io.BytesIO(FASTA), 'reference.fa'), display_name='Bju', polyploid='true'))
        self.assertEqual(response.status_code, 201, response.json)
        filename = response.json['reference']['filename']
        annotated = self.client.post('/references/configure', headers=headers, data=dict(
            filename=filename, display_name='Bju', polyploid='true',
            subgenome_map=(io.BytesIO(b'BjuVA01G1,A\n'), 'groups.csv')))
        self.assertEqual(annotated.status_code, 200, annotated.json)
        self.assertTrue(annotated.json['reference']['subgenome_ready'])
        duplicate = self.client.post('/references/upload', headers=headers, content_type='multipart/form-data', data=dict(
            reference=(io.BytesIO(gzip.compress(FASTA)), 'renamed.fa.gz'), display_name='Duplicate', polyploid='false'))
        self.assertEqual(duplicate.status_code, 200)
        self.assertTrue(duplicate.json['duplicate'])
        configured = self.client.post('/references/configure', headers=headers, json=dict(filename=filename, display_name='Renamed', polyploid=False))
        self.assertEqual(configured.status_code, 200)
        self.assertFalse(configured.json['reference']['polyploid'])
        self.assertEqual(self.client.post('/references/hash', headers=headers, json=dict(filename=filename)).status_code, 200)
        self.assertEqual(self.client.post('/admin/delete', headers=headers, json=dict(keys=[filename])).status_code, 401)
        with self.client.session_transaction() as session:
            session.clear()
        listing = self.client.get('/references', headers={'Accept':'application/json'}).json
        self.assertTrue(listing['can_edit'])
        self.assertEqual(listing['references'][0]['display_name'], 'Renamed')
        headers = self.login()
        removed = self.client.post('/admin/delete', headers=headers, json=dict(keys=[filename]))
        self.assertEqual(removed.status_code, 200)
        self.assertEqual(self.app.extensions['reference_library'].list(), [])
        self.assertFalse((self.root / 'uploads' / (filename + '.subgenomes.tsv')).exists())

    def test_cli_dry_run_and_admin_inventory_show_unused_references(self):
        uploads = self.root / 'uploads'
        uploads.mkdir()
        library = self.app.extensions['reference_library']
        self.assertIs(self.app.extensions['storage_maintenance'].reference_library, library)
        now = datetime.now(timezone.utc)
        with patch('scripts.reference_library.utcnow', return_value=now - timedelta(days=8)):
            stale = library.upload(uploaded(), 'Stale', None)['reference']['filename']
        fresh = library.upload(uploaded('fresh.fa', FASTA.replace(b'G1', b'G2')), 'Fresh', None)['reference']['filename']
        runner = self.app.test_cli_runner()
        preview = runner.invoke(args=['cleanup-results', '--dry-run'])
        self.assertEqual(preview.exit_code, 0, preview.output)
        self.assertEqual(json.loads(preview.output)['deleted'], [stale])
        self.assertTrue((uploads / stale).exists())
        self.login()
        inventory = self.client.get('/admin', headers={'Accept': 'application/json'}).json
        self.assertEqual(inventory['reference_retention_days'], 7)
        entries = {e['key']: e for e in inventory['entries']}
        self.assertTrue(entries[stale]['expired'])
        self.assertFalse(entries[fresh]['expired'])
        expires = datetime.fromisoformat(entries[fresh]['expires_at'])
        self.assertLess(abs(expires - (now + timedelta(days=7))), timedelta(minutes=1))
        actual = runner.invoke(args=['cleanup-results'])
        self.assertEqual(actual.exit_code, 0, actual.output)
        self.assertEqual([e['filename'] for e in library.list()], [fresh])
        self.assertNotIn(stale, json.loads(Path(self.app.config['REFERENCE_CATALOG_FILE']).read_text()))

    def test_optional_admin_only_configuration(self):
        self.app.config['REFERENCE_PUBLIC_EDIT'] = False
        for route in ('/references/upload', '/references/configure', '/references/hash'):
            self.assertEqual(self.client.post(route, json={}).status_code, 401)


if __name__ == '__main__':
    unittest.main()


T0 = datetime(2026, 11, 1, 12, tzinfo=timezone.utc)


class ReferenceRetentionTests(unittest.TestCase):
    """Shared references unused for REFERENCE_RETENTION_DAYS are removed by the retention job."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.uploads = self.root / 'uploads'
        self.uploads.mkdir()
        self.library = ReferenceLibrary(self.uploads, self.root / 'instance' / 'references.json')
        self.storage = StorageMaintenance(self.uploads, reference_library=self.library)

    def tearDown(self):
        self.tmp.cleanup()

    def upload(self, when=T0, **kwargs):
        with patch('scripts.reference_library.utcnow', return_value=when):
            return self.library.upload(uploaded(**kwargs), 'Bju', True)['reference']['filename']

    def entry(self, key, now):
        return next(e for e in self.storage.inventory(now) if e['key'] == key)

    def test_period_is_seven_days_unless_overridden_by_environment(self):
        self.assertEqual(REFERENCE_RETENTION_DAYS, 7)
        self.assertEqual(self.storage.reference_retention, timedelta(days=7))
        for value, days in (('3', 3), ('0', 7), ('-2', 7), ('soon', 7)):
            with patch.dict(os.environ, {'HCR_REFERENCE_RETENTION_DAYS': value}):
                self.assertEqual(StorageMaintenance(self.uploads).reference_retention, timedelta(days=days))

    def test_fresh_upload_is_kept_and_shows_its_expiry(self):
        filename = self.upload()
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=1))['deleted'], [])
        self.assertTrue((self.uploads / filename).exists())
        entry = self.entry(filename, T0 + timedelta(days=1))
        self.assertEqual(entry['modified_at'], T0.isoformat())
        self.assertEqual(entry['expires_at'], (T0 + timedelta(days=7)).isoformat())
        self.assertFalse(entry['expired'])

    def test_seven_day_boundary_and_dry_run(self):
        filename = self.upload()
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=7, seconds=-1))['deleted'], [])
        preview = self.storage.cleanup(dry_run=True, now=T0 + timedelta(days=7))
        self.assertEqual(preview['deleted'], [filename])
        self.assertTrue((self.uploads / filename).exists())
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=7))['deleted'], [filename])
        self.assertFalse((self.uploads / filename).exists())

    def test_uses_reset_the_clock(self):
        filename = self.upload()
        service = SplitWorkflow(self.uploads, {'BP_1': 'A' * 19}, {}, {}, {}, reference_library=self.library)
        with patch('scripts.reference_library.utcnow', return_value=T0 + timedelta(days=5)):
            prepared = service.prepare(dict(input_mode='single', name='test', geneID='BjuVA01G1', bp_id='BP_1',
                                            sequenceType='reference', reference_mode='existing',
                                            existing_ref_genome=filename), {})
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=10))['deleted'], [])
        with patch('scripts.reference_library.utcnow', return_value=T0 + timedelta(days=9)):
            self.assertTrue(service.queue(prepared['id'], {}))
        manifest_path = self.uploads / 'split_results' / prepared['id'] / 'request.json'
        manifest = json.loads(manifest_path.read_text())
        manifest['status'] = 'completed'
        manifest_path.write_text(json.dumps(manifest))
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=15, seconds=-1))['deleted'], [])
        # A re-upload of the same content (also the legacy triplet path) and an explicit touch count as uses.
        self.upload(when=T0 + timedelta(days=14), name='again.fa')
        self.assertTrue(self.library.touch(filename, now=T0 + timedelta(days=20)))
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=27, seconds=-1))['deleted'], [])
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=27))['deleted'], [filename])

    def test_references_in_use_symlinks_and_damaged_manifests_are_protected(self):
        filename = self.upload()
        later = T0 + timedelta(days=30)
        job = self.uploads / 'split_results' / ('a' * 32)
        job.mkdir(parents=True)
        manifest = dict(id='a' * 32, status='queued', created_at=later.isoformat(),
                        reference_path=str(self.uploads / filename))
        (job / 'request.json').write_text(json.dumps(manifest))
        report = self.storage.cleanup(now=later)
        self.assertIn(filename, report['skipped'])
        self.assertTrue((self.uploads / filename).exists())
        (job / 'request.json').write_text('broken')
        self.assertIn(filename, self.storage.cleanup(now=later)['skipped'])
        (job / 'request.json').write_text(json.dumps(dict(manifest, status='completed')))
        self.storage.legacy_active = 1  # a triplet request is in flight
        self.assertIn(filename, self.storage.cleanup(now=later)['skipped'])
        self.storage.legacy_active = 0
        outside = self.root / 'outside.fa'
        outside.write_bytes(FASTA)
        (self.uploads / 'linked.fa').symlink_to(outside)
        report = self.storage.cleanup(now=later)
        self.assertEqual(report['deleted'], [filename])
        self.assertTrue(outside.exists() and (self.uploads / 'linked.fa').is_symlink())

    def test_legacy_reference_without_record_gets_a_grace_period(self):
        legacy = self.uploads / 'legacy.fasta'
        legacy.write_bytes(FASTA)
        os.utime(legacy, (0, 0))  # file times never decide expiry
        self.assertEqual(self.storage.cleanup(now=T0)['deleted'], [])
        self.assertEqual(self.library._read()['legacy.fasta']['last_used_at'], T0.isoformat())
        self.assertEqual(self.storage.cleanup(dry_run=True, now=T0 + timedelta(days=6))['deleted'], [])
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=7))['deleted'], ['legacy.fasta'])

    def test_damaged_catalog_stops_automatic_expiry(self):
        filename = self.upload()
        self.library.catalog_path.write_text('[broken')
        self.assertFalse(self.library.touch(filename))
        self.assertEqual(self.storage.cleanup(now=T0 + timedelta(days=30))['deleted'], [])
        self.assertTrue((self.uploads / filename).exists())

    def test_expiry_removes_companion_files_and_catalog_entry(self):
        with patch('scripts.reference_library.utcnow', return_value=T0):
            filename = self.library.upload(uploaded(), 'Bju', True, uploaded('groups.csv', b'BjuVA01G1,A\n'))['reference']['filename']
        companions = [filename + '.gz', filename + '.nhr', filename + '.00.nsq', filename + '.subgenomes.tsv']
        for name in companions[:-1]:
            (self.uploads / name).write_text('index')
        kept = self.upload(when=T0 + timedelta(days=6), name='other.fa', data=FASTA.replace(b'G1', b'G2'))
        self.assertEqual(self.entry(filename, T0)['file_count'], 5)
        report = self.storage.cleanup(now=T0 + timedelta(days=7))
        self.assertEqual(report['deleted'], [filename])
        self.assertFalse(any((self.uploads / name).exists() for name in [filename] + companions))
        self.assertNotIn(filename, self.library._read())
        self.assertIn(kept, self.library._read())
        self.assertEqual([e['filename'] for e in self.library.list()], [kept])


class ReferenceExpiryListingTests(unittest.TestCase):
    """Pages show when the retention job may delete each reference, without writing to the catalog."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.uploads = self.root / 'uploads'
        self.uploads.mkdir()
        self.library = ReferenceLibrary(self.uploads, self.root / 'instance' / 'references.json')
        with patch('scripts.reference_library.utcnow', return_value=T0):
            self.filename = self.library.upload(uploaded(), 'Bju', True)['reference']['filename']

    def tearDown(self):
        self.tmp.cleanup()

    def entry(self, now):
        return next(e for e in self.library.list(now=now) if e['filename'] == self.filename)

    def test_listing_fields_and_boundary_match_the_cleanup_job(self):
        storage = StorageMaintenance(self.uploads, reference_library=self.library)
        entry = self.entry(T0 + timedelta(days=4, hours=12))
        self.assertEqual(entry['last_used_at'], T0.isoformat())
        self.assertEqual(entry['expires_at'], (T0 + timedelta(days=7)).isoformat())
        self.assertEqual((entry['days_left'], entry['expiring_soon']), (3, False))
        self.assertEqual(self.entry(T0 + timedelta(days=6))['expiring_soon'], True)
        self.assertEqual(self.entry(T0 + timedelta(days=6))['days_left'], 1)
        # Exactly the retention period: the job deletes it, and the listing says so.
        boundary = T0 + timedelta(days=REFERENCE_RETENTION_DAYS)
        self.assertEqual((self.entry(boundary)['days_left'], self.entry(boundary)['expiring_soon']), (0, True))
        self.assertTrue(next(e for e in storage.inventory(boundary) if e['key'] == self.filename)['expired'])
        self.assertEqual(self.entry(boundary - timedelta(seconds=1))['days_left'], 1)

    def test_listing_follows_the_retention_setting_and_never_writes(self):
        catalog = self.library.catalog_path
        records = json.loads(catalog.read_text())
        del records[self.filename]['last_used_at']
        catalog.write_text(json.dumps(records))
        before = catalog.read_bytes()
        with patch.dict(os.environ, {'HCR_REFERENCE_RETENTION_DAYS': '3'}):
            entry = self.entry(T0)
        # No usage record yet: the job starts the clock on its next run, so the full period remains.
        self.assertEqual((entry['last_used_at'], entry['days_left']), (None, 3))
        self.assertEqual(entry['expires_at'], (T0 + timedelta(days=3)).isoformat())
        self.assertEqual(catalog.read_bytes(), before)
        self.library.retention = timedelta(days=10)
        self.assertEqual(self.entry(T0)['days_left'], 10)

    def test_damaged_entries_and_catalog_show_no_expiry(self):
        catalog = self.library.catalog_path
        records = json.loads(catalog.read_text())
        for damaged in ('not a date', 12, 'damaged entry'):
            if damaged == 'damaged entry':
                records[self.filename] = ['not', 'a', 'mapping']
            else:
                records[self.filename]['last_used_at'] = damaged
            catalog.write_text(json.dumps(records))
            with self.subTest(damaged=damaged):
                entry = self.entry(T0)
                self.assertEqual((entry['expires_at'], entry['days_left'], entry['last_used_at']), (None, None, None))
        catalog.write_text('[broken')
        entry = self.entry(T0)
        self.assertEqual((entry['filename'], entry['days_left'], entry['subgenome_ready']), (self.filename, None, False))


class ReferenceExpiryPageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import app as web
        cls.web = web

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / 'uploads').mkdir()
        self.library = ReferenceLibrary(self.root / 'uploads', self.root / 'references.json')
        with patch.dict(os.environ, {'HCR_REFERENCE_RETENTION_DAYS': '5'}):
            self.library.retention = StorageMaintenance(self.root / 'uploads').reference_retention
        patcher = patch.dict(self.web.app.extensions, reference_library=self.library)
        patcher.start()
        self.addCleanup(patcher.stop)
        now = datetime.now(timezone.utc)
        with patch('scripts.reference_library.utcnow', return_value=now - timedelta(days=2, hours=1)):
            self.library.upload(uploaded('fresh.fa'), 'Fresh', False)
        with patch('scripts.reference_library.utcnow', return_value=now - timedelta(days=4, hours=12)):
            self.library.upload(uploaded('old.fa', FASTA.replace(b'G1', b'G2')), 'Old', None)

    def tearDown(self):
        self.tmp.cleanup()

    def render(self, template, language):
        if template != 'references.html':
            return self.web.app.test_client().get({'split_wizard.html': '/split', 'triplet.html': '/triplet'}[template]
                                                  + '?lang=' + language).text
        from flask import render_template
        with self.web.app.test_request_context('/references?lang=' + language):
            return render_template(template, references=self.library.list(), can_edit=False, can_delete=False)

    def test_pickers_show_days_left_and_the_configured_period_in_both_languages(self):
        self.assertEqual(self.library.retention, timedelta(days=5))
        for template in ('split_wizard.html', 'triplet.html'):
            zh, en = self.render(template, 'zh'), self.render(template, 'en')
            with self.subTest(template=template):
                self.assertIn('Fresh · 3 天后过期</option>', zh)
                self.assertIn('class="reference-expiring" data-expiring="true"', zh)
                self.assertIn('Old · ⚠ 1 天内过期</option>', zh)
                self.assertIn('Fresh · expires in 3 days</option>', en)
                self.assertIn('Old · ⚠ expires within a day</option>', en)
        self.assertIn('重新开始 5 天计时；5 天未用的参考会被自动删除', self.render('split_wizard.html', 'zh'))
        self.assertIn("restarts its 5-day clock; references unused for 5 days", self.render('split_wizard.html', 'en'))

    def test_reference_page_rows_show_last_use_and_expiry(self):
        zh, en = self.render('references.html', 'zh'), self.render('references.html', 'en')
        self.assertEqual(zh.count('data-local-time'), 4)
        self.assertIn('最近使用：', zh)
        self.assertIn('过期时间：', zh)
        self.assertIn('reference-expiry is-expiring', zh)
        self.assertIn('Last used: ', en)
        self.assertIn('expires in 3 days', en)
        self.library.catalog_path.write_text('{"fresh.fa": "broken"}')
        self.assertIn('Usage record unreadable; expiry not shown', self.render('references.html', 'en'))
