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
            self.assertEqual(duplicate['reference'], first)
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

    def test_optional_admin_only_configuration(self):
        self.app.config['REFERENCE_PUBLIC_EDIT'] = False
        for route in ('/references/upload', '/references/configure', '/references/hash'):
            self.assertEqual(self.client.post(route, json={}).status_code, 401)


if __name__ == '__main__':
    unittest.main()
