import glob
import json
import os
import random
import subprocess
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import utils
from scripts.storage_maintenance import StorageMaintenance
from scripts.utils import BLAST_INDEX_RECORD, create_blast_db
from src.common.blast_runtime import blast_executable
from src.common.sequence_utils import check_left_probe_gc_junction


class TestGCJunctionCheck(unittest.TestCase):
    """测试Left探针连接处G/C碱基检查功能（影响SplintR酶活性）

    只有前两个碱基都是 G/C 才算硬性过滤；单个 G/C 开头由筛选阶段的排序惩罚处理。
    """

    def test_two_leading_gc_bases_are_a_junction(self):
        """GG/GC/CG/CC 开头：严重抑制SplintR连接，直接过滤"""
        for start in ("GG", "GC", "CG", "CC"):
            with self.subTest(start=start):
                self.assertTrue(check_left_probe_gc_junction(start + "CATGCATGATCG"))

    def test_single_leading_gc_base_is_not_a_junction(self):
        """只有第一个或第二个碱基是 G/C：不做硬性过滤"""
        for start in ("GA", "GT", "CA", "CT", "AG", "AC", "TG", "TC"):
            with self.subTest(start=start):
                self.assertFalse(check_left_probe_gc_junction(start + "CATGCATGATCG"))

    def test_at_junction(self):
        """A/T开头的探针序列（有利于SplintR连接）"""
        for start in ("AT", "AA", "TA", "TT"):
            with self.subTest(start=start):
                self.assertFalse(check_left_probe_gc_junction(start + "CATGCATGATCG"))

    def test_sequences_shorter_than_two_bases(self):
        """不足两个碱基时无法判断连接处，按未命中处理"""
        for sequence in ("", "G", "C", "A", "T"):
            with self.subTest(sequence=sequence):
                self.assertFalse(check_left_probe_gc_junction(sequence))


SINGLE_VOLUME = ('.ndb', '.nhr', '.nin', '.njs', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto')
MULTI_VOLUME = ('.00.nhr', '.00.nin', '.00.nsq', '.01.nhr', '.01.nin', '.01.nsq', '.nal', '.ndb', '.njs')


class FakeMakeblastdb:
    """Stands in for NcbimakeblastdbCommandline and writes index files at ``out``."""

    def __init__(self, suffixes=SINGLE_VOLUME, delay=0.0, fail_after=None):
        self.suffixes, self.delay, self.fail_after = suffixes, delay, fail_after
        self.calls, self.running, self.peak = [], 0, 0
        self.guard = threading.Lock()

    def __call__(self, **kwargs):
        def run():
            with self.guard:
                self.calls.append(kwargs)
                self.running += 1
                self.peak = max(self.peak, self.running)
            try:
                time.sleep(self.delay)
                with open(kwargs['input_file'], 'rb') as handle:
                    content = handle.read()
                for number, suffix in enumerate(self.suffixes):
                    if self.fail_after is not None and number == self.fail_after:
                        raise RuntimeError('makeblastdb interrupted')
                    with open(kwargs['out'] + suffix, 'wb') as handle:
                        handle.write(content + suffix.encode())
            finally:
                with self.guard:
                    self.running -= 1
        return run


def write_fasta(path, body):
    with open(path, 'w') as handle:
        handle.write(body)


def index_names(path):
    base = os.path.basename(path)
    return sorted(name[len(base):] for name in os.listdir(os.path.dirname(path))
                  if name.startswith(base) and name != base)


class TestCreateBlastDb(unittest.TestCase):
    """Index building accepts what the reference library accepts on upload."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.fasta = os.path.join(self.tmp.name, 'ref.fasta')
        write_fasta(self.fasta, '>gene1\nACGTACGT\n')

    def fake(self, **kwargs):
        fake = FakeMakeblastdb(**kwargs)
        for patcher in (patch('scripts.utils.NcbimakeblastdbCommandline', side_effect=fake),
                        patch('scripts.utils.blast_executable', return_value='makeblastdb')):
            patcher.start()
            self.addCleanup(patcher.stop)
        return fake

    def record(self):
        with open(self.fasta + BLAST_INDEX_RECORD) as handle:
            return json.load(handle)

    def test_lowercase_and_ambiguity_codes_are_accepted_without_rewriting_the_file(self):
        fake = self.fake()
        body = ">Mp1g00010.1\natgcNNryswkmbdhvn\n>AT1G01010\nacgtACGTU\n"
        write_fasta(self.fasta, body)
        self.assertEqual(create_blast_db(self.fasta), self.fasta)
        self.assertEqual(len(fake.calls), 1)
        self.assertTrue(fake.calls[0]['parse_seqids'])
        with open(self.fasta) as handle:
            self.assertEqual(handle.read(), body)  # the recorded SHA-256 must stay valid

    def test_non_nucleotide_and_empty_records_raise(self):
        fake = self.fake()
        for body, message in [(">prot\nMEEPQSDL\n", "Sequence prot contains invalid DNA characters"),
                              (">empty\n>next\nACGT\n", "Sequence empty is empty"),
                              ("", "no sequences")]:
            with self.subTest(body=body):
                write_fasta(self.fasta, body)
                with self.assertRaisesRegex(ValueError, message):
                    create_blast_db(self.fasta)
        self.assertEqual(fake.calls, [])
        self.assertEqual(index_names(self.fasta), [])

    def test_unchanged_fasta_reuses_index_without_hashing(self):
        fake = self.fake()
        create_blast_db(self.fasta)
        self.assertEqual(index_names(self.fasta), sorted(SINGLE_VOLUME + (BLAST_INDEX_RECORD,)))
        record = self.record()
        self.assertTrue(record['parse_seqids'])
        self.assertEqual(sorted(record['files']), sorted('ref.fasta' + s for s in SINGLE_VOLUME))
        with patch('scripts.utils._sha256', side_effect=AssertionError('hashed')):
            for _ in range(3):
                create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 1)

    def test_replaced_or_edited_fasta_rebuilds(self):
        """Failure mode 1: the FASTA changes while the old index files remain."""
        fake = self.fake()
        create_blast_db(self.fasta)
        stat = os.stat(self.fasta)
        # Replaced by rename with the same size and the same timestamps (cp -p, rsync -t).
        replacement = os.path.join(self.tmp.name, 'incoming')
        write_fasta(replacement, '>gene2\nTTTTTTTT\n')
        os.utime(replacement, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        os.replace(replacement, self.fasta)
        create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 2)
        # Edited in place, different size.
        write_fasta(self.fasta, '>gene3\nACGTACGTAAAA\n')
        create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 3)
        with open(self.fasta + '.nsq', 'rb') as handle:
            self.assertTrue(handle.read().startswith(b'>gene3'))
        self.assertEqual(self.record()['fasta']['size'], os.path.getsize(self.fasta))

    def test_touched_fasta_with_same_content_is_hashed_once_not_rebuilt(self):
        fake = self.fake()
        create_blast_db(self.fasta)
        os.utime(self.fasta, ns=(1_000_000_000, 1_000_000_000))
        with patch('scripts.utils._sha256', wraps=utils._sha256) as digest:
            create_blast_db(self.fasta)
            create_blast_db(self.fasta)
        self.assertEqual(digest.call_count, 1)
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(self.record()['fasta']['mtime_ns'], 1_000_000_000)

    def test_legacy_index_without_record_is_rebuilt_once(self):
        """Failure mode 2: an older version's index (maybe without -parse_seqids) has no record."""
        for suffix in ('.nhr', '.nin', '.nsq', '.nal', '.00.nsq'):
            write_fasta(self.fasta + suffix, 'legacy')
        fake = self.fake()
        create_blast_db(self.fasta)
        create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(index_names(self.fasta), sorted(SINGLE_VOLUME + (BLAST_INDEX_RECORD,)))
        self.assertNotIn(b'legacy', Path(self.fasta + '.nhr').read_bytes())

    def test_record_without_parse_seqids_or_unknown_format_is_rebuilt(self):
        fake = self.fake()
        create_blast_db(self.fasta)
        for change in ({'parse_seqids': False}, {'format': 0}):
            record = self.record()
            record.update(change)
            write_fasta(self.fasta + BLAST_INDEX_RECORD, json.dumps(record))
            create_blast_db(self.fasta)
        write_fasta(self.fasta + BLAST_INDEX_RECORD, 'not json')
        create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 4)

    def test_interrupted_build_never_looks_complete(self):
        """Failure mode 3: makeblastdb dies part-way; its files must not be trusted later."""
        self.fake(fail_after=4)
        with self.assertRaisesRegex(RuntimeError, 'interrupted'):
            create_blast_db(self.fasta)
        self.assertEqual(index_names(self.fasta), [])  # nothing half-written next to the FASTA
        self.assertEqual([n for n in os.listdir(self.tmp.name) if n.startswith('.')], [])
        # A build directory left by a killed process is cleared by the next build.
        os.mkdir(os.path.join(self.tmp.name, '.ref.fasta.blastdb-build-dead'))
        fake = self.fake()
        create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(sorted(os.listdir(self.tmp.name)),
                         sorted(['ref.fasta'] + ['ref.fasta' + s for s in SINGLE_VOLUME + (BLAST_INDEX_RECORD,)]))

    def test_failed_rebuild_keeps_no_record_for_old_index(self):
        self.fake()
        create_blast_db(self.fasta)
        write_fasta(self.fasta, '>gene2\nTTTTTTTTTT\n')
        self.fake(fail_after=2)
        with self.assertRaises(RuntimeError):
            create_blast_db(self.fasta)
        fake = self.fake()
        create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 1)

    def test_missing_or_resized_index_file_is_rebuilt(self):
        fake = self.fake()
        create_blast_db(self.fasta)
        os.unlink(self.fasta + '.nsq')
        create_blast_db(self.fasta)
        with open(self.fasta + '.nhr', 'ab') as handle:
            handle.write(b'truncated elsewhere')
        create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 3)

    def test_multi_volume_index_is_recognised(self):
        """Failure mode 4: .00.nhr + .nal indexes used to be rebuilt on every call."""
        write_fasta(self.fasta + '.nhr', 'stale single volume')
        fake = self.fake(suffixes=MULTI_VOLUME)
        for _ in range(3):
            create_blast_db(self.fasta)
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(index_names(self.fasta), sorted(MULTI_VOLUME + (BLAST_INDEX_RECORD,)))
        # Rebuilding to a single volume drops the old volumes.
        write_fasta(self.fasta, '>gene2\nTTTTTTTTTT\n')
        self.fake()
        create_blast_db(self.fasta)
        self.assertEqual(index_names(self.fasta), sorted(SINGLE_VOLUME + (BLAST_INDEX_RECORD,)))

    def test_concurrent_requests_build_once(self):
        """Failure mode 5: simultaneous requests ran makeblastdb on the same output at once."""
        fake = self.fake(delay=0.2)
        errors = []

        def call():
            try:
                create_blast_db(self.fasta)
            except Exception as error:  # pragma: no cover - reported below
                errors.append(error)
        threads = [threading.Thread(target=call) for _ in range(4)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(errors, [])
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(fake.peak, 1)

    def test_fasta_changed_during_build_is_not_recorded(self):
        fasta = self.fasta

        class Racing(FakeMakeblastdb):
            def __call__(self, **kwargs):
                run = super().__call__(**kwargs)

                def racing():
                    run()
                    write_fasta(fasta, '>late\nACGTACGTACGTACGT\n')
                return racing
        self.fake()
        patcher = patch('scripts.utils.NcbimakeblastdbCommandline', side_effect=Racing())
        patcher.start()
        self.addCleanup(patcher.stop)
        with self.assertRaisesRegex(RuntimeError, 'changed while'):
            create_blast_db(self.fasta)
        self.assertFalse(os.path.exists(self.fasta + BLAST_INDEX_RECORD))

    def test_storage_inventory_deletes_record_with_reference(self):
        self.fake()
        create_blast_db(self.fasta)
        storage = StorageMaintenance(self.tmp.name)
        entry = next(e for e in storage.inventory() if e['key'] == 'ref.fasta')
        self.assertEqual(entry['file_count'], 1 + len(SINGLE_VOLUME) + 1)
        storage.delete(['ref.fasta'])
        self.assertEqual(os.listdir(self.tmp.name), [])


def blast_available():
    try:
        blast_executable('makeblastdb')
        blast_executable('blastdbcmd')
        return True
    except FileNotFoundError:
        return False


@unittest.skipUnless(blast_available(), 'BLAST+ is not installed')
class TestCreateBlastDbWithBlast(unittest.TestCase):
    """The same failure modes against real makeblastdb/blastdbcmd."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.fasta = os.path.join(self.tmp.name, 'ref.fasta')
        self.random = random.Random(7)

    def write(self, names, length=300):
        with open(self.fasta, 'w') as handle:
            for name in names:
                handle.write(f'>{name}\n' + ''.join(self.random.choices('ACGT', k=length)) + '\n')

    def accessions(self):
        result = subprocess.run([blast_executable('blastdbcmd'), '-db', self.fasta, '-entry', 'all',
                                 '-outfmt', '%a'], capture_output=True, text=True, check=True)
        return sorted(result.stdout.split())

    def test_real_index_tracks_fasta_and_ids(self):
        # Legacy index from an older version: no -parse_seqids, no record.
        self.write(['geneA'])
        subprocess.run([blast_executable('makeblastdb'), '-in', self.fasta, '-dbtype', 'nucl',
                        '-out', self.fasta], capture_output=True, check=True)
        self.assertNotEqual(self.accessions(), ['geneA'])
        create_blast_db(self.fasta)
        self.assertEqual(self.accessions(), ['geneA'])
        # Replaced reference: the index follows the new IDs.
        self.write(['geneB', 'geneC'])
        create_blast_db(self.fasta)
        self.assertEqual(self.accessions(), ['geneB', 'geneC'])
        # Partial index (empty files) with no record is rebuilt, not trusted.
        os.unlink(self.fasta + BLAST_INDEX_RECORD)
        for suffix in ('.nhr', '.nin', '.nsq'):
            open(self.fasta + suffix, 'wb').close()
        create_blast_db(self.fasta)
        self.assertEqual(self.accessions(), ['geneB', 'geneC'])

    def test_real_multi_volume_index_survives_the_move(self):
        real = utils.NcbimakeblastdbCommandline
        calls = []

        def small_volumes(**kwargs):
            calls.append(kwargs)
            return real(max_file_sz='1MB', **kwargs)
        self.write([f'g{i}' for i in range(1500)], length=5000)
        with patch('scripts.utils.NcbimakeblastdbCommandline', side_effect=small_volumes):
            create_blast_db(self.fasta)
            create_blast_db(self.fasta)
        self.assertEqual(len(calls), 1)
        self.assertTrue(os.path.exists(self.fasta + '.nal'))
        self.assertTrue(glob.glob(self.fasta + '.01.*'))
        accessions = self.accessions()
        self.assertEqual(len(accessions), 1500)
        self.assertIn('g1499', accessions)


if __name__ == "__main__":
    unittest.main()
