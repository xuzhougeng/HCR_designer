import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from scripts.utils import create_blast_db
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


class TestCreateBlastDb(unittest.TestCase):
    """Index building accepts what the reference library accepts on upload."""

    @patch("scripts.utils.NcbimakeblastdbCommandline")
    def test_lowercase_and_ambiguity_codes_are_accepted_without_rewriting_the_file(self, mock_makeblastdb):
        mock_makeblastdb.return_value = MagicMock()
        body = ">Mp1g00010.1\natgcNNryswkmbdhvn\n>AT1G01010\nacgtACGTU\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "input.fasta")
            with open(fasta_path, "w") as handle:
                handle.write(body)
            self.assertEqual(create_blast_db(fasta_path), fasta_path)
            mock_makeblastdb.assert_called_once()
            with open(fasta_path) as handle:
                self.assertEqual(handle.read(), body)  # the recorded SHA-256 must stay valid

    @patch("scripts.utils.NcbimakeblastdbCommandline")
    def test_non_nucleotide_and_empty_records_raise(self, mock_makeblastdb):
        mock_makeblastdb.return_value = MagicMock()
        for body, message in [(">prot\nMEEPQSDL\n", "Sequence prot contains invalid DNA characters"),
                              (">empty\n>next\nACGT\n", "Sequence empty is empty")]:
            with self.subTest(body=body), tempfile.TemporaryDirectory() as tmpdir:
                fasta_path = os.path.join(tmpdir, "input.fasta")
                with open(fasta_path, "w") as handle:
                    handle.write(body)
                with self.assertRaisesRegex(ValueError, message):
                    create_blast_db(fasta_path)
        mock_makeblastdb.assert_not_called()


if __name__ == "__main__":
    unittest.main()
