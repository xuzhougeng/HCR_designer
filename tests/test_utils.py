import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from scripts.utils import create_blast_db
from src.common.sequence_utils import check_left_probe_gc_junction


class TestGCJunctionCheck(unittest.TestCase):
    """测试Left探针连接处G/C碱基检查功能（影响SplintR酶活性）"""

    def test_gc_junction_G_first(self):
        """测试G开头的探针序列（不利于SplintR连接）"""
        self.assertTrue(check_left_probe_gc_junction("GGCATGCATGATCG"))
        self.assertTrue(check_left_probe_gc_junction("GCCATGCATGATCG"))
        self.assertTrue(check_left_probe_gc_junction("GACATGCATGATCG"))
        self.assertTrue(check_left_probe_gc_junction("GTCATGCATGATCG"))

    def test_gc_junction_C_first(self):
        """测试C开头的探针序列（不利于SplintR连接）"""
        self.assertTrue(check_left_probe_gc_junction("CGCATGCATGATCG"))
        self.assertTrue(check_left_probe_gc_junction("CCCATGCATGATCG"))
        self.assertTrue(check_left_probe_gc_junction("CACATGCATGATCG"))
        self.assertTrue(check_left_probe_gc_junction("CTCATGCATGATCG"))

    def test_at_junction_A_first(self):
        """测试A开头的探针序列（有利于SplintR连接）"""
        self.assertFalse(check_left_probe_gc_junction("AGCATGCATGATCG"))
        self.assertFalse(check_left_probe_gc_junction("ATCATGCATGATCG"))
        self.assertFalse(check_left_probe_gc_junction("AACATGCATGATCG"))
        self.assertFalse(check_left_probe_gc_junction("ACCATGCATGATCG"))

    def test_at_junction_T_first(self):
        """测试T开头的探针序列（有利于SplintR连接）"""
        self.assertFalse(check_left_probe_gc_junction("TACATGCATGATCG"))
        self.assertFalse(check_left_probe_gc_junction("TTCATGCATGATCG"))
        self.assertFalse(check_left_probe_gc_junction("TGCATGCATGATCG"))
        self.assertFalse(check_left_probe_gc_junction("TCCATGCATGATCG"))

    def test_short_sequence(self):
        """测试序列长度小于1的情况"""
        probe_seq = ""
        self.assertFalse(check_left_probe_gc_junction(probe_seq))

    def test_single_base_G(self):
        """测试单个G碱基（不利于连接）"""
        self.assertTrue(check_left_probe_gc_junction("G"))

    def test_single_base_C(self):
        """测试单个C碱基（不利于连接）"""
        self.assertTrue(check_left_probe_gc_junction("C"))

    def test_single_base_A(self):
        """测试单个A碱基（有利于连接）"""
        self.assertFalse(check_left_probe_gc_junction("A"))

    def test_single_base_T(self):
        """测试单个T碱基（有利于连接）"""
        self.assertFalse(check_left_probe_gc_junction("T"))


class TestCreateBlastDb(unittest.TestCase):

    @patch("scripts.utils.NcbimakeblastdbCommandline")
    def test_non_arabidopsis_lowercase_sequences_allowed(self, mock_makeblastdb):
        mock_makeblastdb.return_value = MagicMock()

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "input.fasta")
            with open(fasta_path, "w") as handle:
                handle.write(
                    ">Mp1g00010.1\natgc\n"  # Liverwort entry with lowercase sequence
                )

            result = create_blast_db(fasta_path)

            self.assertEqual(result, fasta_path)
            mock_makeblastdb.assert_called_once()

            with open(fasta_path, "r") as handle:
                content = handle.read()

        self.assertIn(">Mp1g00010.1\nATGC\n", content)

    @patch("scripts.utils.NcbimakeblastdbCommandline")
    def test_arabidopsis_lowercase_sequences_raise(self, mock_makeblastdb):
        mock_makeblastdb.return_value = MagicMock()

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "input.fasta")
            with open(fasta_path, "w") as handle:
                handle.write(
                    ">AT1G01010\natgc\n"  # Arabidopsis entry should fail due to lowercase
                )

            with self.assertRaisesRegex(
                ValueError,
                "Sequence AT1G01010 contains lowercase characters"
            ):
                create_blast_db(fasta_path)

        mock_makeblastdb.assert_not_called()


if __name__ == "__main__":
    unittest.main()
