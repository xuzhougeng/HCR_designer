import unittest
from src.common.sequence_utils import (
    calculate_gc_content,
    check_poly_n,
    is_complementary,
    has_hairpin,
    check_complementarity,
    has_dimer_issues
)


class TestSequenceUtils(unittest.TestCase):
    def test_calculate_gc_content(self):
        self.assertEqual(calculate_gc_content("ATGC"), 0.5)
        self.assertEqual(calculate_gc_content("GGCC"), 1.0)
        self.assertEqual(calculate_gc_content("AATT"), 0.0)
        self.assertEqual(calculate_gc_content("GCTA"), 0.5)

    def test_check_poly_n(self):
        self.assertTrue(check_poly_n("ATGC"))
        self.assertFalse(check_poly_n("AAAA"))
        self.assertTrue(check_poly_n("ATAAAG"))
        self.assertFalse(check_poly_n("ATGGGG"))
        # Test with different n values
        self.assertFalse(check_poly_n("AAA", n=3))
        self.assertTrue(check_poly_n("AAA", n=4))

    def test_is_complementary(self):
        self.assertTrue(is_complementary("AT", "TA"))
        self.assertTrue(is_complementary("GC", "CG"))
        self.assertTrue(is_complementary("ATGC", "TACG"))
        self.assertTrue(is_complementary("AAAA", "TTTT"))
        self.assertFalse(is_complementary("ATGC", "GCAT"))
        self.assertFalse(is_complementary("AT", "AT"))
        self.assertFalse(is_complementary("AT", "TG"))
        self.assertFalse(is_complementary("AT", "T"))

    def test_has_hairpin(self):
        self.assertTrue(has_hairpin("ATGCGCAT", min_stem_length=4))  # Perfect hairpin
        self.assertFalse(has_hairpin("ATGC"))  # Too short
        self.assertTrue(has_hairpin("ATATATAT"))  # No complementary stems
        # Test with different min_stem_length
        self.assertTrue(has_hairpin("ATCGCGAT", min_stem_length=3))
        self.assertFalse(has_hairpin("ATCGCGAT", min_stem_length=5))


    def test_has_dimer_issues(self):
        # Test self-complementarity
        self.assertTrue(has_dimer_issues(["ATGCGCAT"]))
        
        # Test complementarity between primers
        self.assertTrue(has_dimer_issues(["ATGC", "GCAT"]))
        
        # Test no issues
        self.assertTrue(has_dimer_issues(["AAAA", "TTTT"]))
        
        # Test empty list
        self.assertFalse(has_dimer_issues([]))
        
        # Test with different min_complementary_length
        self.assertTrue(has_dimer_issues(["ATG", "CAT"], min_complementary_length=3))
        self.assertFalse(has_dimer_issues(["ATG", "CAT"], min_complementary_length=4))


if __name__ == '__main__':
    unittest.main()
