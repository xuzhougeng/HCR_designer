"""Subgenome membership comes from a file, never from guessing identifier spelling."""
import os
from pathlib import Path
import tempfile
import unittest

from src.common.reference_annotation import (MAX_LABEL_LENGTH, SubgenomeAnnotation,
                                             format_subgenome_map, load_subgenome_map,
                                             parse_subgenome_map)


class SubgenomeMapTests(unittest.TestCase):
    def test_separators_header_and_comments(self):
        text = ('transcript_id\tsubgenome\n'
                '# a comment\n'
                'BjuVA06G01970\tA\n'
                'BjuVB06G29330,B\n'
                'BjuVA01G00010   A\n'
                '\n')
        mapping = parse_subgenome_map(text)
        annotation = SubgenomeAnnotation(mapping)
        self.assertEqual(len(annotation), 3)
        self.assertEqual(annotation.of('bjuva06g01970'), 'A')
        self.assertEqual(annotation.of('BjuVB06G29330'), 'B')
        self.assertIsNone(annotation.of('BoC01G000510'))
        self.assertEqual(annotation.labels(), ['A', 'B'])

    def test_a_first_row_that_is_not_a_header_is_kept(self):
        annotation = SubgenomeAnnotation(parse_subgenome_map('BjuVA06G01970\tA\nBjuVB06G29330\tB\n'))
        self.assertEqual(annotation.of('BjuVA06G01970'), 'A')
        self.assertEqual(len(annotation), 2)

    def test_malformed_rows_conflicts_and_empty_files_are_rejected_by_line(self):
        for text, needle in [('BjuVA06G01970\n', '第 1 行'),
                             ('A\tB\tC\n', '第 1 行'),
                             ('X\tA\nX\tB\n', '第 2 行'),
                             ('X\t' + 'L' * (MAX_LABEL_LENGTH + 1) + '\n', '第 1 行'),
                             ('# only comments\n', '没有可用记录'),
                             ('', '没有可用记录')]:
            with self.assertRaises(ValueError) as caught:
                parse_subgenome_map(text)
            self.assertIn(needle, str(caught.exception))
        self.assertEqual(len(parse_subgenome_map('X\tA\nx\tA\n')), 1)  # same label twice is fine

    def test_summary_counts_coverage_against_a_reference(self):
        annotation = SubgenomeAnnotation(parse_subgenome_map('t1\tA\nt2\tA\nt3\tB\nt9\tB\n'))
        summary = annotation.summary(['T1', 't2', 't3', 't4'])
        self.assertEqual(summary['subgenome_counts'], {'A': 2, 'B': 1})
        self.assertEqual((summary['annotated'], summary['unannotated'], summary['outside_reference']), (3, 1, 1))
        self.assertEqual(annotation.summary()['annotated'], 4)

    def test_round_trip_and_cache_follows_file_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'subgenomes.tsv'
            path.write_text(format_subgenome_map([('t1', 'A'), ('t2', 'B')]))
            self.assertEqual(load_subgenome_map(path).of('t2'), 'B')
            self.assertIs(load_subgenome_map(path), load_subgenome_map(path))
            # Same size, rewritten at once: a stat-based key cannot tell these apart.
            size, stamp = path.stat().st_size, path.stat().st_mtime_ns
            path.write_text(format_subgenome_map([('t1', 'A'), ('t2', 'C')]))
            os.utime(path, ns=(stamp, stamp))
            self.assertEqual((path.stat().st_size, path.stat().st_mtime_ns), (size, stamp))
            self.assertEqual(load_subgenome_map(path).of('t2'), 'C')


if __name__ == '__main__':
    unittest.main()
