"""The exhaustive scanner must agree with brute-force enumeration and never hide failures."""
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from src.common import paired_scan
from src.common.paired_scan import scan_pairs, ScanError

DNA = 'ACGT'


def mutate(rng, s, positions):
    s = list(s)
    for p in positions:
        s[p] = rng.choice(DNA.replace(s[p], ''))
    return ''.join(s)


def brute_force(sequences, pairs, k):
    expected = set()
    for sid, seq in sequences.items():
        for left, right, gap in pairs:
            for start in range(len(seq) - len(left) - gap - len(right) + 1):
                rs = start + len(left) + gap
                l, r = seq[start:start + len(left)], seq[rs:rs + len(right)]
                if set(l + r) - set(DNA):
                    continue
                ml, mr = sum(a != b for a, b in zip(l, left)), sum(a != b for a, b in zip(r, right))
                if ml <= k and mr <= k:
                    expected.add(((left, right, gap), sid, start + 1, start + len(left), rs + 1, rs + len(right), ml, mr))
    return expected


class ScannerTests(unittest.TestCase):
    def test_matches_brute_force_including_boundaries_repeats_and_ambiguous_bases(self):
        rng = random.Random(194714)
        randseq = lambda n: ''.join(rng.choice(DNA) for _ in range(n))
        pairs = [(randseq(15 + i % 11), randseq(15 + (i * 3) % 11), i % 3) for i in range(22)]
        pairs.append(('ACGT' * 12 + 'AC', 'TTGACCA' * 3, 5))  # long arms: seeds longer than the prototype cap
        sequences = {'background': randseq(2500)}
        for i, (left, right, gap) in enumerate(pairs):
            cases = [(left, right), (mutate(rng, left, [0, 5, len(left) - 1]), mutate(rng, right, [1, 6, len(right) - 2])),
                     (mutate(rng, left, [0, 3, 6, 9]), right), (left[:4] + 'N' + left[5:], right)]
            for j, (l, r) in enumerate(cases):
                sequences[f'q{i}_case{j}'] = ('' if j == 0 else randseq(11)) + l + randseq(gap) + r + ('' if j == 1 else randseq(13))
            sequences[f'q{i}_repeat'] = left + randseq(gap) + right + randseq(20) + left + randseq(gap) + right
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / 'reference.fa'
            reference.write_text(''.join(f'>{sid} description text\n{seq[:60].lower()}\n{seq[60:]}\n' for sid, seq in sequences.items()))
            for k in (0, 3):
                hits, info = scan_pairs(pairs, reference, k)
                observed = {(pair, h['subject_id'], h['left_start'], h['left_end'], h['right_start'], h['right_end'],
                             h['left_mismatches'], h['right_mismatches']) for pair, records in hits.items() for h in records}
                self.assertEqual(observed, brute_force(sequences, pairs, k))
                self.assertEqual(info['status'], 'completed')
                self.assertEqual(info['sequences'], len(sequences))
                self.assertFalse(info['indels_checked'])
            for pair, records in hits.items():
                for h in records:
                    for side, arm in (('left', pair[0]), ('right', pair[1])):
                        positions = sorted(len(arm) - i for i, (a, b) in enumerate(zip(arm, h[f'{side}_matched_target'])) if a != b)
                        self.assertEqual(h[f'{side}_mismatch_positions_probe_5to3'], positions)

    def test_duplicate_pairs_are_scanned_once_and_absent_pairs_return_empty_lists(self):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / 'reference.fa'
            reference.write_text('>only\n' + 'ACGTTGCA' * 10 + '\n')
            pair = ('ACGTTGCAACGTTGCA', 'ACGTTGCAACGTTGCA', 0)
            hits, info = scan_pairs([pair, pair, ('TTTTTTTTTTTTTTTTT', 'GGGGGGGGGGGGGGGGG', 2)], reference, 1)
            self.assertEqual(info['pairs_scanned'], 2)
            self.assertEqual(len(hits[pair]), 7)  # exact period-8 repeats; shifted windows exceed 1 mismatch
            self.assertEqual(hits[('TTTTTTTTTTTTTTTTT', 'GGGGGGGGGGGGGGGGG', 2)], [])

    def test_failures_raise_instead_of_reporting_no_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / 'reference.fa'
            reference.write_text('>only\nACGTACGTACGTACGTACGTACGT\n')
            pair = ('ACGTACGTACGTACGT', 'ACGTACGTACGTACGT', 0)
            with self.assertRaises(ValueError):
                scan_pairs([('ACGN' * 4, 'ACGT' * 4, 0)], reference, 3)
            with self.assertRaises(FileNotFoundError):
                scan_pairs([pair], Path(directory) / 'missing.fa', 3)
            with patch.dict(os.environ, {'HCR_SCAN_BIN': str(Path(directory) / 'absent')}):
                with self.assertRaises(FileNotFoundError):
                    scan_pairs([pair], reference, 3)
            (Path(directory) / 'empty.fa').write_text('')
            with self.assertRaises(ScanError):
                scan_pairs([pair], Path(directory) / 'empty.fa', 3)
            with patch.object(paired_scan, 'scanner_executable', return_value='/bin/sh'):
                with self.assertRaises(ScanError):
                    scan_pairs([pair], reference, 3)


if __name__ == '__main__':
    unittest.main()
