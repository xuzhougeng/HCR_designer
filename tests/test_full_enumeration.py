"""Independent geometry oracle and regression checks for complete arm search."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
import gzip
import json
from unittest.mock import patch

from Bio.Seq import Seq
from src.dual_probe import DualProbeConfig, DualProbeDesigner, ProbeSet, ProbeOutputHandler
from src.common.sequence_utils import calculate_kmer_count, is_valid_probe, has_low_complexity, has_dimer_issues
from src.common.paired_specificity import TargetPolicy


SEQ = 'ACGTTGCAAGGCTTACCGATGCATTAGCCGGATCCAATTGGCATCGATCGGATTACCGTAGCTAGGCTAACG'


def coordinates(pairs):
    return {(p.left_probe[1], p.left_probe[2], p.right_probe[1], p.right_probe[2]) for p in pairs}


class FullEnumerationTests(unittest.TestCase):
    def config(self, **kwargs):
        values = dict(min_length=15, max_length=18, gc_min=0, gc_max=100,
                      tm_min=0, tm_max=100, poly_n=50, kmer_size=10,
                      min_kmer_count=50, min_complementary_length=30, candidate_mode='exhaustive')
        values.update(kwargs)
        return DualProbeConfig(**values)

    def test_all_lengths_and_last_possible_start_with_nonzero_inner_gap(self):
        for gap in (0, 2):
            config = self.config(inner_gap=gap)
            pairs = DualProbeDesigner(config).design_probes(SEQ)
            expected = {(s, s+l, s+l+gap, s+l+gap+r)
                        for l in range(15, 19) for r in range(15, 19)
                        for s in range(len(SEQ)-l-r-gap+1)}
            self.assertEqual(coordinates(pairs), expected)
            self.assertEqual(len(pairs), len(expected))  # no duplicate geometries
            for p in pairs:
                for arm in (p.left_probe, p.right_probe):
                    self.assertEqual(arm[0], str(Seq(SEQ[arm[1]:arm[2]]).reverse_complement()))

    def test_alternative_lengths_survive_when_first_pair_is_incompatible(self):
        config = self.config()
        compatible = lambda left, right: len(left[0]) == 16 and len(right[0]) == 17
        with patch.object(DualProbeDesigner, '_check_probe_compatibility', side_effect=compatible):
            old = DualProbeDesigner(replace(config, candidate_mode='first_valid')).design_probes(SEQ)
            full = DualProbeDesigner(config).design_probes(SEQ)
        self.assertEqual(old, [])
        self.assertEqual(coordinates(full), {(s, s+16, s+16, s+33) for s in range(len(SEQ)-32)})

    def test_real_physical_filters_match_independent_pair_loop(self):
        config = self.config(gc_min=40, gc_max=65, tm_min=45, tm_max=55,
                             poly_n=4, kmer_size=8, min_kmer_count=2, min_complementary_length=5)
        expected = set(); repeats = calculate_kmer_count(SEQ, 8, 2)
        for start in range(len(SEQ)):
            for left_length in range(15, 19):
                for right_length in range(15, 19):
                    end = start+left_length+right_length
                    if end > len(SEQ): continue
                    segments = [SEQ[start:start+left_length], SEQ[start+left_length:end]]
                    if any(has_low_complexity(s, repeats) or not is_valid_probe(
                        s, min_length=15, max_length=18, gc_min=40, gc_max=65,
                        tm_min=45, tm_max=55, poly_n=4) for s in segments): continue
                    antisense = [str(Seq(s).reverse_complement()) for s in segments]
                    if has_dimer_issues(antisense, 5): continue
                    expected.add((start, start+left_length, start+left_length, end))
        actual = DualProbeDesigner(config).design_probes(SEQ)
        self.assertTrue(expected)
        self.assertEqual(coordinates(actual), expected)
        legacy = DualProbeDesigner(replace(config, candidate_mode='first_valid')).design_probes(SEQ)
        self.assertTrue(coordinates(legacy) <= expected)

    def test_fixed_length_modes_identical_and_configuration_round_trips(self):
        full = self.config(max_length=15)
        legacy = replace(full, candidate_mode='first_valid')
        self.assertEqual(DualProbeDesigner(full).design_probes(SEQ), DualProbeDesigner(legacy).design_probes(SEQ))
        self.assertEqual(DualProbeConfig(**full.to_dict()).candidate_mode, 'exhaustive')
        with self.assertRaises(ValueError): self.config(candidate_mode='invalid')

    def test_CT_map_allows_both_brassica_juncea_subgenomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'map.tsv'
            path.write_text('BjuVA_target\tT\nBjuVB_copy\tT\nBjuV_contig\tT\nBoC_copy\tC\n')
            policy = TargetPolicy(mode='subgenome', subgenome_map=path)
            for subject in ('BjuVB_copy', 'BjuV_contig'):
                self.assertEqual(policy.relation(subject, 'BjuVA_target', 'Bju'), 'allowed')
            self.assertIn(policy.relation('BoC_copy', 'BjuVA_target', 'Bju'), policy.non_allowed_relations())

    def test_compact_arm_evidence_keeps_full_hits_once_and_preserves_score(self):
        left, right = SEQ[:15], SEQ[15:30]
        probe = ProbeSet((left, 0, 15, 50, 50), (right, 15, 30, 50, 50))
        evidence = {'status': 'completed', 'mismatch_stats': {0: 2},
                    'detailed_matches': [{'subject_id': 'BjuTarget'}, {'subject_id': 'BoOther'}]}
        cache = {left: evidence, right: evidence}
        geometry = (str(Seq(left).reverse_complement()), str(Seq(right).reverse_complement()), 0)
        scan = ({geometry: []}, {'status': 'completed', 'max_substitutions_per_arm': 3})
        with tempfile.TemporaryDirectory() as tmp:
            output = ProbeOutputHandler(Path(tmp))
            with patch('src.dual_probe.output.analyze_probe_sequences', return_value=cache), \
                 patch('src.common.paired_specificity.scan_pairs', return_value=scan):
                output.save_probe_sets([probe, probe], 'test', 'dual_probe', 'reference.fa',
                                       target_gene_id='BjuTarget', compact_arm_evidence=True)
            with gzip.open(Path(tmp)/'single_arm_blast.json.gz', 'rt') as handle: saved = json.load(handle)
            self.assertEqual(set(saved), {left, right})
            self.assertEqual(saved[left]['detailed_matches'], evidence['detailed_matches'])
            for p in output.probe_data['probe_sets']:
                self.assertEqual(p['evaluation']['quality_score'], 10)
                arm = p['left_probe']['blast_results']
                self.assertEqual(arm['details_key'], left)
                self.assertEqual(arm['detailed_match_count'], 2)
                self.assertNotIn('detailed_matches', arm)


if __name__ == '__main__':
    unittest.main()
