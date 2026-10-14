import copy
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from src.common.paired_specificity import (
    analyze_probe_sequences, assess_pair, assess_probe_sets, parse_seed_hits, reverse_complement,
    specificity_rejection, paired_specificity_score, TargetPolicy, pair_geometry,
)
from src.common.probeset_select import evaluate_specificity, select_probe_sets
from src.dual_probe import DualProbeConfig, DualProbeDesigner, ProbeSet
from scripts.create_split_probe import generate_dual_probe


LEFT = 'ACGTTGCACTAGCTA'
RIGHT = 'TGCAGATCGACTGAC'
SEARCH = {'status': 'completed', 'method': 'test', 'max_substitutions_per_arm': 3}


def hit(subject, left_start, gap=0, left_mismatches=0, right_mismatches=0, length=15):
    right_start = left_start + length + gap
    return {'subject_id': subject, 'strand': '+', 'gap': gap,
            'left_start': left_start, 'left_end': left_start + length - 1,
            'right_start': right_start, 'right_end': right_start + length - 1,
            'left_mismatches': left_mismatches, 'right_mismatches': right_mismatches,
            'left_mismatch_positions_probe_5to3': list(range(1, left_mismatches + 1)),
            'right_mismatch_positions_probe_5to3': list(range(1, right_mismatches + 1)),
            'left_matched_target': 'N' * length, 'right_matched_target': 'N' * length}


def candidate(gap=0, start=1, identifier=1):
    return {'id': identifier,
            'left_probe': {'sequence': reverse_complement(LEFT), 'position': {'start': start, 'end': start + 14},
                           'blast_results': {'status': 'completed', 'mismatch_stats': {}, 'detailed_matches': []},
                           'tm': 50, 'gc_content': 50},
            'right_probe': {'sequence': reverse_complement(RIGHT), 'position': {'start': start + 15 + gap, 'end': start + 29 + gap},
                            'blast_results': {'status': 'completed', 'mismatch_stats': {}, 'detailed_matches': []},
                            'tm': 50, 'gc_content': 50},
            'evaluation': {'quality_score': 0, 'has_gc_junction': False}}


class SeedVerificationTests(unittest.TestCase):
    def test_clipped_end_is_extended_and_reported_in_probe_orientation(self):
        changed = 'TA' + LEFT[2:]
        reference = {'BjuTarget': 'GGG' + changed}
        line = '\t'.join(map(str, ['q0', 'BjuTarget', 15, 3, 15, 6, 18, 13, 13, 0, 0, 0, LEFT[2:], changed[2:]]))
        hits = parse_seed_hits(line + '\n' + line, {'q0': LEFT}, reference)['q0']
        self.assertEqual(len(hits), 1)
        self.assertEqual((hits[0]['subject_start'], hits[0]['subject_end']), (4, 18))
        self.assertEqual(hits[0]['mismatches'], 2)
        self.assertEqual(hits[0]['mismatch_positions_probe_5to3'], [14, 15])

    def test_reverse_strand_and_out_of_reference_extensions_are_ignored(self):
        reverse = '\t'.join(map(str, ['q0', 'BjuTarget', 15, 1, 15, 15, 1, 15, 15, 0, 0, 0, LEFT, LEFT]))
        clipped = '\t'.join(map(str, ['q0', 'BjuTarget', 15, 3, 15, 1, 13, 13, 13, 0, 0, 0, LEFT[2:], LEFT[2:]]))
        self.assertEqual(parse_seed_hits(reverse + '\n' + clipped, {'q0': LEFT}, {'BjuTarget': LEFT})['q0'], [])

    def test_blast_error_is_not_reported_as_no_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / 'ref.fasta'
            fasta.write_text('>BjuTarget\n' + LEFT + RIGHT + '\n')
            with patch('src.common.paired_specificity.blast_executable', return_value='blastn'), \
                 patch('src.common.paired_specificity.subprocess.run', side_effect=subprocess.CalledProcessError(1, ['blastn'])):
                with self.assertRaises(subprocess.CalledProcessError):
                    analyze_probe_sequences([reverse_complement(LEFT)], fasta)


class PairedAssessmentTests(unittest.TestCase):
    def test_geometry_is_derived_from_coordinates_and_hits_must_match_it(self):
        self.assertEqual(pair_geometry(candidate(gap=2)), (LEFT, RIGHT, 2))
        with self.assertRaises(ValueError):
            assess_pair(candidate(gap=2), 'BjuTarget', [hit('BjuTarget', 1, gap=0)], search=SEARCH)
        bad = candidate()
        bad['right_probe']['position']['start'] = 10
        with self.assertRaises(ValueError):
            pair_geometry(bad)

    def test_relations_counts_and_repeated_sites(self):
        hits = [hit('BjuTarget', 1), hit('BjuOther', 10, left_mismatches=1), hit('BoOther', 50),
                hit('BoOther', 100, right_mismatches=3), hit('Unknown', 7, left_mismatches=2)]
        assessment = assess_pair(candidate(), 'BJUTARGET', hits, search=SEARCH)
        self.assertEqual(assessment['status'], 'completed')
        self.assertTrue(assessment['on_target_exact'])
        self.assertEqual(assessment['counts']['target'], [1, 1, 1, 1])
        self.assertEqual(assessment['counts']['same_species_other'], [0, 1, 1, 1])
        self.assertEqual(assessment['counts']['cross_species'], [1, 1, 1, 1])
        self.assertEqual(assessment['counts']['unknown_species'], [0, 0, 1, 1])
        self.assertEqual(len(assessment['paired_matches']), 5)
        self.assertEqual(assessment['rejection'], 'cross_species_exact_pair')
        self.assertEqual(paired_specificity_score(assessment), 100 + 2 + 6)

    def test_missing_scan_is_not_evaluated_and_cannot_pass_required_screen(self):
        p = candidate()
        self.assertIsNone(evaluate_specificity(p))
        p['specificity'] = assess_pair(p, 'BjuTarget')
        self.assertEqual(p['specificity']['status'], 'not_evaluated')
        p['specificity'] = assess_pair(p, 'BjuTarget', [], search={'status': 'failed'})
        self.assertEqual(p['specificity']['status'], 'not_evaluated')
        select_probe_sets({'specificity_required': True, 'probe_sets': [p]})
        self.assertFalse(p['keep'])
        self.assertEqual(p['reason_code'], 'specificity_not_evaluated')
        self.assertEqual(assess_pair(p, '', [], search=SEARCH)['status'], 'missing_target_id')

    def test_exact_target_has_zero_penalty_but_absent_target_is_rejected(self):
        p = candidate()
        p['specificity'] = assess_pair(p, 'BjuTarget', [hit('BjuTarget', 1)], search=SEARCH)
        self.assertEqual(evaluate_specificity(p), 0)
        self.assertIsNone(specificity_rejection(p['specificity']))
        self.assertEqual(specificity_rejection(assess_pair(candidate(), 'BjuTarget', [], search=SEARCH)),
                         'no_exact_on_target_pair')

    def test_species_mode_ranks_same_species_but_excludes_cross_and_unknown(self):
        same = assess_pair(candidate(), 'BjuTarget', [hit('BjuTarget', 1), hit('BjuOther', 50)], search=SEARCH)
        self.assertIsNone(specificity_rejection(same))
        self.assertEqual(paired_specificity_score(same), 10)
        unknown = assess_pair(candidate(), 'BjuTarget', [hit('BjuTarget', 1), hit('Unclassified', 50)], search=SEARCH)
        self.assertEqual(specificity_rejection(unknown), 'unclassified_exact_pair')

    def test_transcript_and_allowlist_modes_exclude_non_allowed_copies(self):
        hits = [hit('BjuTarget', 1), hit('BjuVB01G0001', 50), hit('BoC01G0002', 80, left_mismatches=1)]
        strict = assess_pair(candidate(), 'BjuTarget', hits, policy={'mode': 'transcript'}, search=SEARCH)
        self.assertEqual(specificity_rejection(strict), 'same_species_exact_pair')
        allowed = assess_pair(candidate(), 'BjuTarget', hits,
                              policy={'mode': 'allowlist', 'allowed_ids': ['bjuvb01g0001']}, search=SEARCH)
        self.assertEqual([h['relation'] for h in allowed['paired_matches']], ['target', 'allowed', 'cross_species'])
        self.assertIsNone(specificity_rejection(allowed))
        self.assertEqual(paired_specificity_score(allowed), 20)
        near = assess_pair(candidate(), 'BjuTarget', hits,
                           policy={'mode': 'allowlist', 'allowed_ids': ['BjuVB01G0001'], 'exclusion_max_mismatches': 1},
                           search=SEARCH)
        self.assertEqual(specificity_rejection(near), 'cross_species_near_pair')
        self.assertEqual(near['policy'], {'mode': 'allowlist', 'allowed_ids': ['BjuVB01G0001'], 'exclusion_max_mismatches': 1})
        for bad in [{'mode': 'copy'}, {'mode': 'allowlist'}, {'mode': 'species', 'allowed_ids': ['x']},
                    {'exclusion_max_mismatches': 4}]:
            with self.assertRaises(ValueError):
                TargetPolicy.from_dict(bad)

    def test_subgenome_mode_allows_the_target_subgenome_and_excludes_the_others(self):
        hits = [hit('BjuTarget', 1), hit('BjuSameSub', 40), hit('BjuOtherSub', 70),
                hit('BjuNoLabel', 90, left_mismatches=1), hit('BoOther', 120)]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'subgenomes.tsv'
            path.write_text('transcript_id\tsubgenome\nbjutarget\tA\nBjuSameSub\tA\nBjuOtherSub\tB\n')
            policy = {'mode': 'subgenome', 'subgenome_map': str(path)}
            assessment = assess_pair(candidate(), 'BjuTarget', hits, policy=policy, search=SEARCH)
            self.assertEqual([h['relation'] for h in assessment['paired_matches']],
                             ['target', 'allowed', 'same_species_other', 'same_species_other', 'cross_species'])
            self.assertEqual([h['subject_subgenome'] for h in assessment['paired_matches']],
                             ['A', 'A', 'B', None, None])
            self.assertEqual(assessment['target_subgenome'], 'A')
            self.assertEqual(assessment['policy']['target_subgenome'], 'A')
            # Cross-species is reported first when both kinds of non-allowed pair exist.
            self.assertEqual(specificity_rejection(assessment), 'cross_species_exact_pair')
            same_species_only = assess_pair(candidate(), 'BjuTarget', hits[:4], policy=policy, search=SEARCH)
            self.assertEqual(specificity_rejection(same_species_only), 'same_species_exact_pair')
            self.assertEqual(same_species_only['counts']['allowed'], [1, 1, 1, 1])
            # An unlabelled transcript is never allowed by membership, but an explicit ID still is.
            extra = assess_pair(candidate(), 'BjuTarget', hits,
                                policy=dict(policy, allowed_ids=['BjuOtherSub', 'BjuNoLabel']), search=SEARCH)
            self.assertEqual([h['relation'] for h in extra['paired_matches']],
                             ['target', 'allowed', 'allowed', 'allowed', 'cross_species'])
            self.assertEqual(specificity_rejection(extra), 'cross_species_exact_pair')
            # One policy object reused across targets resolves each target's own subgenome.
            shared = TargetPolicy.from_dict(policy)
            self.assertEqual(shared.relation('BjuSameSub', 'BjuTarget', 'Bju'), 'allowed')
            self.assertEqual(shared.relation('BjuSameSub', 'BjuOtherSub', 'Bju'), 'same_species_other')

    def test_subgenome_mode_requires_an_annotated_target_and_a_map(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'subgenomes.tsv'
            path.write_text('BjuSameSub\tA\nBjuOtherSub\tB\n')
            with self.assertRaises(ValueError) as caught:
                assess_pair(candidate(), 'BjuTarget', [hit('BjuTarget', 1)],
                            policy={'mode': 'subgenome', 'subgenome_map': str(path)}, search=SEARCH)
            self.assertIn('no subgenome', str(caught.exception))
        for bad in [{'mode': 'subgenome'}, {'mode': 'species', 'subgenome_map': 'x.tsv'},
                    {'mode': 'transcript', 'target_subgenome': 'A'}]:
            with self.assertRaises(ValueError):
                TargetPolicy.from_dict(bad)
        self.assertNotIn('subgenome_map', TargetPolicy().to_dict())

    def test_exact_cross_species_is_excluded_before_quality_or_gc(self):
        p = candidate()
        p['evaluation']['has_gc_junction'] = True
        p['specificity'] = assess_pair(p, 'BjuTarget', [hit('BjuTarget', 1), hit('BoOther', 50)], search=SEARCH)
        selected = select_probe_sets({'specificity_required': True, 'probe_sets': [p]})['probe_sets'][0]
        self.assertFalse(selected['keep'])
        self.assertEqual(selected['reason_code'], 'cross_species_exact_pair')
        self.assertGreater(selected['evaluation']['specificity_score'], 0)

    def test_assess_probe_sets_scans_each_geometry_once(self):
        probes = [candidate(identifier=1), candidate(identifier=2, start=40), candidate(identifier=3, gap=1)]
        scanned = {}
        def fake_scan(pairs, reference, k, timeout=None):
            scanned['pairs'] = list(pairs)
            return ({(LEFT, RIGHT, 0): [hit('BjuTarget', 1)], (LEFT, RIGHT, 1): []}, dict(SEARCH, reference=str(reference)))
        with patch('src.common.paired_specificity.scan_pairs', side_effect=fake_scan):
            info = assess_probe_sets(probes, 'BjuTarget', 'ref.fa')
        self.assertEqual(scanned['pairs'], [(LEFT, RIGHT, 0), (LEFT, RIGHT, 0), (LEFT, RIGHT, 1)])
        self.assertEqual(info['status'], 'completed')
        self.assertEqual([p['specificity']['on_target_exact'] for p in probes], [True, True, False])

    def test_spacing_and_zero_selection_limit(self):
        first = candidate(start=1, identifier=1)
        second = candidate(start=33, identifier=2)  # two free bases between pairs
        for distance, count in [(2, 2), (3, 1)]:
            data = {'probe_sets': copy.deepcopy([first, second])}
            select_probe_sets(data, max_selected=2, min_distance=distance)
            self.assertEqual(sum(p['keep'] for p in data['probe_sets']), count)
        data = {'probe_sets': [first]}
        select_probe_sets(data, max_selected=0)
        self.assertFalse(first['keep'])


class DesignerGeometryTests(unittest.TestCase):
    def test_all_geometric_positions_are_enumerated_without_occupancy(self):
        # 40 nt, two 15 nt arms, gap 0: starts 0..10 give 11 candidate pairs.
        config = DualProbeConfig(min_length=15, max_length=15, gc_min=0, gc_max=100, tm_min=0, tm_max=100,
                                 poly_n=50, kmer_size=10, min_kmer_count=50, min_complementary_length=16)
        designer = DualProbeDesigner(config)
        with patch('src.dual_probe.designer.has_low_complexity', return_value=False), \
             patch('src.dual_probe.designer.is_valid_probe', return_value=True), \
             patch('src.dual_probe.designer.has_dimer_issues', return_value=False):
            sets = designer.design_probes('ACGTTGCACTAGCTATGCAGATCGACTGACAATTGGCCAT')
        self.assertEqual([(s.left_probe[1], s.right_probe[2]) for s in sets], [(i, i + 30) for i in range(11)])
        self.assertEqual(sets[-1].right_probe[2], 40)


class PipelineEvidenceTests(unittest.TestCase):
    def test_output_preserves_evidence_scores_ids_and_coordinate_convention(self):
        left, right = reverse_complement(LEFT), reverse_complement(RIGHT)
        probe = ProbeSet((left, 0, 15, 50., 50.), (right, 15, 30, 50., 50.), False)
        arm_hit = lambda subject, start: {'subject_id': subject, 'subject_start': start, 'subject_end': start + 14,
                                          'subject_strand': '+', 'query_coverage': 1.0, 'gaps': 0, 'mismatches': 0}
        blast = lambda hits: {'status': 'completed', 'mismatch_stats': {0: len(hits)}, 'detailed_matches': hits}
        cache = {left: blast([arm_hit('BjuTarget', 1), arm_hit('BoOther', 50)]),
                 right: blast([arm_hit('BjuTarget', 16), arm_hit('BoOther', 65)])}
        scan = ({(LEFT, RIGHT, 0): [hit('BjuTarget', 1), hit('BoOther', 50)]}, dict(SEARCH, paired_hits=2))
        with tempfile.TemporaryDirectory() as directory:
            config = DualProbeConfig(output_dir=directory, blast_db='reference.fasta', min_gap=2)
            with patch('scripts.create_split_probe.DualProbeDesigner.design_probes', return_value=[probe]), \
                 patch('src.dual_probe.output.analyze_probe_sequences', return_value=cache), \
                 patch('src.common.paired_specificity.scan_pairs', return_value=scan):
                generate_dual_probe(LEFT + RIGHT, 'sample', 'BP_0001', 'A' * 19, config, target_gene_id='BjuTarget')
            detailed = json.loads((Path(directory) / 'dual_probe_detailed.json').read_text())
            filtered = json.loads((Path(directory) / 'sample_filtered.json').read_text())
            raw, selected = detailed['probe_sets'][0], filtered['probe_sets'][0]
            self.assertEqual(raw['left_probe'], selected['left_probe'])
            self.assertEqual(raw['right_probe'], selected['right_probe'])
            self.assertEqual(raw['evaluation']['quality_score'], 10.0)  # two real extra single-arm hits, not rank 0
            self.assertEqual(raw['evaluation']['quality_score'], selected['evaluation']['quality_score'])
            self.assertEqual(selected['left_probe']['position'], {'start': 1, 'end': 15})
            self.assertEqual(selected['right_probe']['position'], {'start': 16, 'end': 30})
            self.assertEqual(selected['reason_code'], 'cross_species_exact_pair')
            self.assertEqual(filtered['selection_parameters']['min_distance'], 2)
            self.assertEqual(filtered['search']['paired_hits'], 2)
            self.assertEqual(filtered['target_policy']['mode'], 'species')
            self.assertFalse(selected['keep'])
            self.assertEqual(len((Path(directory) / 'sample_filtered_dual.csv').read_text().splitlines()), 1)

    def test_scan_failure_aborts_the_design_instead_of_keeping_candidates(self):
        left, right = reverse_complement(LEFT), reverse_complement(RIGHT)
        probe = ProbeSet((left, 0, 15, 50., 50.), (right, 15, 30, 50., 50.), False)
        empty = {'status': 'completed', 'mismatch_stats': {}, 'detailed_matches': []}
        with tempfile.TemporaryDirectory() as directory:
            config = DualProbeConfig(output_dir=directory, blast_db='reference.fasta', min_gap=2)
            with patch('scripts.create_split_probe.DualProbeDesigner.design_probes', return_value=[probe]), \
                 patch('src.dual_probe.output.analyze_probe_sequences', return_value={left: empty, right: empty}), \
                 patch('src.common.paired_specificity.scan_pairs', side_effect=RuntimeError('scan died')):
                with self.assertRaises(RuntimeError):
                    generate_dual_probe(LEFT + RIGHT, 'sample', 'BP_0001', 'A' * 19, config, target_gene_id='BjuTarget')
            self.assertFalse((Path(directory) / 'sample_filtered.json').exists())


if __name__ == '__main__':
    unittest.main()
