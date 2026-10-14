import copy
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from src.common.paired_specificity import (
    analyze_probe_sequences, assess_pair, parse_seed_hits, reverse_complement,
    specificity_rejection,
)
from src.common.probeset_select import evaluate_specificity, select_probe_sets
from src.dual_probe import DualProbeConfig, ProbeSet
from scripts.create_split_probe import generate_dual_probe


LEFT = 'ACGTTGCACTAGCTA'
RIGHT = 'TGCAGATCGACTGAC'


def hit(subject, start, length=15, mismatches=0, **overrides):
    result = {'subject_id': subject, 'subject_start': start, 'subject_end': start + length - 1,
              'subject_strand': '+', 'query_start': 1, 'query_end': length,
              'query_length': length, 'query_coverage': 1.0, 'gaps': 0,
              'mismatches': mismatches, 'mismatch_positions_probe_5to3': list(range(1, mismatches + 1))}
    result.update(overrides)
    return result


def blast(hits):
    return {'status': 'completed', 'mismatch_stats': {0: len(hits)}, 'detailed_matches': hits}


def candidate(left_hits=None, right_hits=None, gap=0, start=1, identifier=1):
    return {'id': identifier,
            'left_probe': {'sequence': reverse_complement(LEFT), 'position': {'start': start, 'end': start + 14},
                           'blast_results': blast(left_hits or []), 'tm': 50, 'gc_content': 50},
            'right_probe': {'sequence': reverse_complement(RIGHT), 'position': {'start': start + 15 + gap, 'end': start + 29 + gap},
                            'blast_results': blast(right_hits or []), 'tm': 50, 'gc_content': 50},
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


class PairedSelectionTests(unittest.TestCase):
    def test_pair_requires_same_transcript_strand_and_gap(self):
        left = [hit('BjuTarget', 1), hit('BoWrongGap', 100), hit('BoOnlyLeft', 20),
                hit('BoWrongStrand', 50), hit('BoCorrect', 80)]
        right = [hit('BjuTarget', 18), hit('BoWrongGap', 116), hit('BoOnlyRight', 37),
                 hit('BoWrongStrand', 67, subject_strand='-'), hit('BoCorrect', 97)]
        p = candidate(left, right, gap=2)
        assessment = assess_pair(p, 'BJUTARGET')
        self.assertTrue(assessment['on_target_exact'])
        self.assertEqual({h['subject_id'] for h in assessment['paired_matches']}, {'BjuTarget', 'BoCorrect'})
        self.assertEqual(assessment['counts']['cross_species'][0], 1)

    def test_target_same_species_cross_species_and_repeated_sites_are_separate(self):
        left = [hit('BjuTarget', 1), hit('BjuOther', 10), hit('BoOther', 50), hit('BoOther', 100)]
        right = [hit('BjuTarget', 16), hit('BjuOther', 25), hit('BoOther', 65), hit('BoOther', 115)]
        assessment = assess_pair(candidate(left, right), 'BjuTarget')
        self.assertEqual(assessment['counts']['target'][0], 1)
        self.assertEqual(assessment['counts']['same_species_other'][0], 1)
        self.assertEqual(assessment['counts']['cross_species'][0], 1)
        self.assertEqual(len(assessment['paired_matches']), 4)

    def test_exact_cross_species_is_excluded_before_quality_or_gc(self):
        p = candidate([hit('BjuTarget', 1), hit('BoOther', 50)],
                      [hit('BjuTarget', 16), hit('BoOther', 65)])
        p['evaluation']['has_gc_junction'] = True
        p['specificity'] = assess_pair(p, 'BjuTarget')
        selected = select_probe_sets({'specificity_required': True, 'probe_sets': [p]})['probe_sets'][0]
        self.assertFalse(selected['keep'])
        self.assertEqual(selected['reason_code'], 'cross_species_exact_pair')
        self.assertGreater(selected['evaluation']['specificity_score'], 0)

    def test_missing_evidence_is_not_zero_and_cannot_pass_required_screen(self):
        p = candidate()
        del p['left_probe']['blast_results']
        self.assertIsNone(evaluate_specificity(p))
        p['specificity'] = assess_pair(p, 'BjuTarget')
        select_probe_sets({'specificity_required': True, 'probe_sets': [p]})
        self.assertFalse(p['keep'])
        self.assertEqual(p['reason_code'], 'specificity_not_evaluated')

    def test_exact_target_has_zero_penalty_but_absent_target_is_rejected(self):
        p = candidate([hit('BjuTarget', 1)], [hit('BjuTarget', 16)])
        p['specificity'] = assess_pair(p, 'BjuTarget')
        self.assertEqual(evaluate_specificity(p), 0)
        self.assertIsNone(specificity_rejection(p['specificity']))
        self.assertEqual(specificity_rejection(assess_pair(candidate(), 'BjuTarget')), 'no_exact_on_target_pair')

    def test_unknown_non_target_species_does_not_silently_pass_exact_match(self):
        p = candidate([hit('BjuTarget', 1), hit('Unclassified', 50)],
                      [hit('BjuTarget', 16), hit('Unclassified', 65)])
        self.assertEqual(specificity_rejection(assess_pair(p, 'BjuTarget')), 'unclassified_exact_pair')

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


class PipelineEvidenceTests(unittest.TestCase):
    def test_output_preserves_evidence_scores_ids_and_coordinate_convention(self):
        left, right = reverse_complement(LEFT), reverse_complement(RIGHT)
        probe = ProbeSet((left, 0, 15, 50., 50.), (right, 15, 30, 50., 50.), False)
        cache = {left: blast([hit('BjuTarget', 1), hit('BoOther', 50)]),
                 right: blast([hit('BjuTarget', 16), hit('BoOther', 65)])}
        with tempfile.TemporaryDirectory() as directory:
            config = DualProbeConfig(output_dir=directory, blast_db='reference.fasta', min_gap=2)
            with patch('scripts.create_split_probe.DualProbeDesigner.design_probes', return_value=[probe]), \
                 patch('src.dual_probe.output.analyze_probe_sequences', return_value=cache):
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
            self.assertFalse(selected['keep'])
            self.assertEqual(len((Path(directory) / 'sample_filtered_dual.csv').read_text().splitlines()), 1)


if __name__ == '__main__':
    unittest.main()
