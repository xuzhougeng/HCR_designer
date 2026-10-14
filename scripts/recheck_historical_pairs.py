"""Re-evaluate historical candidates without designing or renumbering probes.

Run from the project root with ``python -m scripts.recheck_historical_pairs``.
An optional independent audit is a regression oracle, never a search input.
"""
import argparse
from collections import defaultdict
import copy
import hashlib
import json
from pathlib import Path

from scripts.utils import create_blast_db
from src.common.paired_specificity import analyze_probe_sequences, assess_pair, METHOD
from src.common.probeset_select import select_probe_sets
from src.common.sequence_utils import check_left_probe_gc_junction
from src.dual_probe.output import ProbeOutputHandler


def compare_independent_audit(records, audit):
    """Measure detection against the exhaustive audit without using it as input."""
    old = {p['key']: p for p in audit['candidates']}
    current = {p['key']: p for p in records}
    if set(old) != set(current):
        raise ValueError('Independent audit must contain exactly the same candidate keys')
    comparison = []
    for bound in range(4):
        expected = {key for key, p in old.items() if p['cross_counts'][bound]}
        observed = {key for key, p in current.items()
                    if p['specificity']['counts']['cross_species'][bound]}
        missing = sorted(expected - observed)
        comparison.append({'max_substitutions_per_arm': bound,
                           'audit_candidates_with_cross_match': len(expected),
                           'blast_candidates_with_cross_match': len(observed),
                           'missing_keys': missing, 'unexpected_keys': sorted(observed - expected),
                           'selected_missing_keys': [key for key in missing if current[key]['keep']]})
    return comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probes', type=Path, default=Path('data/probe_audit_20261014/original_397_pairs.json'))
    parser.add_argument('--reference', type=Path, default=Path('uploads/bol_bju_cdna.fasta'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase1_regression'))
    parser.add_argument('--audit-file', type=Path, help='Independent earlier audit to check exact matches against')
    parser.add_argument('--max-selected', type=int, default=5)
    parser.add_argument('--min-distance', type=int, default=2)
    args = parser.parse_args()
    candidates = json.loads(args.probes.read_text())
    if len({c['key'] for c in candidates}) != len(candidates):
        raise ValueError('Historical candidate keys must be unique')
    create_blast_db(str(args.reference))
    arms = [c[side]['sequence'] for c in candidates for side in ('left', 'right')]
    cache = analyze_probe_sequences(arms, str(args.reference))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'arm_matches.json').write_text(json.dumps(cache, ensure_ascii=False, indent=2))
    handler = ProbeOutputHandler(args.output)
    groups = defaultdict(list)
    for c in candidates:
        p = {key: copy.deepcopy(c[key]) for key in ('id', 'key', 'gene_id', 'symbol', 'purchased',
                                                   'original_selected', 'full_left', 'full_right', 'bridge_id')}
        for side in ('left', 'right'):
            p[side + '_probe'] = copy.deepcopy(c[side])
            p[side + '_probe']['blast_results'] = copy.deepcopy(cache[c[side]['sequence']])
        p['evaluation'] = {
            'quality_score': handler._calculate_probe_set_score({side: cache[c[side]['sequence']] for side in ('left', 'right')}),
            'has_gc_junction': check_left_probe_gc_junction(c['left']['sequence']),
            'gap': c['gap']}
        p['specificity'] = assess_pair(p, c['gene_id'])
        groups[c['symbol']].append(p)
    records = []
    for symbol, probes in groups.items():
        data = {'target_gene_id': probes[0]['gene_id'], 'coordinate_system': '1-based inclusive',
                'specificity_required': True, 'reference': str(args.reference.resolve()), 'probe_sets': probes}
        select_probe_sets(data, args.max_selected, min_distance=args.min_distance)
        records.extend(data['probe_sets'])
        (args.output / f'{symbol}_filtered.json').write_text(json.dumps(data, ensure_ascii=False, indent=2))
    exact = {p['key'] for p in records if p['specificity']['counts']['cross_species'][0]}
    exact_purchased = {p['key'] for p in records if p['purchased'] and p['key'] in exact}
    exact_old_selected = {p['key'] for p in records if p['original_selected'] and p['key'] in exact}
    selected_exact = [p['key'] for p in records if p['keep'] and p['key'] in exact]
    no_target = [p['key'] for p in records if not p['specificity']['on_target_exact']]
    summary = {
        'method': METHOD, 'candidates': len(records), 'purchased_pairs': sum(p['purchased'] for p in records),
        'cross_species_exact_all': len(exact), 'cross_species_exact_purchased': len(exact_purchased),
        'cross_species_exact_original_selected': len(exact_old_selected),
        'cross_species_exact_purchased_keys': sorted(exact_purchased),
        'new_selected_count': sum(p['keep'] for p in records), 'new_selected_cross_exact': selected_exact,
        'missing_exact_target': no_target,
        'nonzero_specificity_scores': sum(p['evaluation']['specificity_score'] > 0 for p in records),
        'gene_summary': {symbol: {'candidates': len(probes), 'selected': sum(p['keep'] for p in probes)}
                         for symbol, probes in groups.items()},
        'reference_sha256': hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        'indel_checked': False, 'experimental_validation': False,
        'note': 'Selections are computational candidates, not a synthesis recommendation.'}
    if args.audit_file:
        audit = json.loads(args.audit_file.read_text())
        comparison = compare_independent_audit(records, audit)
        summary['cross_match_detection_comparison'] = comparison
        summary['exact_match_regression'] = {'missing': comparison[0]['missing_keys'],
                                           'unexpected': comparison[0]['unexpected_keys']}
    (args.output / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if selected_exact or no_target:
        raise AssertionError('Exact-pair regression failed; see saved summary')
    if args.audit_file and any(summary['exact_match_regression'].values()):
        raise AssertionError('Exact matches differ from the independent audit')


if __name__ == '__main__':
    main()
