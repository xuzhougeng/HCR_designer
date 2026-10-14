"""Re-evaluate historical candidates without designing or renumbering probes.

Run from the project root with ``python -m scripts.recheck_historical_pairs``.
Paired evidence comes from the exhaustive scanner; the optional independent
audit is a regression oracle compared site by site, never a search input.
"""
import argparse
from collections import defaultdict
import copy
import csv
import hashlib
import json
from pathlib import Path

from scripts.utils import create_blast_db
from src.common.paired_specificity import (analyze_probe_sequences, assess_probe_sets, METHOD,
                                           MAX_SUBSTITUTIONS, TargetPolicy)
from src.common.probeset_select import select_probe_sets
from src.common.sequence_utils import check_left_probe_gc_junction
from src.dual_probe.output import ProbeOutputHandler


def site_key(key, hit):
    return (key, hit['subject_id'], hit['left_start'], hit['left_end'], hit['right_start'], hit['right_end'],
            hit['left_mismatches'], hit['right_mismatches'],
            ','.join(map(str, hit['left_mismatch_positions_probe_5to3'])),
            ','.join(map(str, hit['right_mismatch_positions_probe_5to3'])),
            hit['left_matched_target'], hit['right_matched_target'])


def compare_independent_audit(records, audit, audit_sites=None):
    """Measure detection against the exhaustive audit without using it as input."""
    old = {p['key']: p for p in audit['candidates']}
    current = {p['key']: p for p in records}
    if set(old) != set(current):
        raise ValueError('Independent audit must contain exactly the same candidate keys')
    comparison = []
    for bound in range(MAX_SUBSTITUTIONS + 1):
        expected = {key for key, p in old.items() if p['cross_counts'][bound]}
        observed = {key for key, p in current.items()
                    if p['specificity']['counts']['cross_species'][bound]}
        missing = sorted(expected - observed)
        comparison.append({'max_substitutions_per_arm': bound,
                           'audit_candidates_with_cross_match': len(expected),
                           'scan_candidates_with_cross_match': len(observed),
                           'missing_keys': missing, 'unexpected_keys': sorted(observed - expected),
                           'selected_missing_keys': [key for key in missing if current[key]['keep']]})
    sites = None
    if audit_sites is not None:
        observed = {site_key(p['key'], h) for p in records for h in p['specificity']['paired_matches']}
        sites = {'audit_sites': len(audit_sites), 'scan_sites': len(observed),
                 'missing': sorted(audit_sites - observed), 'unexpected': sorted(observed - audit_sites)}
    return comparison, sites


def read_audit_sites(path):
    with open(path) as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        return {(r['pair_key'], r['subject'], int(r['l_start']), int(r['l_end']), int(r['r_start']), int(r['r_end']),
                 int(r['mm_l']), int(r['mm_r']), r['mm_pos_l_5to3'], r['mm_pos_r_5to3'],
                 r['matched_l_target'], r['matched_r_target']) for r in rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probes', type=Path, default=Path('data/probe_audit_20261014/original_397_pairs.json'))
    parser.add_argument('--reference', type=Path, default=Path('uploads/bol_bju_cdna.fasta'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase2_regression'))
    parser.add_argument('--audit-file', type=Path, help='Independent earlier audit to check candidate sets against')
    parser.add_argument('--audit-sites', type=Path, help='Independent per-site TSV (paired_hits_max3.tsv) to compare')
    parser.add_argument('--policy-mode', choices=('species', 'transcript', 'subgenome'), default='species')
    parser.add_argument('--subgenome-map', type=Path, help='Two-column transcript/subgenome file for subgenome mode')
    parser.add_argument('--exclusion-max-mismatches', type=int, default=0)
    parser.add_argument('--max-selected', type=int, default=5)
    parser.add_argument('--min-distance', type=int, default=2)
    args = parser.parse_args()
    if (args.policy_mode == 'subgenome') != bool(args.subgenome_map):
        parser.error('--policy-mode subgenome requires --subgenome-map, and vice versa')
    policy = TargetPolicy(args.policy_mode, exclusion_max_mismatches=args.exclusion_max_mismatches,
                          subgenome_map=args.subgenome_map)
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
        groups[c['symbol']].append(p)
    records, searches = [], {}
    for symbol, probes in groups.items():
        searches[symbol] = assess_probe_sets(probes, probes[0]['gene_id'], args.reference, policy)
        data = {'target_gene_id': probes[0]['gene_id'], 'coordinate_system': '1-based inclusive',
                'specificity_required': True, 'reference': str(args.reference.resolve()),
                'target_policy': policy.to_dict(), 'search': searches[symbol], 'probe_sets': probes}
        select_probe_sets(data, args.max_selected, min_distance=args.min_distance)
        records.extend(data['probe_sets'])
        (args.output / f'{symbol}_filtered.json').write_text(json.dumps(data, ensure_ascii=False, indent=2))
    exact = {p['key'] for p in records if p['specificity']['counts']['cross_species'][0]}
    exact_purchased = {p['key'] for p in records if p['purchased'] and p['key'] in exact}
    exact_old_selected = {p['key'] for p in records if p['original_selected'] and p['key'] in exact}
    selected_exact = [p['key'] for p in records if p['keep'] and p['key'] in exact]
    no_target = [p['key'] for p in records if not p['specificity']['on_target_exact']]
    summary = {
        'method': METHOD, 'target_policy': policy.to_dict(),
        'candidates': len(records), 'purchased_pairs': sum(p['purchased'] for p in records),
        'cross_species_exact_all': len(exact), 'cross_species_exact_purchased': len(exact_purchased),
        'cross_species_exact_original_selected': len(exact_old_selected),
        'cross_species_exact_purchased_keys': sorted(exact_purchased),
        'cross_species_candidates_by_bound': [sum(bool(p['specificity']['counts']['cross_species'][b]) for p in records)
                                             for b in range(MAX_SUBSTITUTIONS + 1)],
        'new_selected_count': sum(p['keep'] for p in records), 'new_selected_cross_exact': selected_exact,
        'new_selected_keys': sorted(p['key'] for p in records if p['keep']),
        'missing_exact_target': no_target,
        'rejection_counts': dict(sorted(__import__('collections').Counter(p['reason_code'] for p in records).items())),
        'nonzero_specificity_scores': sum(p['evaluation']['specificity_score'] > 0 for p in records),
        'gene_summary': {symbol: {'candidates': len(probes), 'selected': sum(p['keep'] for p in probes)}
                         for symbol, probes in groups.items()},
        'search_runs': {symbol: {k: v for k, v in info.items() if k not in ('reference', 'executable')}
                        for symbol, info in searches.items()},
        'total_scan_seconds': round(sum(info['seconds'] for info in searches.values()), 1),
        'peak_scan_rss_mb': max(info['peak_rss_mb'] for info in searches.values()),
        'reference_sha256': hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        'indel_checked': False, 'experimental_validation': False,
        'note': 'Selections are computational candidates, not a synthesis recommendation.'}
    if args.audit_file:
        audit = json.loads(args.audit_file.read_text())
        sites = read_audit_sites(args.audit_sites) if args.audit_sites else None
        comparison, site_comparison = compare_independent_audit(records, audit, sites)
        summary['cross_match_detection_comparison'] = comparison
        summary['site_comparison'] = site_comparison
    (args.output / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if selected_exact or no_target:
        raise AssertionError('Exact-pair regression failed; see saved summary')
    if args.audit_file:
        if any(c['missing_keys'] or c['unexpected_keys'] for c in summary['cross_match_detection_comparison']):
            raise AssertionError('Cross-species candidate sets differ from the independent audit at some bound')
        sites = summary['site_comparison']
        if sites and (sites['missing'] or sites['unexpected']):
            raise AssertionError('Paired sites differ from the independent audit')


if __name__ == '__main__':
    main()
