"""Full-arm, same-transcript checks for dual probes.

BLAST supplies seeds, not a proof that all substitution matches were found.
We extend ungapped HSPs to the full query and verify against the saved FASTA.
All public coordinates are 1-based inclusive in the sense transcript direction.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import os
import subprocess
import tempfile

from .blast_runtime import blast_executable


METHOD = 'BLAST-seeded, full-arm ungapped verification; substitution search is not exhaustive'
MAX_SUBSTITUTIONS = 3


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def infer_species(identifier):
    identifier = (identifier or '').casefold()
    if identifier.startswith('bju'):
        return 'Bju'
    if identifier.startswith('bo'):
        return 'Bo'
    return None


@lru_cache(maxsize=1)
def _reference_sequences(path, mtime_ns, size):
    sequences = {}
    name, parts = None, []
    def retain():
        if name is not None:
            if name in sequences:
                raise ValueError(f'Duplicate reference ID: {name}')
            sequences[name] = ''.join(parts).upper()
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                retain()
                name, parts = line[1:].split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError('FASTA sequence precedes its header')
                parts.append(line.strip())
        retain()
    if not sequences:
        raise ValueError('Reference FASTA is empty')
    return sequences


def load_reference(path):
    path = Path(path).resolve()
    stat = path.stat()
    return _reference_sequences(str(path), stat.st_mtime_ns, stat.st_size)


def parse_seed_hits(output, queries, reference, max_substitutions=MAX_SUBSTITUTIONS):
    """Recover full-length hits from HSP coordinates, including clipped ends."""
    hits = {key: {} for key in queries}
    for line in output.splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) != 14:
            raise ValueError('Unexpected BLAST output fields')
        qid, subject = fields[:2]
        if qid not in queries or subject not in reference:
            raise ValueError(f'Unknown BLAST query/reference ID: {qid}, {subject}')
        qlen, qstart, qend, sstart, send, length, nident, mismatches, gaps, gapopen = map(int, fields[2:12])
        query = queries[qid]
        if qlen != len(query):
            raise ValueError('BLAST query length disagrees with input')
        # Antisense transcript matches and gapped HSPs are outside this check.
        if sstart > send or qstart > qend or gaps or gapopen:
            continue
        if not (1 <= qstart <= qend <= qlen) or qend - qstart != send - sstart:
            raise ValueError('Inconsistent ungapped HSP coordinates')
        start = sstart - qstart + 1
        end = start + qlen - 1
        if start < 1 or end > len(reference[subject]):
            continue
        matched = reference[subject][start - 1:end]
        if set(matched) - set('ACGT'):
            continue
        differences = [i for i, (a, b) in enumerate(zip(query, matched)) if a != b]
        if len(differences) > max_substitutions:
            continue
        hits[qid][(subject, start, end)] = {
            'subject_id': subject, 'subject_start': start, 'subject_end': end,
            'subject_strand': '+', 'query_start': 1, 'query_end': qlen,
            'query_length': qlen, 'query_coverage': 1.0,
            'query_sequence': query, 'subject_sequence': matched,
            'mismatches': len(differences), 'gaps': 0, 'gap_opens': 0,
            'total_mismatches': len(differences),
            'mismatch_positions_probe_5to3': sorted(qlen - i for i in differences),
            'seed_hsp': {'query_start': qstart, 'query_end': qend,
                         'subject_start': sstart, 'subject_end': send,
                         'alignment_length': length, 'identities': nident,
                         'mismatches': mismatches, 'gaps': gaps, 'gap_opens': gapopen,
                         'query_alignment': fields[12], 'subject_alignment': fields[13]}
        }
    return {qid: sorted(records.values(), key=lambda h: (h['subject_id'], h['subject_start']))
            for qid, records in hits.items()}


def analyze_probe_sequences(probe_sequences, reference_path, num_threads=4):
    """Batch BLAST all distinct targeting arms; errors propagate to the caller."""
    sequences = list(dict.fromkeys(probe_sequences))
    if any(len(s) < 7 or set(s) - set('ACGT') for s in sequences):
        raise ValueError('Targeting probes must contain A/C/G/T only and be at least 7 nt')
    if not sequences:
        return {}
    reference = load_reference(reference_path)
    queries = {f'q{i}': reverse_complement(s) for i, s in enumerate(sequences)}
    with tempfile.TemporaryDirectory(prefix='hcr-paired-') as directory:
        fasta = Path(directory) / 'queries.fasta'
        fasta.write_text(''.join(f'>{key}\n{sequence}\n' for key, sequence in queries.items()))
        command = [blast_executable('blastn'), '-query', str(fasta), '-db', str(reference_path),
                   '-task', 'blastn-short', '-word_size', '7', '-dust', 'no',
                   '-soft_masking', 'false', '-strand', 'plus', '-ungapped',
                   '-evalue', '1000', '-max_target_seqs', str(max(1, len(reference))),
                   '-num_threads', str(max(1, min(num_threads, os.cpu_count() or 1))),
                   '-outfmt', '6 qseqid sseqid qlen qstart qend sstart send length nident mismatch gaps gapopen qseq sseq']
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        hits = parse_seed_hits(result.stdout, queries, reference)
    return {sequence: {
        'status': 'completed', 'method': METHOD,
        'query_orientation': 'reverse_complement_of_targeting_probe',
        'max_substitutions': MAX_SUBSTITUTIONS,
        'mismatch_stats': dict(Counter(h['mismatches'] for h in hits[f'q{i}'])),
        'detailed_matches': hits[f'q{i}']
    } for i, sequence in enumerate(sequences)}


def assess_pair(probe_set, target_gene_id, target_species=None):
    """Join complete arm hits by transcript, orientation and the designed gap."""
    species = target_species or infer_species(target_gene_id)
    result = {'status': 'not_evaluated', 'method': METHOD,
              'target_gene_id': target_gene_id, 'target_species': species,
              'classification_scope': 'species' if species else 'transcript',
              'paired_matches': [], 'counts': {}, 'on_target_exact': False}
    if not target_gene_id:
        result['status'] = 'missing_target_id'
        return result
    arms = [probe_set[name] for name in ('left_probe', 'right_probe')]
    if any(a.get('blast_results', {}).get('status') != 'completed' for a in arms):
        return result
    gap = arms[1]['position']['start'] - arms[0]['position']['end'] - 1
    if gap < 0:
        raise ValueError('Dual probe coordinates must be 1-based, ordered and non-overlapping')
    for arm in arms:
        if arm['position']['end'] - arm['position']['start'] + 1 != len(arm['sequence']):
            raise ValueError('Probe coordinate span differs from its sequence length')
    result['gap'] = gap
    right_index = defaultdict(list)
    def usable(hit):
        return (hit.get('subject_strand') == '+' and hit.get('query_coverage') == 1.0
                and hit.get('gaps') == 0 and hit.get('query_start') == 1
                and hit.get('query_end') == hit.get('query_length'))
    for hit in arms[1]['blast_results']['detailed_matches']:
        if usable(hit):
            right_index[(hit['subject_id'], hit['subject_start'])].append(hit)
    seen = set()
    for left in arms[0]['blast_results']['detailed_matches']:
        if not usable(left):
            continue
        for right in right_index[(left['subject_id'], left['subject_end'] + gap + 1)]:
            subject = left['subject_id']
            identity = (subject, left['subject_start'], right['subject_start'])
            if identity in seen:
                continue
            seen.add(identity)
            other_species = infer_species(subject)
            if subject.casefold() == target_gene_id.casefold():
                relation = 'target'
            elif species and other_species:
                relation = 'same_species_other' if species == other_species else 'cross_species'
            else:
                relation = 'unknown_species'
            result['paired_matches'].append({
                'subject_id': subject, 'relation': relation,
                'subject_species': other_species, 'strand': '+', 'gap': gap,
                'left_start': left['subject_start'], 'left_end': left['subject_end'],
                'right_start': right['subject_start'], 'right_end': right['subject_end'],
                'left_mismatches': left['mismatches'], 'right_mismatches': right['mismatches'],
                'left_mismatch_positions_probe_5to3': left['mismatch_positions_probe_5to3'],
                'right_mismatch_positions_probe_5to3': right['mismatch_positions_probe_5to3']})
    for relation in ('target', 'same_species_other', 'cross_species', 'unknown_species'):
        matches = [h for h in result['paired_matches'] if h['relation'] == relation]
        # Counts are unique transcripts at each cumulative mismatch bound.
        result['counts'][relation] = [len({h['subject_id'] for h in matches
                                         if max(h['left_mismatches'], h['right_mismatches']) <= bound})
                                      for bound in range(MAX_SUBSTITUTIONS + 1)]
    result['on_target_exact'] = result['counts']['target'][0] > 0
    result['status'] = 'completed'
    return result


def specificity_rejection(assessment):
    """Conservative hard exclusions, separate from ranking heuristics."""
    if not assessment or assessment.get('status') != 'completed':
        return 'specificity_not_evaluated'
    if not assessment.get('on_target_exact'):
        return 'no_exact_on_target_pair'
    if assessment['counts']['cross_species'][0]:
        return 'cross_species_exact_pair'
    if assessment['counts']['unknown_species'][0]:
        return 'unclassified_exact_pair'
    return None


def paired_specificity_score(assessment):
    """Ranking weights, not probabilities or predictions of experimental signal."""
    if not assessment or assessment.get('status') != 'completed':
        return None
    weights = {'cross_species': (100, 20, 6, 2),
               'same_species_other': (10, 2, 1, 0.5),
               'unknown_species': (100, 20, 6, 2)}
    score = 0.0
    for relation, values in weights.items():
        counts = assessment['counts'][relation]
        previous = 0
        for cumulative, weight in zip(counts, values):
            score += (cumulative - previous) * weight
            previous = cumulative
    return score
