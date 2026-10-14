"""Full-arm, same-transcript checks for dual probes.

Paired evidence comes from the exhaustive scanner in ``paired_scan``: every
sense-strand, ungapped site where both arms match with at most
``MAX_SUBSTITUTIONS`` substitutions each, at the designed gap. Single-arm BLAST
results are kept only as ranking evidence; they do not decide exclusions.
All public coordinates are 1-based inclusive in the sense transcript direction.
"""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import os
import subprocess
import tempfile

from .blast_runtime import blast_executable
from .paired_scan import scan_pairs, METHOD as SCAN_METHOD
from .reference_annotation import load_subgenome_map


METHOD = SCAN_METHOD
MAX_SUBSTITUTIONS = 3
ARM_METHOD = 'BLAST-seeded, full-arm ungapped verification; single-arm ranking evidence only, not exhaustive'
RELATIONS = ('target', 'allowed', 'same_species_other', 'cross_species', 'unknown_species')
POLICY_MODES = ('species', 'transcript', 'allowlist', 'subgenome')


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def infer_species(identifier):
    identifier = (identifier or '').casefold()
    if identifier.startswith('bju'):
        return 'Bju'
    if identifier.startswith('bo'):
        return 'Bo'
    return None


class TargetPolicy:
    """Which transcripts may be co-detected, and how close a non-allowed pair may be.

    ``species``: the historical rule; other transcripts of the target species are
    ranked but never excluded. ``transcript``: only the target transcript is
    allowed. ``allowlist``: the target plus ``allowed_ids`` are allowed. ``subgenome``:
    the target, every transcript that ``subgenome_map`` puts in the target's own
    subgenome, and any extra ``allowed_ids`` are allowed; a transcript the file does
    not annotate is never allowed by membership. In the last three modes every other
    transcript, same species or not, is a non-target. A non-allowed pair whose worse
    arm has at most ``exclusion_max_mismatches`` substitutions is excluded.
    """
    def __init__(self, mode='species', allowed_ids=(), exclusion_max_mismatches=0,
                 subgenome_map=None, target_subgenome=None):
        if mode not in POLICY_MODES:
            raise ValueError(f'Unknown target policy mode: {mode}')
        if not 0 <= int(exclusion_max_mismatches) <= MAX_SUBSTITUTIONS:
            raise ValueError('exclusion_max_mismatches must be between 0 and MAX_SUBSTITUTIONS')
        self.mode = mode
        self.allowed_ids = sorted({str(i).strip() for i in allowed_ids if str(i).strip()})
        self._allowed = {i.casefold() for i in self.allowed_ids}
        if mode == 'allowlist' and not self._allowed:
            raise ValueError('allowlist mode needs at least one allowed transcript ID')
        if mode in ('species', 'transcript') and self._allowed:
            raise ValueError('allowed transcript IDs require allowlist or subgenome mode')
        if mode == 'subgenome':
            if not subgenome_map:
                raise ValueError('subgenome mode needs a subgenome annotation file')
        elif subgenome_map or target_subgenome:
            raise ValueError('a subgenome annotation requires subgenome mode')
        self.subgenome_map = str(subgenome_map) if subgenome_map else None
        self.target_subgenome = str(target_subgenome) if target_subgenome else None
        self.exclusion_max_mismatches = int(exclusion_max_mismatches)
        self._resolved = {}

    @classmethod
    def from_dict(cls, data):
        if data is None:
            return cls()
        if isinstance(data, TargetPolicy):
            return data
        return cls(data.get('mode', 'species'), data.get('allowed_ids', ()),
                   data.get('exclusion_max_mismatches', 0), data.get('subgenome_map'),
                   data.get('target_subgenome'))

    def to_dict(self):
        result = {'mode': self.mode, 'allowed_ids': list(self.allowed_ids),
                  'exclusion_max_mismatches': self.exclusion_max_mismatches}
        if self.mode == 'subgenome':
            result['subgenome_map'] = self.subgenome_map
            result['target_subgenome'] = self.target_subgenome
        return result

    @property
    def annotation(self):
        """The loaded subgenome file, or ``None`` outside subgenome mode."""
        return load_subgenome_map(self.subgenome_map) if self.subgenome_map else None

    def resolve_target(self, target_gene_id):
        """The target's own subgenome; a target the file does not annotate is an error."""
        if self.mode != 'subgenome':
            return None
        key = (target_gene_id or '').casefold()
        if key not in self._resolved:
            label = self.target_subgenome or self.annotation.of(target_gene_id)
            if not label:
                raise ValueError(f'{target_gene_id} has no subgenome in {self.subgenome_map}; '
                                 'subgenome mode cannot define its allowed set')
            self._resolved[key] = label
        return self._resolved[key]

    def subgenome_of(self, subject):
        annotation = self.annotation
        return annotation.of(subject) if annotation else None

    def relation(self, subject, target_gene_id, target_species):
        if subject.casefold() == target_gene_id.casefold():
            return 'target'
        if subject.casefold() in self._allowed:
            return 'allowed'
        if self.mode == 'subgenome':
            label = self.subgenome_of(subject)
            if label and label == self.resolve_target(target_gene_id):
                return 'allowed'
        other = infer_species(subject)
        if target_species and other:
            return 'same_species_other' if target_species == other else 'cross_species'
        return 'unknown_species'

    def non_allowed_relations(self):
        if self.mode == 'species':
            return ('cross_species', 'unknown_species')
        return ('cross_species', 'unknown_species', 'same_species_other')


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
    """Recover full-length single-arm hits from HSP coordinates, including clipped ends."""
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
    """Batch BLAST all distinct targeting arms for single-arm ranking evidence; errors propagate."""
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
        'status': 'completed', 'method': ARM_METHOD,
        'query_orientation': 'reverse_complement_of_targeting_probe',
        'max_substitutions': MAX_SUBSTITUTIONS,
        'mismatch_stats': dict(Counter(h['mismatches'] for h in hits[f'q{i}'])),
        'detailed_matches': hits[f'q{i}']
    } for i, sequence in enumerate(sequences)}


def pair_geometry(probe_set):
    """Transcript-strand targets and gap of a dual probe; validates coordinates."""
    arms = [probe_set[name] for name in ('left_probe', 'right_probe')]
    for arm in arms:
        if arm['position']['end'] - arm['position']['start'] + 1 != len(arm['sequence']):
            raise ValueError('Probe coordinate span differs from its sequence length')
    gap = arms[1]['position']['start'] - arms[0]['position']['end'] - 1
    if gap < 0:
        raise ValueError('Dual probe coordinates must be 1-based, ordered and non-overlapping')
    return reverse_complement(arms[0]['sequence']), reverse_complement(arms[1]['sequence']), gap


def assess_pair(probe_set, target_gene_id, hits=None, policy=None, search=None, target_species=None):
    """Classify the exhaustive paired hits of one candidate under a target policy.

    ``hits`` are the scanner records for this pair's geometry; ``None`` means the
    scan was not run, which is reported as ``not_evaluated`` rather than as clean.
    """
    policy = TargetPolicy.from_dict(policy)
    species = target_species or infer_species(target_gene_id)
    result = {'status': 'not_evaluated', 'method': METHOD,
              'target_gene_id': target_gene_id, 'target_species': species,
              'policy': policy.to_dict(),
              'classification_scope': 'species' if species else 'transcript',
              'search': search, 'paired_matches': [], 'counts': {}, 'on_target_exact': False}
    if not target_gene_id:
        result['status'] = 'missing_target_id'
        return result
    result['gap'] = pair_geometry(probe_set)[2]
    if policy.mode == 'subgenome':
        # Resolve per target: one policy object may serve several targets.
        result['target_subgenome'] = policy.resolve_target(target_gene_id)
        result['policy'] = dict(policy.to_dict(), target_subgenome=result['target_subgenome'])
    if hits is None or not search or search.get('status') != 'completed':
        return result
    for hit in hits:
        if hit.get('strand') != '+' or hit.get('gap') != result['gap']:
            raise ValueError('Scanner hit does not match the candidate geometry')
        record = dict(hit, relation=policy.relation(hit['subject_id'], target_gene_id, species),
                      subject_species=infer_species(hit['subject_id']))
        if policy.mode == 'subgenome':
            record['subject_subgenome'] = policy.subgenome_of(hit['subject_id'])
        result['paired_matches'].append(record)
    bounds = range(search.get('max_substitutions_per_arm', MAX_SUBSTITUTIONS) + 1)
    for relation in RELATIONS:
        matches = [h for h in result['paired_matches'] if h['relation'] == relation]
        # Counts are unique transcripts at each cumulative per-arm mismatch bound.
        result['counts'][relation] = [len({h['subject_id'] for h in matches
                                         if max(h['left_mismatches'], h['right_mismatches']) <= bound})
                                      for bound in bounds]
    result['on_target_exact'] = result['counts']['target'][0] > 0
    result['status'] = 'completed'
    result['rejection'] = specificity_rejection(result)
    return result


def assess_probe_sets(probe_sets, target_gene_id, reference_path, policy=None,
                      max_substitutions=MAX_SUBSTITUTIONS, timeout=None):
    """Scan every candidate once and attach ``specificity`` to each; scan failures propagate."""
    policy = TargetPolicy.from_dict(policy)
    geometries = [pair_geometry(p) for p in probe_sets]
    hits, run_info = scan_pairs(geometries, reference_path, max_substitutions, timeout=timeout)
    for probe_set, geometry in zip(probe_sets, geometries):
        probe_set['specificity'] = assess_pair(probe_set, target_gene_id, hits[geometry], policy, run_info)
    return run_info


REJECTION_LABELS = {'cross_species': 'cross_species', 'unknown_species': 'unclassified',
                    'same_species_other': 'same_species'}


def specificity_rejection(assessment):
    """Conservative hard exclusions, separate from ranking heuristics."""
    if not assessment or assessment.get('status') != 'completed':
        return 'specificity_not_evaluated'
    if not assessment.get('on_target_exact'):
        return 'no_exact_on_target_pair'
    policy = TargetPolicy.from_dict(assessment.get('policy'))
    for relation in policy.non_allowed_relations():
        offending = [max(h['left_mismatches'], h['right_mismatches'])
                     for h in assessment['paired_matches'] if h['relation'] == relation]
        if offending and min(offending) <= policy.exclusion_max_mismatches:
            return f'{REJECTION_LABELS[relation]}_{"exact" if min(offending) == 0 else "near"}_pair'
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
        counts = assessment['counts'].get(relation, [])
        previous = 0
        for cumulative, weight in zip(counts, values):
            score += (cumulative - previous) * weight
            previous = cumulative
    return score
