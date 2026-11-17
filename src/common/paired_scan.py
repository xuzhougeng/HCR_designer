"""Run the exhaustive paired-arm scanner (``paired_scan.cpp``) from Python.

The binary is taken from ``HCR_SCAN_BIN`` when set, otherwise it is compiled
on demand into ``.tools/paired_scan/<source hash>/``. A missing compiler,
a failed run or a truncated output raise; they are never reported as "no hits".
"""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

SOURCE = Path(__file__).resolve().with_name('paired_scan.cpp')
PROJECT_ROOT = Path(__file__).resolve().parents[2]
METHOD = 'exhaustive ungapped substitution scan (pigeonhole seeds on the longer arm, full verification of both arms)'
COLUMNS = ['key', 'subject', 'l_start', 'l_end', 'r_start', 'r_end', 'mm_l', 'mm_r',
           'mm_pos_l_5to3', 'mm_pos_r_5to3', 'matched_l', 'matched_r']


class ScanError(RuntimeError):
    """The scan did not complete; callers must not treat this as a clean result."""


def scanner_executable():
    configured = os.environ.get('HCR_SCAN_BIN')
    if configured:
        candidate = Path(configured)
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
        raise FileNotFoundError(f'HCR_SCAN_BIN={configured} is not an executable file')
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()[:16]
    binary = PROJECT_ROOT / '.tools' / 'paired_scan' / digest / 'paired_scan'
    if binary.is_file() and os.access(binary, os.X_OK):
        return str(binary)
    compiler = next((c for c in (os.environ.get('CXX'), 'c++', 'g++', 'clang++') if c and shutil.which(c)), None)
    if not compiler:
        raise FileNotFoundError('No C++ compiler found to build paired_scan; set HCR_SCAN_BIN to a prebuilt binary')
    binary.parent.mkdir(parents=True, exist_ok=True)
    temporary = binary.with_name('paired_scan.partial')
    subprocess.run([compiler, '-std=c++17', '-O2', '-o', str(temporary), str(SOURCE)],
                   check=True, capture_output=True, text=True)
    temporary.replace(binary)
    return str(binary)


def _parse_positions(text):
    return [int(x) for x in text.split(',')] if text else []


def scan_pairs(pairs, reference_path, max_substitutions=3, timeout=None):
    """Find every sense-strand, ungapped paired match for ``pairs``.

    ``pairs`` is an iterable of ``(left_target, right_target, gap)`` where the
    targets are the transcript segments complementary to the arms (5'->3' on the
    transcript) and ``gap`` the number of transcript bases between them.
    Returns ``(hits, run_info)``: ``hits`` maps each distinct pair to a list of
    hit records sorted by subject and start; ``run_info`` records the search scope.
    """
    reference_path = Path(reference_path).resolve()
    if not reference_path.is_file():
        raise FileNotFoundError(f'Reference FASTA not found: {reference_path}')
    distinct = list(dict.fromkeys((str(l), str(r), int(g)) for l, r, g in pairs))
    for left, right, gap in distinct:
        if gap < 0 or set(left + right) - set('ACGT') or min(len(left), len(right)) < max_substitutions + 1:
            raise ValueError('Pair targets must be A/C/G/T, longer than max_substitutions and non-overlapping')
    run_info = {'method': METHOD, 'reference': str(reference_path),
                'max_substitutions_per_arm': max_substitutions, 'strand': 'sense_transcript_only',
                'indels_checked': False, 'ambiguous_reference_bases': 'windows containing non-ACGT bases are skipped',
                'pairs_scanned': len(distinct)}
    if not distinct:
        return {}, dict(run_info, status='completed', sequences=0, bases=0, paired_hits=0, seconds=0.0)
    binary = scanner_executable()
    try:
        workers = int(os.environ.get('HCR_SCAN_WORKERS', '1'))
    except ValueError as error:
        raise ValueError('HCR_SCAN_WORKERS must be an integer between 1 and 32') from error
    if not 1 <= workers <= 32:
        raise ValueError('HCR_SCAN_WORKERS must be an integer between 1 and 32')
    workers = min(workers, len(distinct))
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='hcr-scan-') as directory:
        def scan_partition(part):
            queries = Path(directory) / f'queries-{part}.tsv'
            output = Path(directory) / f'hits-{part}.tsv'
            indices = range(part, len(distinct), workers)
            queries.write_text(''.join(f'p{i}\t{distinct[i][0]}\t{distinct[i][1]}\t{distinct[i][2]}\n' for i in indices))
            try:
                completed = subprocess.run([binary, str(queries), str(reference_path), str(output), str(max_substitutions)],
                                           capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired as error:
                raise ScanError(f'paired_scan exceeded {timeout} s on {reference_path.name}') from error
            if completed.returncode != 0:
                raise ScanError(f'paired_scan failed (exit {completed.returncode}): {completed.stderr.strip()}')
            lines = output.read_text().splitlines()
            if len(lines) < 2 or lines[0] != '\t'.join(COLUMNS) or not lines[-1].startswith('#DONE'):
                raise ScanError('paired_scan output is incomplete; the scan did not finish')
            trailer = dict(item.split('=') for item in lines[-1].split('\t')[1:])
            if int(trailer['queries']) != len(indices):
                raise ScanError('paired_scan did not account for every query')
            return lines[1:-1], trailer
        if workers == 1:
            partitions = [scan_partition(0)]
        else:
            # Partition queries, never the reference: every pair still searches
            # every transcript. Results and failure handling match the serial run.
            with ThreadPoolExecutor(max_workers=workers) as pool:
                partitions = list(pool.map(scan_partition, range(workers)))
    trailer = dict(partitions[0][1])
    for field in ('sequences', 'bases', 'ambiguous_bases'):
        if any(stats[field] != trailer[field] for _, stats in partitions):
            raise ScanError('paired_scan partitions disagree about the reference')
    for field in ('verified_candidates', 'paired_hits', 'peak_rss_kb'):
        trailer[field] = sum(int(stats[field]) for _, stats in partitions)
    hits = defaultdict(list)
    for line in (line for records, _ in partitions for line in records):
        fields = line.split('\t')
        if len(fields) != len(COLUMNS):
            raise ScanError('paired_scan output has a malformed record')
        left, right, gap = distinct[int(fields[0][1:])]
        hits[(left, right, gap)].append({
            'subject_id': fields[1], 'strand': '+', 'gap': gap,
            'left_start': int(fields[2]), 'left_end': int(fields[3]),
            'right_start': int(fields[4]), 'right_end': int(fields[5]),
            'left_mismatches': int(fields[6]), 'right_mismatches': int(fields[7]),
            'left_mismatch_positions_probe_5to3': _parse_positions(fields[8]),
            'right_mismatch_positions_probe_5to3': _parse_positions(fields[9]),
            'left_matched_target': fields[10], 'right_matched_target': fields[11]})
    for records in hits.values():
        records.sort(key=lambda h: (h['subject_id'], h['left_start']))
    run_info.update(status='completed', sequences=int(trailer['sequences']), bases=int(trailer['bases']),
                    ambiguous_bases=int(trailer['ambiguous_bases']), verified_candidates=int(trailer['verified_candidates']),
                    paired_hits=int(trailer['paired_hits']), seconds=round(time.perf_counter() - started, 3),
                    peak_rss_mb=round(int(trailer['peak_rss_kb']) / 1024, 1),
                    peak_rss_aggregation='sum of worker peaks' if workers > 1 else 'single process peak',
                    workers=workers, executable=binary)
    return {pair: hits.get(pair, []) for pair in distinct}, run_info
