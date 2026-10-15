"""On-demand, gapped transcript alignments for reviewing annotated polyploid jobs.

This is a display aid, not a replacement for the exhaustive paired-arm scan.
Membership comes from the job's annotation snapshot. Within other annotated
groups, known cross-species hits take priority; BLAST bit score orders each tier.
"""
from functools import lru_cache
import hashlib
from pathlib import Path
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq

from .blast_runtime import blast_executable
from .paired_specificity import infer_species
from .reference_annotation import load_subgenome_map


class AlignmentUnavailable(ValueError):
    pass


def fingerprint(path):
    stat = Path(path).stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino


@lru_cache(maxsize=8)
def _digest(path, stamp):
    digest = hashlib.sha256()
    with open(path, 'rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    if fingerprint(path) != stamp:
        raise AlignmentUnavailable('参考库正在变化，请稍后重试。')
    return digest.hexdigest()


def parse_alignments(output, subjects, query):
    """Validate BLAST coordinates and normalize the query to its saved 5′→3′ strand."""
    hits = []
    for line in output.splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) != 10 or fields[0] not in subjects:
            raise ValueError('Invalid homolog alignment output')
        token, qs, qe, ss, se, score, expect, identity, qa, sa = fields
        qs, qe, ss, se = map(int, (qs, qe, ss, se))
        qa, sa = qa.upper(), sa.upper()
        if qs > qe:
            qs, qe, ss, se = qe, qs, se, ss
            qa, sa = str(Seq(qa).reverse_complement()), str(Seq(sa).reverse_complement())
        subject = subjects[token]
        if (len(qa) != len(sa) or not 1 <= qs <= qe <= len(query)
                or qa.replace('-', '') != query[qs - 1:qe]
                or len(sa.replace('-', '')) != abs(se - ss) + 1
                or not 1 <= min(ss, se) <= max(ss, se) <= subject['subject_length']):
            raise ValueError('Inconsistent homolog alignment coordinates')
        hits.append(dict(subject, query_start=qs, query_end=qe, subject_start=ss, subject_end=se,
                         strand='+' if ss <= se else '-', bit_score=float(score), evalue=float(expect),
                         identity=float(identity), query_coverage=100 * (qe - qs + 1) / len(query),
                         query_aligned=qa, subject_aligned=sa))
    return sorted(hits, key=lambda hit: (hit['relation'] != 'cross_species', -hit['bit_score'],
                                       -hit['query_coverage'], -hit['identity'], hit['subject_id']))


def find_homolog(reference_path, annotation_path, expected_sha256, target_id, query):
    reference_path, annotation_path = str(Path(reference_path).resolve()), str(Path(annotation_path).resolve())
    stamp = fingerprint(reference_path)
    if not expected_sha256 or _digest(reference_path, stamp) != expected_sha256:
        raise AlignmentUnavailable('参考库与本次设计时的内容不一致，无法补充可靠比对；请重新设计。')
    annotation_sha = hashlib.sha256(Path(annotation_path).read_bytes()).hexdigest()
    result = _find_homolog(reference_path, stamp, annotation_path, annotation_sha, target_id, query)
    if fingerprint(reference_path) != stamp or hashlib.sha256(Path(annotation_path).read_bytes()).hexdigest() != annotation_sha:
        raise AlignmentUnavailable('参考或归属文件在比对期间发生变化，请重试。')
    return result


@lru_cache(maxsize=24)
def _find_homolog(reference_path, stamp, annotation_path, annotation_sha, target_id, query):
    annotation = load_subgenome_map(annotation_path)
    target_group = annotation.of(target_id)
    if not target_group:
        raise AlignmentUnavailable('本次任务的归属文件未标注目标转录本，无法确定另一分组。')
    species = infer_species(target_id)
    subjects = {}
    with tempfile.TemporaryDirectory(prefix='hcr-homolog-') as directory:
        query_file, subject_file = Path(directory) / 'query.fa', Path(directory) / 'subjects.fa'
        query_file.write_text('>target\n' + query + '\n')
        with subject_file.open('w') as output:
            for record in SeqIO.parse(reference_path, 'fasta'):
                group = annotation.of(record.id)
                if not group or group == target_group or record.id.casefold() == target_id.casefold():
                    continue
                token = 's' + str(len(subjects))
                other_species = infer_species(record.id)
                subjects[token] = dict(subject_id=record.id, subject_group=group, subject_length=len(record.seq),
                                       relation='cross_species' if species and other_species and species != other_species else 'other_subgenome')
                output.write('>' + token + '\n' + str(record.seq).upper() + '\n')
        common = dict(target_id=target_id, target_group=target_group, searched_sequences=len(subjects),
                      method='BLASTN gapped local alignment; best HSP per subject; cross-species first, then bit score',
                      coordinate_system='1-based inclusive')
        if not subjects:
            return dict(common, status='no_match', message='归属文件中没有可用于比对的其他分组序列。')
        command = [blast_executable('blastn'), '-query', str(query_file), '-subject', str(subject_file),
                   '-task', 'blastn', '-dust', 'no', '-soft_masking', 'false', '-evalue', '1e-5',
                   '-max_hsps', '1', '-max_target_seqs', str(len(subjects)),
                   '-outfmt', '6 sseqid qstart qend sstart send bitscore evalue pident qseq sseq']
        try:
            completed = subprocess.run(command, text=True, capture_output=True, check=True, timeout=90)
        except subprocess.TimeoutExpired as error:
            raise AlignmentUnavailable('同源序列比对超时，请稍后重试。') from error
        except subprocess.CalledProcessError as error:
            raise AlignmentUnavailable('同源序列比对未完成，请检查 BLAST 环境后重试。') from error
        hits = parse_alignments(completed.stdout, subjects, query)
        if not hits:
            return dict(common, status='no_match', message='在其他已标注分组中未检出显著相似的序列（E ≤ 1e-5）。')
        return dict(common, status='completed', **hits[0])
