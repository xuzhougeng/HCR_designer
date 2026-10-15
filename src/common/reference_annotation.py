"""Subgenome membership for reference transcripts.

Two equivalent ways to declare it are supported:

1. One two-column text file: transcript ID, then its subgenome label
   (tab, comma or whitespace separated; an optional header and ``#`` comments
   are skipped).
2. Several FASTA files, one per subgenome. The workflow merges them and writes
   the same two-column file, so everything downstream reads one format only.

Labels come from the file, never from guessing at identifier spelling; a
transcript absent from the file has no subgenome and is counted as unknown.
"""
from collections import Counter
import hashlib
from pathlib import Path
import re

MAX_LABEL_LENGTH = 64
HEADER_FIRST = {'transcript_id', 'transcript', 'id', 'seqid', 'seq_id', 'name', 'gene_id'}
HEADER_SECOND = {'subgenome', 'sub_genome', 'genome', 'group', 'label'}
_SEPARATOR = re.compile(r'[\t,;]|\s+')


def parse_subgenome_map(text):
    """Return ``{folded transcript id: (original id, label)}``; errors name the line."""
    mapping = {}
    for number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = [f for f in _SEPARATOR.split(line) if f]
        if len(fields) != 2:
            raise ValueError(f'第 {number} 行需为“转录本编号 亚基因组”两列')
        identifier, label = fields
        if number == 1 and identifier.casefold() in HEADER_FIRST and label.casefold() in HEADER_SECOND:
            continue
        if len(label) > MAX_LABEL_LENGTH or any(ord(c) < 32 for c in label):
            raise ValueError(f'第 {number} 行的亚基因组名称无效，需为 1–{MAX_LABEL_LENGTH} 个可见字符')
        folded = identifier.casefold()
        if folded in mapping and mapping[folded][1] != label:
            raise ValueError(f'第 {number} 行的 {identifier} 已归入 {mapping[folded][1]}，不能同时归入 {label}')
        mapping[folded] = (identifier, label)
    if not mapping:
        raise ValueError('亚基因组归属文件没有可用记录')
    return mapping


def format_subgenome_map(pairs):
    """Serialize ``(identifier, label)`` pairs in the format this module reads."""
    return 'transcript_id\tsubgenome\n' + ''.join(f'{identifier}\t{label}\n' for identifier, label in pairs)


_CACHE = {}
_CACHE_LIMIT = 4


def load_subgenome_map(path):
    """Load an annotation file, cached by path and a hash of its content.

    Modification time and size are not a safe key: filesystems with coarse
    timestamps give a same-size rewrite within one tick the same stat result.
    """
    path = Path(path).resolve()
    data = path.read_bytes()
    key = (str(path), hashlib.sha256(data).hexdigest())
    annotation = _CACHE.get(key)
    if annotation is None:
        annotation = SubgenomeAnnotation(parse_subgenome_map(data.decode('utf-8-sig')), source=str(path))
        while len(_CACHE) >= _CACHE_LIMIT:
            _CACHE.pop(next(iter(_CACHE)))
        _CACHE[key] = annotation
    return annotation


class SubgenomeAnnotation:
    def __init__(self, mapping, source=None):
        self._mapping = dict(mapping)
        self.source = str(source) if source else None

    def of(self, identifier):
        """The subgenome of a transcript, or ``None`` when it is not annotated."""
        record = self._mapping.get((identifier or '').casefold())
        return record[1] if record else None

    def labels(self):
        return sorted({label for _, label in self._mapping.values()})

    def __len__(self):
        return len(self._mapping)

    def summary(self, reference_ids=None):
        """Per-label counts, restricted to the reference when one is given."""
        if reference_ids is None:
            counts = Counter(label for _, label in self._mapping.values())
            return {'source': self.source, 'annotated': len(self._mapping),
                    'subgenome_counts': dict(sorted(counts.items()))}
        labels = [self.of(key) for key in reference_ids]
        counts = Counter(label for label in labels if label)
        folded = {key.casefold() for key in reference_ids}
        return {'source': self.source, 'annotated': sum(1 for label in labels if label),
                'unannotated': sum(1 for label in labels if not label),
                'outside_reference': sum(1 for key in self._mapping if key not in folded),
                'subgenome_counts': dict(sorted(counts.items()))}
