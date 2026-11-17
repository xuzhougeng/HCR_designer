import contextlib
from datetime import datetime, timezone
import gzip
import hashlib
import json
import re
import shutil
import tempfile
import threading
import uuid
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from src.common.blast_runtime import blast_executable
import os

try:
    import fcntl
except ImportError:  # Windows: builds are still serialised within one process
    fcntl = None

def load_fasta(file_path):
    opener = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    with opener(file_path, mode) as fasta_file:
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, 'fasta')}
    
    return sequences

def load_alias(file_path):
    gene_alias = {}
    with gzip.open(file_path, 'r') as f:
        for line in f:
            line = line.decode('latin-1').strip()
            parts = line.strip().split("\t")
            locus_name = parts[0]
            symbol = parts[1]
            symbol = symbol.upper()
            gene_alias[symbol] = locus_name.upper()
    return gene_alias

def save_sequence(gene, sequence, file_path):
    """
    保存序列到FASTA文件
    
    Parameters:
    -----------
    gene : str
        基因ID
    sequence : str
        序列
    file_path : str
        输出文件路径
    """
    with open(file_path, 'w') as f:
        f.write(f">{gene}\n")
        # 每行输出60个碱基
        for i in range(0, len(sequence), 60):
            f.write(f"{sequence[i:i+60]}\n")

def generate_unique_id():
    unique_id = uuid.uuid4().hex
    return unique_id

# Reference FASTA alphabet: IUPAC nucleotide codes in either case, the same set the
# reference library accepts on upload. Case never matters downstream (BLAST and the
# pair scan both ignore it, and the file itself is left untouched so its recorded
# SHA-256 stays valid); windows containing anything but A/C/G/T are skipped there.
REFERENCE_ALPHABET = frozenset('ACGTURYSWKMBDHVN')

# Every file makeblastdb may write for a nucleotide database, single- or multi-volume
# (``ref.fa.nhr`` or ``ref.fa.00.nhr`` + ``ref.fa.nal``). Matches storage_maintenance.
BLAST_INDEX_SUFFIX = re.compile(r'(?:\.\d+)?\.n(?:db|hr|in|og|os|ot|sq|tf|to|al|js|sd|si|nd|ni)$')
# Written last, after a complete build: ties the index files to the FASTA they came from.
BLAST_INDEX_RECORD = '.blastdb.json'
BLAST_INDEX_FORMAT = 1
_BUILD_PREFIX = '.blastdb-build-'
_THREAD_LOCKS = {}
_THREAD_LOCKS_GUARD = threading.Lock()


def _fasta_fingerprint(path):
    stat = os.stat(path)
    return dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns, inode=stat.st_ino)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _index_files(fasta_path):
    """Current BLAST index files next to the FASTA: {name: size}."""
    directory, base = os.path.split(os.path.abspath(fasta_path))
    found = {}
    for name in os.listdir(directory):
        if name.startswith(base) and BLAST_INDEX_SUFFIX.fullmatch(name[len(base):]):
            path = os.path.join(directory, name)
            if os.path.isfile(path) and not os.path.islink(path):
                found[name] = os.path.getsize(path)
    return found


def _read_index_record(fasta_path):
    try:
        with open(fasta_path + BLAST_INDEX_RECORD, encoding='utf-8') as handle:
            record = json.load(handle)
    except (OSError, ValueError):
        return None
    if not isinstance(record, dict) or record.get('format') != BLAST_INDEX_FORMAT \
            or record.get('parse_seqids') is not True or not isinstance(record.get('fasta'), dict) \
            or not isinstance(record.get('files'), dict):
        return None
    return record


def _write_index_record(fasta_path, record):
    target = fasta_path + BLAST_INDEX_RECORD
    fd, temporary = tempfile.mkstemp(dir=os.path.dirname(target), prefix='.' + os.path.basename(target) + '.')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
        os.replace(temporary, target)
    except BaseException:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


def _index_is_current(fasta_path, record, verify_content=False):
    """True when the recorded index still matches the FASTA and is complete on disk.

    Size/mtime/inode unchanged means the content is trusted without hashing. Otherwise,
    with ``verify_content`` (only under the build lock), the FASTA is hashed once and,
    if only the timestamps moved, the record is refreshed instead of rebuilding.
    """
    if record is None or _index_files(fasta_path) != record['files']:
        return False
    fingerprint = _fasta_fingerprint(fasta_path)
    recorded = record['fasta']
    if all(recorded.get(key) == value for key, value in fingerprint.items()):
        return True
    if not verify_content or recorded.get('size') != fingerprint['size'] \
            or recorded.get('sha256') != _sha256(fasta_path):
        return False
    record['fasta'] = dict(fingerprint, sha256=recorded['sha256'])
    _write_index_record(fasta_path, record)
    return True


@contextlib.contextmanager
def _build_lock(fasta_path):
    """Serialise index builds of one FASTA across threads and processes."""
    key = os.path.realpath(fasta_path)
    with _THREAD_LOCKS_GUARD:
        thread_lock = _THREAD_LOCKS.setdefault(key, threading.Lock())
    with thread_lock:
        if fcntl is None:
            yield
            return
        # Lock files live outside the reference folder so they never show up as uploads.
        lock_dir = os.path.join(tempfile.gettempdir(), 'hcr-blastdb-locks')
        os.makedirs(lock_dir, exist_ok=True)
        lock_path = os.path.join(lock_dir, hashlib.sha256(key.encode('utf-8')).hexdigest()[:32] + '.lock')
        with open(lock_path, 'a') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)


def _validate_reference(fasta_path):
    sequences = load_fasta(fasta_path)
    if not sequences:
        raise ValueError("FASTA file contains no sequences")
    for seq_id, seq in sequences.items():
        if not seq:
            raise ValueError(f"Sequence {seq_id} is empty")
        # Check if sequence contains only nucleotide codes (either case)
        invalid = set(seq.upper()) - REFERENCE_ALPHABET
        if invalid:
            raise ValueError(f"Sequence {seq_id} contains invalid DNA characters: {''.join(sorted(invalid))}")


def _build_blast_db(fasta_path):
    directory, base = os.path.split(fasta_path)
    # Leftovers of an interrupted build (we hold the lock, so none is in progress).
    for name in os.listdir(directory):
        if name.startswith('.' + base + _BUILD_PREFIX):
            shutil.rmtree(os.path.join(directory, name), ignore_errors=True)

    before = _fasta_fingerprint(fasta_path)
    _validate_reference(fasta_path)
    sha256 = _sha256(fasta_path)

    # Same basename in a sibling folder: multi-volume .nal/.ndb files name their
    # volumes by basename, so the files stay valid after moving next to the FASTA.
    build_dir = tempfile.mkdtemp(dir=directory, prefix='.' + base + _BUILD_PREFIX)
    try:
        cline = NcbimakeblastdbCommandline(
            cmd=blast_executable('makeblastdb'),
            dbtype="nucl",
            input_file=fasta_path,
            out=os.path.join(build_dir, base),
            parse_seqids=True
        )
        cline()
        built = {name: os.path.getsize(os.path.join(build_dir, name)) for name in os.listdir(build_dir)
                 if name.startswith(base) and BLAST_INDEX_SUFFIX.fullmatch(name[len(base):])}
        if not ({base + '.nin', base + '.nhr', base + '.nsq'} <= built.keys() or base + '.nal' in built):
            raise RuntimeError(f"makeblastdb did not produce a complete index for {fasta_path}")
        if _fasta_fingerprint(fasta_path) != before:
            raise RuntimeError(f"{fasta_path} changed while its BLAST index was being built; retry")

        # The record goes first and comes back last, so a crash in between leaves an
        # index without a record, which the next call rebuilds instead of trusting.
        record_path = fasta_path + BLAST_INDEX_RECORD
        if os.path.lexists(record_path):
            os.unlink(record_path)
        for name in _index_files(fasta_path):
            os.unlink(os.path.join(directory, name))
        for name in built:
            os.replace(os.path.join(build_dir, name), os.path.join(directory, name))
        _write_index_record(fasta_path, dict(
            format=BLAST_INDEX_FORMAT, parse_seqids=True, files=built,
            fasta=dict(before, sha256=sha256),
            built_at=datetime.now(timezone.utc).isoformat()))
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)


def create_blast_db(fasta_path):
    """Make sure an up-to-date ``makeblastdb -parse_seqids`` index sits next to the FASTA.

    The index is trusted only when its ``<fasta>.blastdb.json`` record matches the FASTA
    (size/mtime/inode, or SHA-256 when those moved) and the index files on disk are exactly
    the recorded ones. Otherwise it is rebuilt: indexes from older versions (no record,
    possibly without -parse_seqids), indexes of a replaced or edited FASTA, partial indexes
    of an interrupted build and multi-volume indexes are all handled. Builds go to a
    temporary folder first and are serialised per FASTA. The FASTA itself is never modified.

    Sequences may be upper or lower case and may contain IUPAC ambiguity codes
    such as N; anything else (for example protein letters) is rejected.

    Parameters:
    -----------
    fasta_path : str
        Path to FASTA file

    Returns:
    --------
    str: Path to BLAST database

    Raises:
    -------
    ValueError: If sequence validation fails
    """
    path = os.path.abspath(fasta_path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")
    if _index_is_current(path, _read_index_record(path)):
        return fasta_path
    with _build_lock(path):
        # Another thread or process may have finished the same build while we waited.
        if not _index_is_current(path, _read_index_record(path), verify_content=True):
            _build_blast_db(path)
    return fasta_path
