"""Locate BLAST without changing the user's system PATH."""
import os
from pathlib import Path
import shutil


def blast_executable(name):
    configured = os.environ.get('HCR_BLAST_BIN')
    if configured:
        candidate = Path(configured) / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
        raise FileNotFoundError(f'{name} is not executable in HCR_BLAST_BIN={configured}')
    installed = shutil.which(name)
    if installed:
        return installed
    root = Path(__file__).resolve().parents[2]
    for candidate in sorted((root / '.tools').glob(f'ncbi-blast-*/bin/{name}'), reverse=True):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    raise FileNotFoundError(f'{name} not found; install NCBI BLAST+ or set HCR_BLAST_BIN')
