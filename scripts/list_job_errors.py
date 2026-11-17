"""Print the distinct failure messages of saved dual-probe jobs (read-only).

Usage (from the project folder, any Python 3.8+, no dependencies):

    python3 scripts/list_job_errors.py [UPLOAD_FOLDER_OR_split_results]

Defaults to ``uploads/split_results``. Only ``request.json`` and ``progress.json``
are read; nothing is written, moved or deleted. Each line shows how many jobs share
one (status, error_code, stage, error) combination, most common first.
"""
from collections import Counter
import json
from pathlib import Path
import sys


def job_folder(argument):
    path = Path(argument)
    return path / 'split_results' if (path / 'split_results').is_dir() else path


def read_json(path):
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None


def main(argv):
    root = job_folder(argv[1] if len(argv) > 1 else 'uploads/split_results')
    if not root.is_dir():
        print(f'Not a folder: {root}', file=sys.stderr)
        return 2
    statuses, failures, unreadable = Counter(), Counter(), 0
    for request in sorted(root.glob('*/request.json')):
        manifest = read_json(request)
        if not isinstance(manifest, dict):
            unreadable += 1
            continue
        status = manifest.get('status') or '?'
        statuses[status] += 1
        if manifest.get('error'):
            progress = read_json(request.parent / 'progress.json') or {}
            stage = manifest.get('failed_stage') or progress.get('stage') or '-'
            failures[(status, manifest.get('error_code') or '-', stage, manifest['error'])] += 1
    print('Jobs by status:', ', '.join(f'{k}={v}' for k, v in statuses.most_common()) or 'none')
    if unreadable:
        print(f'Unreadable request.json files: {unreadable}')
    print()
    print('count\tstatus\terror_code\tstage\terror')
    for (status, code, stage, error), count in failures.most_common():
        print(f'{count}\t{status}\t{code}\t{stage}\t' + error.replace('\n', ' / '))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
