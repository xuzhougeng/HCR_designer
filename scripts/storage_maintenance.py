"""Inventory and retention for user uploads. Never follows symbolic links."""
import calendar
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import threading


STORAGE_LOCK = threading.RLock()
ACTIVE_STATES = {'queued', 'running'}
RESULT_FOLDERS = {'split_results': '双探针结果', 'triplet_results': '三探针结果'}
BLAST_SUFFIX = re.compile(r'(?:\.\d+)?\.n(?:db|hr|in|og|os|ot|sq|tf|to|al|js|sd|si|nd|ni)$')


def utcnow():
    return datetime.now(timezone.utc)


def three_months_ago(now):
    """Calendar months, clamping the day at month ends (not a fixed 90 days)."""
    month = now.year * 12 + now.month - 1 - 3
    year, month = divmod(month, 12)
    month += 1
    return now.replace(year=year, month=month, day=min(now.day, calendar.monthrange(year, month)[1]))


def read_manifest(path):
    value = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(value, dict):
        raise ValueError('Invalid manifest')
    return value


def tree_info(path):
    """Size, latest mtime and whether every descendant is safe to remove."""
    if path.is_symlink():
        return 0, 0, False
    try:
        stat = path.stat()
    except OSError:
        return 0, 0, False
    if path.is_file():
        return stat.st_size, stat.st_mtime, True
    if not path.is_dir():
        return 0, stat.st_mtime, False
    size, latest, safe = 0, stat.st_mtime, True
    try:
        children = list(path.iterdir())
    except OSError:
        return 0, stat.st_mtime, False
    for child in children:
        child_size, modified, child_safe = tree_info(child)
        size += child_size
        latest, safe = max(latest, modified), safe and child_safe
    return size, latest, safe


class StorageMaintenance:
    def __init__(self, upload_root):
        self.root = Path(upload_root).resolve()
        self.legacy_active = 0
        self.last_cleanup = None

    def _active_references(self):
        references = set()
        uncertain = False
        folder = self.root / 'split_results'
        if folder.is_dir() and not folder.is_symlink():
            for directory in folder.iterdir():
                if not directory.is_dir() or directory.is_symlink():
                    continue
                try:
                    manifest = read_manifest(directory / 'request.json')
                    status = manifest.get('status')
                    if not isinstance(status, str):
                        uncertain = True
                        continue
                    if status in ACTIVE_STATES or (directory / 'running.lock').exists():
                        reference = manifest.get('reference_path')
                        if reference:
                            references.add(Path(reference).resolve())
                except (OSError, ValueError, TypeError):
                    # A damaged manifest cannot establish which reference is in use.
                    uncertain = True
        return references, uncertain

    def inventory(self, now=None):
        with STORAGE_LOCK:
            return self._inventory(now or utcnow())

    def _inventory(self, now):
        if not self.root.is_dir():
            return []
        entries, consumed = [], set()
        children = sorted(self.root.iterdir())
        references, uncertain = self._active_references()

        def add(key, label, kind, paths, auto=False, manifest=None, blocked=None):
            size, latest, safe = 0, 0, True
            for path in paths:
                amount, modified, path_safe = tree_info(path)
                size += amount
                latest, safe = max(latest, modified), safe and path_safe
            modified = datetime.fromtimestamp(latest, timezone.utc)
            if manifest:
                stamp = next((manifest.get(k) for k in ('completed_at', 'failed_at', 'queued_at', 'created_at')
                              if manifest.get(k)), None)
                if stamp:
                    try:
                        modified = datetime.fromisoformat(stamp.replace('Z', '+00:00'))
                        if modified.tzinfo is None:
                            modified = modified.replace(tzinfo=timezone.utc)
                    except (ValueError, TypeError, AttributeError):
                        blocked = '任务时间记录无效，请检查后再清理'
            if not safe:
                blocked = '含符号链接或特殊文件，需在服务器检查'
            status = (manifest or {}).get('status')
            if not isinstance(status, str):
                status = None
            if manifest is not None and status not in {'prepared', 'queued', 'running', 'completed', 'failed', 'interrupted'}:
                blocked = '任务状态记录无效，请检查后再清理'
            if status in ACTIVE_STATES or any(p.is_dir() and (p / 'running.lock').exists() for p in paths):
                blocked = '任务正在排队或运行'
            expired = auto and modified < three_months_ago(now)
            entries.append(dict(key=key, name=label, kind=kind, size=size, modified_at=modified.isoformat(),
                                status=status, expired=expired, blocked=blocked,
                                file_count=len(paths), _paths=paths))

        # A FASTA, its annotation, compressed source and BLAST volumes are one deletion unit.
        for path in children:
            base = path.name[:-3] if path.name.lower().endswith('.gz') else path.name
            if Path(base).suffix.lower() not in ('.fa', '.fasta', '.fna') or path.name in consumed:
                continue
            paths = [p for p in children if p.name in (base, base + '.subgenomes.tsv') or
                     (p.name.lower().endswith('.gz') and p.name[:-3] == base) or
                     (p.name.startswith(base) and BLAST_SUFFIX.fullmatch(p.name[len(base):]))]
            consumed.update(p.name for p in paths)
            blocked = '参考库正在被任务使用' if self.legacy_active or uncertain or any(
                p.resolve() in references for p in paths) else None
            add(base, base, '参考基因组', paths, blocked=blocked)

        for path in children:
            if path.name in consumed:
                continue
            if path.name in RESULT_FOLDERS and path.is_dir() and not path.is_symlink():
                result_children = sorted(path.iterdir())
                grouped = set()
                # Historical triplet ZIPs live next to the task directories.
                for directory in result_children:
                    if not directory.is_dir() or directory.is_symlink():
                        continue
                    paths = [directory]
                    if re.fullmatch('[a-f0-9]{32}', directory.name):
                        paths += [p for p in result_children if re.fullmatch(
                            r'(?:split|triplet)(?:_batch)?_results_' + directory.name + r'\.zip', p.name)]
                    grouped.update(p.name for p in paths)
                    manifest, blocked = None, None
                    manifest_path = directory / 'request.json'
                    if manifest_path.exists():
                        try:
                            manifest = read_manifest(manifest_path)
                        except (OSError, ValueError):
                            blocked = '任务记录损坏，请在服务器检查'
                    if self.legacy_active and path.name == 'triplet_results':
                        blocked = '三探针任务正在运行'
                    add(str(directory.relative_to(self.root)), directory.name, RESULT_FOLDERS[path.name],
                        paths, auto=True, manifest=manifest, blocked=blocked)
                for result in result_children:
                    if result.name not in grouped:
                        add(str(result.relative_to(self.root)), result.name, RESULT_FOLDERS[path.name], [result],
                            auto=result.suffix.lower() == '.zip',
                            blocked='三探针任务正在运行' if self.legacy_active and path.name == 'triplet_results' else None)
            else:
                auto = bool(re.fullmatch(r'probe_recommend_[a-f0-9]{32}\.csv', path.name))
                add(path.name, path.name, '探针推荐结果' if auto else '其他上传文件', [path], auto=auto,
                    blocked='三探针任务正在运行' if self.legacy_active else None)
        return sorted(entries, key=lambda e: e['size'], reverse=True)

    def delete(self, keys=None, expired_only=False, dry_run=False, now=None):
        """Rebuild the inventory under the same lock used when preparing/queuing jobs."""
        with STORAGE_LOCK:
            entries = self._inventory(now or utcnow())
            by_key = {entry['key']: entry for entry in entries}
            if expired_only:
                selected = [entry for entry in entries if entry['expired']]
            else:
                if not isinstance(keys, list) or not keys or any(not isinstance(key, str) for key in keys):
                    raise ValueError('请选择需要删除的文件。')
                if any(key not in by_key for key in keys):
                    raise ValueError('文件已发生变化，请刷新列表后重试。')
                selected = [by_key[key] for key in dict.fromkeys(keys)]
                if any(entry['blocked'] for entry in selected):
                    raise ValueError('所选文件正在使用或需要人工检查，不能删除。')
            report = dict(deleted=[], skipped=[], errors=[], bytes=0, dry_run=dry_run)
            for entry in selected:
                if entry['blocked']:
                    report['skipped'].append(entry['key'])
                    continue
                try:
                    if not dry_run:
                        for path in entry['_paths']:
                            # Paths come exclusively from the current inventory, never from HTTP input.
                            if path.is_symlink() or self.root not in path.resolve().parents:
                                raise ValueError('Unsafe storage path')
                            amount, _, safe = tree_info(path)
                            if not safe:
                                raise ValueError('Storage entry changed during deletion')
                            shutil.rmtree(path) if path.is_dir() else path.unlink()
                            report['bytes'] += amount
                    else:
                        report['bytes'] += entry['size']
                    report['deleted'].append(entry['key'])
                except (OSError, ValueError):
                    report['errors'].append(entry['key'])
            return report

    def cleanup(self, dry_run=False, now=None):
        report = self.delete(expired_only=True, dry_run=dry_run, now=now)
        if not dry_run:
            self.last_cleanup = dict(report, ran_at=(now or utcnow()).isoformat())
        return report


class CleanupScheduler:
    """Single-process service: run once when started, then every 24 hours."""
    def __init__(self, service, logger, interval=24 * 60 * 60):
        self.service, self.logger, self.interval = service, logger, interval
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._run, name='storage-cleanup', daemon=True)

    def start(self):
        self.thread.start()

    def _run(self):
        while not self.stop.is_set():
            try:
                report = self.service.cleanup()
                self.logger.info('Result retention cleanup: %s', report)
            except Exception:
                self.logger.exception('Result retention cleanup failed; will retry next interval')
            if self.stop.wait(self.interval):
                return
