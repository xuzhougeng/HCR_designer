"""Shared FASTA catalog, ploidy configuration and content-based upload deduplication."""
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

from werkzeug.utils import secure_filename

from scripts.storage_maintenance import STORAGE_LOCK
from src.common.reference_annotation import SubgenomeAnnotation, format_subgenome_map, parse_subgenome_map


class ReferenceLibrary:
    def __init__(self, upload_root, catalog_path):
        self.root = Path(upload_root).resolve()
        self.catalog_path = Path(catalog_path)

    def _read(self):
        if not self.catalog_path.exists():
            return {}
        records = json.loads(self.catalog_path.read_text(encoding='utf-8'))
        if not isinstance(records, dict):
            raise ValueError('参考库配置文件损坏，请联系管理员检查。')
        return records

    def _save(self, records):
        self.catalog_path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(dir=self.catalog_path.parent, prefix='.references-', suffix='.tmp')
        try:
            with os.fdopen(fd, 'w', encoding='utf-8') as handle:
                json.dump(records, handle, ensure_ascii=False, indent=2)
            os.replace(temporary, self.catalog_path)
        finally:
            Path(temporary).unlink(missing_ok=True)

    def _files(self):
        return sorted(path for path in self.root.iterdir() if not path.is_symlink() and path.is_file()
                      and path.suffix.lower() in ('.fa', '.fasta', '.fna')) if self.root.is_dir() else []

    @staticmethod
    def _fingerprint(path):
        stat = path.stat()
        return dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns, inode=stat.st_ino)

    def list(self):
        with STORAGE_LOCK:
            records = self._read()
            listing = []
            for path in self._files():
                metadata = records.get(path.name, {})
                current = self._fingerprint(path)
                verified = metadata.get('fingerprint') == current
                annotation = metadata.get('subgenome')
                annotation_path = self.root / (path.name + '.subgenomes.tsv')
                annotation_ready = bool(annotation and annotation.get('reference_fingerprint') == current
                                        and annotation_path.is_file() and not annotation_path.is_symlink())
                listing.append(dict(filename=path.name, display_name=metadata.get('display_name') or path.stem,
                                    polyploid=metadata.get('polyploid'), sha256=metadata.get('sha256') if verified else None,
                                    sequence_count=metadata.get('sequence_count') if verified else None,
                                    size=current['size'], uploaded_at=metadata.get('uploaded_at'),
                                    subgenome={k: v for k, v in annotation.items() if k != 'reference_fingerprint'} if annotation else None,
                                    subgenome_ready=annotation_ready))
            return listing

    def path(self, filename):
        if not isinstance(filename, str) or filename not in {p.name for p in self._files()}:
            raise ValueError('参考库不存在或已被清理，请重新选择。')
        return self.root / filename

    @staticmethod
    def _digest(path):
        digest = hashlib.sha256()
        with path.open('rb') as source:
            for block in iter(lambda: source.read(1024 * 1024), b''):
                digest.update(block)
        return digest.hexdigest()

    @staticmethod
    def _validate(path):
        seen, sequence = set(), False
        with path.open('r', encoding='utf-8-sig') as source:
            for number, line in enumerate(source, 1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if seen and not sequence:
                        raise ValueError(f'第 {number} 行之前有缺少序列内容的 FASTA 条目。')
                    fields = line[1:].split()
                    if not fields:
                        raise ValueError(f'第 {number} 行缺少序列编号。')
                    folded = fields[0].casefold()
                    if folded in seen:
                        raise ValueError(f'第 {number} 行的编号 {fields[0]} 重复（编号不区分大小写）。')
                    seen.add(folded)
                    sequence = False
                else:
                    if not seen:
                        raise ValueError(f'第 {number} 行的序列出现在 FASTA 标题之前。')
                    if re.search(r'[^ACGTURYSWKMBDHVNacgturyswkmbdhvn]', line):
                        raise ValueError(f'第 {number} 行含非核酸字符，请上传核酸 FASTA。')
                    sequence = True
        if not seen or not sequence:
            raise ValueError('FASTA 没有完整的序列记录。')
        return len(seen)

    def _index(self, path, records):
        fingerprint = self._fingerprint(path)
        metadata = records.setdefault(path.name, {})
        if metadata.get('fingerprint') != fingerprint or not metadata.get('sha256'):
            metadata.pop('sequence_count', None)
            metadata.update(sha256=self._digest(path), fingerprint=fingerprint)
        return metadata

    def index(self, filename):
        with STORAGE_LOCK:
            path = self.path(filename)
            records = self._read()
            metadata = self._index(path, records)
            if metadata.get('sequence_count') is None:
                metadata['sequence_count'] = self._validate(path)
            self._save(records)
            return next(record for record in self.list() if record['filename'] == filename)

    @staticmethod
    def _settings(display_name, polyploid):
        display_name = display_name.strip() if isinstance(display_name, str) else ''
        if not display_name or len(display_name) > 120 or any(ord(c) < 32 for c in display_name):
            raise ValueError('参考库名称需为 1–120 个可见字符。')
        if polyploid is not None and not isinstance(polyploid, bool):
            raise ValueError('请选择多倍体、非多倍体或尚未配置。')
        return dict(display_name=display_name, polyploid=polyploid)

    def _prepare_annotation(self, upload, reference):
        if not upload or not upload.filename:
            return None
        name = Path(upload.filename).name
        if Path(name).suffix.lower() not in ('.csv', '.tsv', '.txt', '.map'):
            raise ValueError('亚基因组归属文件需为 CSV、TSV、TXT 或 MAP 两列文件。')
        try:
            mapping = parse_subgenome_map(upload.read().decode('utf-8-sig'))
        except (UnicodeError, ValueError) as error:
            raise ValueError(f'亚基因组归属文件有误：{error}') from error
        with reference.open(encoding='utf-8-sig') as source:
            identifiers = [line[1:].split()[0] for line in source if line.startswith('>') and line[1:].split()]
        summary = SubgenomeAnnotation(mapping).summary(identifiers)
        if not summary['annotated']:
            raise ValueError('亚基因组归属文件中的编号与参考库没有交集，请检查文件。')
        data = format_subgenome_map(mapping.values()).encode('utf-8')
        summary.update(source=name, sha256=hashlib.sha256(data).hexdigest())
        return data, summary

    def _save_with_annotation(self, records, path, annotation):
        if annotation is None:
            self._save(records)
            return
        destination = self.root / (path.name + '.subgenomes.tsv')
        if destination.is_symlink():
            raise ValueError('归属文件路径异常，请联系管理员检查。')
        previous = destination.read_bytes() if destination.exists() else None
        data, summary = annotation
        records[path.name]['subgenome'] = dict(summary, reference_fingerprint=self._fingerprint(path))
        fd, temporary = tempfile.mkstemp(dir=self.catalog_path.parent, prefix='.annotation-', suffix='.tmp')
        try:
            with os.fdopen(fd, 'wb') as target:
                target.write(data)
            os.replace(temporary, destination)
            try:
                self._save(records)
            except Exception:
                if previous is None:
                    destination.unlink()
                else:
                    destination.write_bytes(previous)
                raise
        finally:
            Path(temporary).unlink(missing_ok=True)

    def annotation(self, filename):
        """Return a verified copy of the library annotation for a new task snapshot."""
        with STORAGE_LOCK:
            path = self.path(filename)
            entry = next(entry for entry in self.list() if entry['filename'] == filename)
            if not entry['subgenome_ready']:
                raise ValueError('参考库尚未上传有效的亚基因组归属文件，请先到参考库管理页配置。')
            data = (self.root / (path.name + '.subgenomes.tsv')).read_bytes()
            if hashlib.sha256(data).hexdigest() != entry['subgenome']['sha256']:
                raise ValueError('参考库的亚基因组归属文件已发生变化，请到参考库管理页重新上传。')
            return data, entry['subgenome']

    def configure(self, filename, display_name, polyploid, subgenome_map=None):
        settings = self._settings(display_name, polyploid)
        with STORAGE_LOCK:
            path = self.path(filename)
            annotation = self._prepare_annotation(subgenome_map, path)
            records = self._read()
            records.setdefault(filename, {}).update(settings)
            self.catalog_path.parent.mkdir(parents=True, exist_ok=True)
            self._save_with_annotation(records, path, annotation)
            return next(record for record in self.list() if record['filename'] == filename)

    def upload(self, upload, display_name, polyploid, subgenome_map=None):
        settings = self._settings(display_name, polyploid)
        original = Path(upload.filename or '').name
        compressed = original.lower().endswith('.gz')
        stem = original[:-3] if compressed else original
        if Path(stem).suffix.lower() not in ('.fa', '.fasta', '.fna'):
            raise ValueError('请上传 .fa、.fasta、.fna 或其 .gz 压缩文件。')
        self.root.mkdir(parents=True, exist_ok=True)
        # Stage outside uploads so a concurrent storage cleanup cannot delete an in-flight upload.
        self.catalog_path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(dir=self.catalog_path.parent, prefix='.reference-upload-', suffix='.tmp')
        staged = Path(temporary)
        try:
            source = gzip.GzipFile(fileobj=upload.stream) if compressed else upload.stream
            digest = hashlib.sha256()
            try:
                with os.fdopen(fd, 'wb') as target:
                    for block in iter(lambda: source.read(1024 * 1024), b''):
                        digest.update(block)
                        target.write(block)
            finally:
                if source is not upload.stream:
                    source.close()
            count = self._validate(staged)
            annotation = self._prepare_annotation(subgenome_map, staged)
            sha256 = digest.hexdigest()
            with STORAGE_LOCK:
                records = self._read()
                for path in self._files():
                    metadata = self._index(path, records)
                    if metadata['sha256'] == sha256:
                        metadata['sequence_count'] = count
                        self._save(records)
                        entry = next(record for record in self.list() if record['filename'] == path.name)
                        return dict(reference=entry, duplicate=True)
                filename = (secure_filename(Path(stem).stem) or ('reference_' + sha256[:12])) + Path(stem).suffix.lower()
                path = self.root / filename
                if path.exists():
                    path = self.root / (Path(filename).stem + '_' + sha256[:12] + Path(filename).suffix)
                if path.exists():
                    raise ValueError('参考文件名冲突，请修改上传文件名后重试。')
                staged.replace(path)
                records[path.name] = dict(settings, sha256=sha256, sequence_count=count,
                                          fingerprint=self._fingerprint(path),
                                          uploaded_at=datetime.now(timezone.utc).isoformat())
                try:
                    self._save_with_annotation(records, path, annotation)
                except Exception:
                    path.unlink()
                    raise
                entry = next(record for record in self.list() if record['filename'] == path.name)
                return dict(reference=entry, duplicate=False)
        except (OSError, EOFError, UnicodeError) as error:
            raise ValueError('参考文件无法读取或保存，请检查 gzip 完整性、UTF-8 编码和存储权限。') from error
        finally:
            staged.unlink(missing_ok=True)
