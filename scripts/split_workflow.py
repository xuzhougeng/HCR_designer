"""Validated inputs and saved results for the guided dual-probe workflow."""
from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import uuid
from zipfile import ZipFile, ZIP_DEFLATED

from scripts.create_split_probe import main as design_split
from scripts.utils import save_sequence
from src.common.paired_specificity import load_reference, infer_species


DEFAULTS = dict(min_length=15, max_length=20, min_gc=40., max_gc=60.,
                min_tm=47., max_tm=53., min_gap=2, max_selected=5, poly_n=4,
                kmer_size=8, min_kmer_count=2, min_complementary_length=5,
                gc_junction_penalty=20.)
LIMITS = dict(min_length=(15, 20), max_length=(15, 50), min_gc=(0, 100),
              max_gc=(0, 100), min_tm=(0, 100), max_tm=(0, 100), min_gap=(0, 10000),
              max_selected=(1, 10), poly_n=(2, 20), kmer_size=(6, 10),
              min_kmer_count=(2, 20), min_complementary_length=(4, 15),
              gc_junction_penalty=(0, 1000))
LABELS = dict(min_length='最小臂长', max_length='最大臂长', min_gc='最低 GC',
              max_gc='最高 GC', min_tm='最低 Tm', max_tm='最高 Tm', min_gap='组合间距',
              max_selected='保留数量', poly_n='连续同碱基阈值', kmer_size='重复片段长度',
              min_kmer_count='重复次数', min_complementary_length='互补长度',
              gc_junction_penalty='连接位点惩罚')


class InputError(ValueError):
    def __init__(self, message, field=None, status=400):
        super().__init__(message)
        self.field = field
        self.status = status


def number(value, kind, low, high, label, field):
    try:
        if isinstance(value, bool):
            raise ValueError()
        result = kind(str(value).strip())
        if not math.isfinite(result) or not low <= result <= high:
            raise ValueError()
        return result
    except (ValueError, TypeError, OverflowError):
        raise InputError(f'{label}需为 {low}–{high} 之间的' + ('整数。' if kind is int else '数值。'), field)


def parse_parameters(values):
    if not isinstance(values, dict):
        raise InputError('请提供完整的设计参数。')
    result = {key: number(values.get(key, default), type(default), *LIMITS[key], LABELS[key], key)
              for key, default in DEFAULTS.items()}
    for low, high in [('min_length', 'max_length'), ('min_gc', 'max_gc'), ('min_tm', 'max_tm')]:
        if result[low] > result[high]:
            raise InputError(f'{LABELS[low]}不能大于{LABELS[high]}。', low)
    return result


def task_name(value):
    value = (value or '').strip()
    if (not value or len(value) > 80 or value in ('.', '..') or
            re.search(r'[\\/:*?"<>|,\x00-\x1f]', value)):
        raise InputError('任务名称需为 1–80 个字符，不可含路径符号、逗号或换行。', 'name')
    return value


def sequence_text(value):
    lines = value.strip().splitlines()
    if lines and lines[0].startswith('>'):
        lines = lines[1:]
    sequence = ''.join(''.join(lines).split()).upper()
    if not sequence or re.search('[^ACGT]', sequence):
        raise InputError('请输入纯 A/T/C/G 序列或单条 FASTA，不能包含多条序列或模糊碱基。', 'seq')
    return sequence


def write_json(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
    temporary.replace(path)


class SplitWorkflow:
    def __init__(self, upload_root, bridges, cdna, cds, aliases):
        self.upload_root = Path(upload_root).resolve()
        self.root = self.upload_root / 'split_results'
        self.bridges = bridges
        self.cdna, self.cds, self.aliases = cdna, cds, aliases

    def job_dir(self, identifier):
        if not re.fullmatch('[a-f0-9]{32}', identifier):
            raise InputError('找不到这个设计任务。', status=404)
        directory = self.root / identifier
        if not (directory / 'request.json').is_file():
            raise InputError('找不到这个设计任务。', status=404)
        return directory

    def _reference(self, directory, form, files):
        existing = form.get('existing_ref_genome', '').strip()
        upload = files.get('ref_genome')
        has_upload = bool(upload and upload.filename)
        mode = form.get('reference_mode') or ('upload' if has_upload else 'existing' if existing else 'none')
        if mode not in ('existing', 'upload', 'none'):
            raise InputError('请选择有效的参考库来源。', 'reference_mode')
        if existing and has_upload:
            raise InputError('已有参考库和上传参考文件只能选择一项。', 'existing_ref_genome')
        if mode == 'none':
            if existing or has_upload:
                raise InputError('参考库选项与提交的文件不一致，请重新选择。', 'reference_mode')
            return None, {}, None
        if mode == 'existing':
            path = (self.upload_root / existing).resolve()
            if (not existing or path.parent != self.upload_root or
                    path.suffix.lower() not in ('.fa', '.fasta', '.fna') or not path.is_file()):
                raise InputError('请选择存在的参考序列库。', 'existing_ref_genome')
            label = path.name
        else:
            if not has_upload:
                raise InputError('请选择参考序列文件。', 'ref_genome')
            label = Path(upload.filename).name
            stem = label[:-3] if label.lower().endswith('.gz') else label
            if Path(stem).suffix.lower() not in ('.fa', '.fasta', '.fna'):
                raise InputError('参考文件需为 .fa、.fasta、.fna 或它们的 .gz 压缩文件。', 'ref_genome')
            path = directory / 'reference.fasta'
            try:
                with path.open('wb') as target:
                    if label.lower().endswith('.gz'):
                        with gzip.GzipFile(fileobj=upload.stream) as source:
                            shutil.copyfileobj(source, target)
                    else:
                        shutil.copyfileobj(upload.stream, target)
            except (OSError, EOFError):
                raise InputError('无法解压参考文件，请检查文件是否完整。', 'ref_genome')
        try:
            references = load_reference(path)
        except (ValueError, IndexError, UnicodeError) as error:
            raise InputError(f'参考 FASTA 格式有误：{error}', 'ref_genome')
        stat = path.stat()
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        info = {'name': label, 'sequence_count': len(references),
                'species_counts': dict(Counter(infer_species(key) or '其他/未知' for key in references)),
                'sha256': digest.hexdigest()}
        fingerprint = {'mtime_ns': stat.st_mtime_ns, 'size': stat.st_size}
        return str(path), references, dict(info, fingerprint=fingerprint)

    @staticmethod
    def _lookup(identifier, sequences, aliases):
        index = {}
        for key in sequences:
            folded = key.casefold()
            if folded in index:
                raise InputError('参考库含仅大小写不同的重复编号，请先消除歧义。')
            index[folded] = key
        identifier = aliases.get(identifier.upper(), identifier)
        found = index.get(identifier.casefold())
        if not found:
            # A locus may resolve to a sole isoform, never to the first fuzzy hit.
            matches = [key for key in sequences if key.casefold().startswith(identifier.casefold() + '.')]
            if len(matches) == 1:
                found = matches[0]
            elif len(matches) > 1:
                raise InputError(f'{identifier} 对应多个转录本，请填写完整编号，例如 {matches[0]}。', 'geneID')
        if not found:
            raise InputError(f'在所选序列库中找不到 {identifier}，请核对完整转录本编号。', 'geneID')
        return found, sequences[found]

    def prepare(self, form, files, mode=None):
        mode = mode or form.get('input_mode', 'single')
        if mode not in ('single', 'batch'):
            raise InputError('请选择单条或批量输入。', 'input_mode')
        identifier = uuid.uuid4().hex
        directory = self.root / identifier
        directory.mkdir(parents=True)
        try:
            reference_path, references, reference_info = self._reference(directory, form, files)
            source = form.get('sequenceType', 'reference')
            if source not in ('reference', 'custom', 'cdna', 'cds'):
                raise InputError('请选择有效的目标序列来源。', 'sequenceType')
            rows = []
            if mode == 'single':
                rows.append((None, form.get('name'), form.get('bp_id'), form.get('geneID', ''), None))
            else:
                upload = files.get('batch_file')
                raw = form.get('batch_text', '').strip()
                if upload and upload.filename:
                    if raw:
                        raise InputError('批量 CSV 文件与粘贴内容只能选择一项。', 'batch_file')
                    try:
                        raw = upload.read().decode('utf-8-sig')
                    except UnicodeError:
                        raise InputError('CSV 需使用 UTF-8 编码。', 'batch_file')
                for line, row in enumerate(csv.reader(io.StringIO(raw.lstrip('\ufeff'))), 1):
                    if not row or not any(cell.strip() for cell in row):
                        continue
                    row = [cell.strip() for cell in row]
                    if not rows and [v.lower() for v in row[:3]] == ['task_name', 'bp_id', 'gene_id']:
                        continue
                    if len(row) not in (3, 4):
                        raise InputError(f'CSV 第 {line} 行需为 3 或 4 列。', 'batch_file')
                    rows.append((line, *row[:3], row[3] if len(row) == 4 and row[3] else None))
                if not rows:
                    raise InputError('请上传或粘贴至少一行批量任务。', 'batch_file')
                source = 'reference' if references else 'cdna'
            tasks, errors, names = [], [], set()
            for line, name, bridge, gene, count in rows:
                try:
                    name = task_name(name)
                    if name.casefold() in names:
                        raise InputError(f'任务名称 {name} 重复，请为每行使用不同名称。', 'name')
                    names.add(name.casefold())
                    bridge = (bridge or '').strip().upper()
                    if bridge not in self.bridges:
                        raise InputError(f'桥接编号 {bridge or "（空）"} 不在当前探针库中。', 'bp_id')
                    gene = (gene or '').strip()
                    if '\n' in gene or '\r' in gene:
                        raise InputError('目标编号不能包含换行。', 'geneID')
                    if not gene and (references or source != 'custom'):
                        raise InputError('请填写目标转录本编号；粘贴序列时也需要它来判断目标匹配。', 'geneID')
                    if source == 'reference':
                        if not references:
                            raise InputError('从参考库提取序列时必须选择参考库。', 'existing_ref_genome')
                        gene, sequence = self._lookup(gene, references, self.aliases)
                    elif source == 'custom':
                        sequence = sequence_text(form.get('seq', ''))
                        if references:
                            gene, _ = self._lookup(gene, references, self.aliases)
                    else:
                        gene, sequence = self._lookup(gene, self.cdna if source == 'cdna' else self.cds, self.aliases)
                        if references:
                            gene, _ = self._lookup(gene, references, self.aliases)
                    sequence = sequence_text(sequence)
                    if len(sequence) < 30:
                        raise InputError('目标序列至少需要 30 nt。', 'seq')
                    count = None if count is None else number(count, int, 1, 10, '该行保留数量', 'max_selected')
                    tasks.append({'name': name, 'gene_id': gene, 'bridge_id': bridge,
                                  'bridge_sequence': self.bridges[bridge], 'sequence': sequence,
                                  'sequence_length': len(sequence), 'max_selected': count})
                except InputError as error:
                    if mode == 'single':
                        raise
                    errors.append(f'第 {line} 行：{error}')
            if errors:
                raise InputError('\n'.join(errors), 'batch_file')
            manifest = {'id': identifier, 'status': 'prepared', 'mode': mode,
                        'created_at': datetime.now(timezone.utc).isoformat(),
                        'reference_path': reference_path, 'reference': reference_info, 'tasks': tasks}
            write_json(directory / 'request.json', manifest)
            return self.preview(manifest)
        except Exception:
            shutil.rmtree(directory)
            raise

    @staticmethod
    def preview(manifest):
        reference = manifest['reference']
        return {'id': manifest['id'], 'mode': manifest['mode'],
                'reference': {k: v for k, v in reference.items() if k != 'fingerprint'} if reference else None,
                'tasks': [{k: v for k, v in task.items() if k != 'sequence'} for task in manifest['tasks']]}

    def run(self, identifier, values):
        parameters = parse_parameters(values)
        directory = self.job_dir(identifier)
        lock = directory / 'running.lock'
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
            os.close(fd)
        except FileExistsError:
            raise InputError('该任务正在计算，请等待当前请求完成，不要重复提交。', status=409)
        manifest = json.loads((directory / 'request.json').read_text())
        try:
            if manifest['status'] == 'completed':
                result = self.result(identifier)
                if result['parameters'] != parameters:
                    raise InputError('这个任务已完成；修改参数后请重新检查输入以创建新任务。', status=409)
                return result
            if manifest['reference_path']:
                stat = Path(manifest['reference_path']).stat()
                if {'mtime_ns': stat.st_mtime_ns, 'size': stat.st_size} != manifest['reference']['fingerprint']:
                    raise InputError('参考库在输入检查后发生变化，请返回第一步重新检查。')
            manifest['status'] = 'running'
            write_json(directory / 'request.json', manifest)
            output = directory / 'output'
            if output.exists():
                shutil.rmtree(output)
            output.mkdir()
            results = []
            for task in manifest['tasks']:
                task_dir = output / task['name']
                task_dir.mkdir()
                count = task['max_selected'] if task['max_selected'] is not None else parameters['max_selected']
                design_split(name=task['name'], sequence=task['sequence'], gene_id=task['gene_id'],
                             min_length=parameters['min_length'], max_length=parameters['max_length'],
                             gc_min=parameters['min_gc'], gc_max=parameters['max_gc'],
                             tm_min=parameters['min_tm'], tm_max=parameters['max_tm'],
                             min_gap=parameters['min_gap'], min_complementary_length=parameters['min_complementary_length'],
                             poly_n=parameters['poly_n'], kmer_size=parameters['kmer_size'],
                             min_kmer_count=parameters['min_kmer_count'], ref_genome=manifest['reference_path'],
                             bridge_probe=task['bridge_sequence'], bridge_probe_id=task['bridge_id'],
                             output_dir=str(task_dir), max_selected=count,
                             gc_junction_penalty=parameters['gc_junction_penalty'])
                save_sequence(task['gene_id'] or task['name'], task['sequence'], str(task_dir / 'target.fasta'))
                data = json.loads((task_dir / (task['name'] + '_filtered.json')).read_text())
                probes = data['probe_sets']
                compact = []
                for probe in probes:
                    item = {k: v for k, v in probe.items() if k not in ('left_probe', 'right_probe')}
                    for side in ('left_probe', 'right_probe'):
                        item[side] = {k: v for k, v in probe[side].items() if k != 'blast_results'}
                    compact.append(item)
                results.append(dict({k: v for k, v in task.items() if k != 'sequence'}, max_selected=count,
                                    candidates=len(probes), selected=sum(p['keep'] for p in probes),
                                    cross_exact=sum(bool(p.get('specificity', {}).get('counts', {}).get('cross_species', [0])[0]) for p in probes),
                                    unassessed=sum(p.get('specificity', {}).get('status') != 'completed' for p in probes),
                                    probe_sets=compact))
            result = dict(self.preview(manifest), status='completed', parameters=parameters, tasks=results,
                          created_at=manifest['created_at'], completed_at=datetime.now(timezone.utc).isoformat(),
                          coordinate_system='1-based inclusive', indel_checked=False,
                          search_exhaustive=False, experimental_validation=False)
            write_json(output / 'run_summary.json', result)
            with ZipFile(directory / 'results.zip', 'w', ZIP_DEFLATED) as archive:
                for path in sorted(output.rglob('*')):
                    if path.is_file():
                        relative = path.relative_to(output)
                        if manifest['mode'] == 'single' and len(relative.parts) > 1:
                            relative = Path(*relative.parts[1:])
                        archive.write(path, str(relative))
            write_json(directory / 'result.json', result)
            manifest['status'] = 'completed'
            write_json(directory / 'request.json', manifest)
            return result
        except Exception:
            if manifest['status'] != 'completed':
                manifest['status'] = 'failed'
                write_json(directory / 'request.json', manifest)
            raise
        finally:
            lock.unlink(missing_ok=True)

    def result(self, identifier):
        path = self.job_dir(identifier) / 'result.json'
        if not path.is_file():
            raise InputError('该任务尚未生成完整结果。', status=409)
        return json.loads(path.read_text())
