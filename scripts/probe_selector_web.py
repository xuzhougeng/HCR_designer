"""Web page for ProbeSelector: pick dual-probe pairs by TPM and download the Excel sheet.

Candidates come from a completed job in the job centre or from an uploaded results ZIP.
Either way each gene is a folder named like the gene table's Name column, as in the
original script; only the four files the selector reads are taken from an upload.
"""
from functools import wraps
import io
import json
import math
from pathlib import Path, PurePosixPath
import shutil
import tempfile
from zipfile import BadZipFile, ZipFile

from flask import jsonify, render_template, request, send_file, url_for

from scripts.i18n import register_i18n, translate
from scripts.split_workflow import InputError, SplitWorkflow
from scripts.storage_maintenance import STORAGE_LOCK
from src.dual_probe.probe_selector import (DEFAULT_OUTPUT_NAME, ProbeSelectorError, build_workbook,
                                           collect_gene_results, discover_gene_dirs, format_warning,
                                           parse_tpm_table, tpm_rule_rows, warning, workbook_bytes)

SELECTOR_FILES = ('dual_probe.bed', 'dual_probe_dual.csv', 'dual_probe.csv', 'dual_probe_detailed.txt')
MAX_EXTRACTED_BYTES = 1 << 30   # only the selector's own files count; detailed reports can be large
MAX_ZIP_MEMBERS = 200_000
MAX_SUMMARY_ROWS = 500
MAX_SUMMARY_WARNINGS = 200
XLSX = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

# Problems that only arise from an uploaded archive (the selector's own ones live in WARNINGS).
ZIP_WARNINGS = {
    'zip_root_unnamed': 'ZIP 根目录的结果文件无法确定基因名（缺少唯一的 *_filtered.json），已跳过。',
    'zip_duplicate_folder': 'ZIP 中有多个名为 {name} 的结果文件夹，只使用 {path}。',
}


def _format(item):
    template = ZIP_WARNINGS.get(item['code'])
    if template is not None:
        return translate(template, **item['values'])
    return format_warning(item, translate)


def _extract_results(upload, workdir):
    """Copy each folder's selector files out of the ZIP; returns ({gene name: folder}, warnings)."""
    archive_path = workdir / 'results.zip'
    upload.save(str(archive_path))
    try:
        archive = ZipFile(archive_path)
    except BadZipFile:
        raise InputError('结果文件需为 ZIP 压缩包。', 'results_zip')
    with archive:
        members = archive.infolist()
        if len(members) > MAX_ZIP_MEMBERS:
            raise InputError('ZIP 中文件过多，请只打包需要挑选的基因结果文件夹。', 'results_zip')
        folders, root_names = {}, []
        for info in members:
            if info.is_dir():
                continue
            parts = [part for part in PurePosixPath(info.filename.replace('\\', '/')).parts if part not in ('', '/', '.')]
            if not parts or '..' in parts or parts[0] == '__MACOSX' or parts[-1].startswith('._'):
                continue
            if parts[-1] in SELECTOR_FILES:
                folders.setdefault(tuple(parts[:-1]), {})[parts[-1]] = info
            elif len(parts) == 1 and parts[0].endswith('_filtered.json'):
                root_names.append(parts[0][:-len('_filtered.json')])
        if sum(info.file_size for files in folders.values() for info in files.values()) > MAX_EXTRACTED_BYTES:
            raise InputError('ZIP 中的结果文件解压后超过 1 GB，请分批上传。', 'results_zip')

        gene_dirs, warnings, budget = {}, [], MAX_EXTRACTED_BYTES
        for index, parent in enumerate(sorted(folders)):
            if parent:
                name = parent[-1]
            elif len(root_names) == 1:  # a single-target download keeps its files at the ZIP root
                name = root_names[0]
            else:
                warnings.append(warning('zip_root_unnamed'))
                continue
            if name in gene_dirs:
                warnings.append(warning('zip_duplicate_folder', name=name, path='/'.join(gene_dirs[name][1]) or '/'))
                continue
            target = workdir / f'gene_{index}'
            target.mkdir()
            for filename, info in folders[parent].items():
                with archive.open(info) as source, open(target / filename, 'wb') as destination:
                    while True:
                        chunk = source.read(1 << 20)
                        if not chunk:
                            break
                        budget -= len(chunk)
                        if budget < 0:
                            raise InputError('ZIP 中的结果文件解压后超过 1 GB，请分批上传。', 'results_zip')
                        destination.write(chunk)
            gene_dirs[name] = (str(target), parent)
    return {name: folder for name, (folder, _) in gene_dirs.items()}, warnings


def _gene_table_text(form, files):
    upload = files.get('gene_file')
    text = form.get('gene_text', '').strip()
    if upload and upload.filename:
        if text:
            raise InputError('基因表文件与粘贴内容只能选择一项。', 'gene_file')
        try:
            return upload.read().decode('utf-8-sig')
        except UnicodeError:
            raise InputError('CSV 需使用 UTF-8 编码。', 'gene_file')
    if not text:
        raise InputError('请上传或粘贴基因表。', 'gene_file')
    return text


def _number(value):
    return value if isinstance(value, (int, float)) and math.isfinite(value) else str(value)


def register_probe_selector_routes(app):
    register_i18n(app)

    def workflow():
        # Only job lookup and listing are used, so no sequence collections are needed.
        return SplitWorkflow(app.config['UPLOAD_FOLDER'], {}, {}, {}, {})

    def api_errors(function):
        @wraps(function)
        def decorated(*args, **kwargs):
            try:
                return function(*args, **kwargs)
            except InputError as error:
                return jsonify(error=str(error), field=error.field), error.status
            except ProbeSelectorError as error:
                return jsonify(error=str(error), field='gene_file'), 400
            except Exception:
                app.logger.exception('Probe selection failed')
                return jsonify(error='挑选未完成。请检查基因表和结果文件后重试。'), 500
        return decorated

    def completed_jobs():
        with STORAGE_LOCK:
            jobs = workflow().list_jobs()
        return [dict(id=job['id'], mode=job['mode'], reference=job.get('reference'),
                     submitted_at=job.get('queued_at') or job.get('created_at'),
                     result_url=url_for('split_result', identifier=job['id']),
                     tasks=[dict(name=t['name'], gene_id=t['gene_id'], bridge_id=t['bridge_id']) for t in job['tasks']])
                for job in jobs if job['status'] == 'completed']

    @app.route('/split/selector', methods=['GET', 'POST'], endpoint='probe_selector')
    @api_errors
    def probe_selector():
        if request.method == 'GET':
            config = dict(jobs=completed_jobs(), selected=request.args.get('job', ''),
                          submit_url=url_for('probe_selector'), filename=DEFAULT_OUTPUT_NAME)
            return render_template('probe_selector.html', config=config, rules=tpm_rule_rows())

        tpm_data, warnings = parse_tpm_table(_gene_table_text(request.form, request.files))
        source = request.form.get('source', 'job')
        workdir = Path(tempfile.mkdtemp(prefix='probe-selector-'))
        try:
            if source == 'job':
                identifier = request.form.get('job_id', '').strip()
                if not identifier:
                    raise InputError('请选择一个已完成的双探针任务。', 'job_id')
                with STORAGE_LOCK:  # retention cleanup must not remove the folders mid-read
                    service = workflow()
                    directory = service.job_dir(identifier)
                    if service.status(identifier)['status'] != 'completed':
                        raise InputError('该任务尚未完成，完成后才能挑选。', 'job_id', status=409)
                    output = directory / 'output'
                    gene_dirs = discover_gene_dirs(str(output)) if output.is_dir() else {}
                    gene_results, more = collect_gene_results(gene_dirs, tpm_data)
            elif source == 'zip':
                upload = request.files.get('results_zip')
                if not upload or not upload.filename:
                    raise InputError('请上传双探针结果 ZIP。', 'results_zip')
                gene_dirs, archive_warnings = _extract_results(upload, workdir)
                gene_results, more = collect_gene_results(gene_dirs, tpm_data)
                more = archive_warnings + more
            else:
                raise InputError('请选择探针候选的来源。', 'source')
        finally:
            shutil.rmtree(workdir, ignore_errors=True)

        warnings = [_format(item) for item in warnings + more]
        if not gene_results:
            raise InputError('\n'.join(['没有可输出的基因：基因表的 Name 需与结果文件夹名（任务名称）一致。'] + warnings[:20]),
                             'gene_file')
        summary = dict(genes=len(gene_results), pairs=sum(g['SelectedPairs'] for g in gene_results),
                       probes=sum(len(g['Probes']) for g in gene_results),
                       corrected=sum(g['Corrected'] for g in gene_results),
                       rows=[[g['GeneName'], g['GeneID'], _number(g['TPM']), g['TargetCount'], g['SelectedPairs'], g['Corrected']]
                             for g in gene_results[:MAX_SUMMARY_ROWS]],
                       warnings=warnings[:MAX_SUMMARY_WARNINGS], warning_count=len(warnings))
        response = send_file(io.BytesIO(workbook_bytes(build_workbook(gene_results))), mimetype=XLSX,
                             as_attachment=True, download_name=DEFAULT_OUTPUT_NAME)
        response.headers['X-Selection-Summary'] = json.dumps(summary)  # ASCII-escaped, header safe
        return response
