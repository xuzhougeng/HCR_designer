"""HTTP adapters for the shared single/batch guided interface, background jobs and the job centre."""
from concurrent.futures import ThreadPoolExecutor
from functools import wraps
import io
import threading

from flask import jsonify, render_template, request, send_file, url_for
from scripts.split_workflow import SplitWorkflow, InputError, DEFAULTS
from scripts.storage_maintenance import STORAGE_LOCK


def storage_locked(function):
    @wraps(function)
    def decorated(*args, **kwargs):
        with STORAGE_LOCK:
            return function(*args, **kwargs)
    return decorated


class JobRunner:
    """One background worker: pair scans are CPU and memory heavy, so jobs run one at a time in order.

    With inline=True (tests) a submitted job runs before submit() returns.
    """

    def __init__(self, logger, inline=False):
        self.logger = logger
        self.inline = inline
        self.executor = None if inline else ThreadPoolExecutor(max_workers=1, thread_name_prefix='split-job')
        self.pending = []
        self.lock = threading.Lock()

    def submit(self, service, identifier):
        with self.lock:
            self.pending.append(identifier)
        if self.inline:
            self._work(service, identifier)
        else:
            self.executor.submit(self._work, service, identifier)

    def position(self, identifier):
        """0 while running, n when n jobs are ahead in the queue, None when not tracked."""
        with self.lock:
            return self.pending.index(identifier) if identifier in self.pending else None

    def _work(self, service, identifier):
        try:
            service.execute(identifier)
        except Exception:
            self.logger.exception('Dual-probe job %s failed', identifier)
        finally:
            with self.lock:
                if identifier in self.pending:
                    self.pending.remove(identifier)


def register_split_routes(app, bridges, cdna, cds, aliases, available_genomes):
    def workflow():
        return SplitWorkflow(app.config['UPLOAD_FOLDER'], bridges, cdna, cds, aliases,
                             reference_library=app.extensions.get('reference_library'))

    runner = JobRunner(app.logger, inline=app.config.get('SPLIT_JOBS_INLINE', False))
    app.extensions['split_jobs'] = runner
    recovered = workflow().recover()
    if recovered:
        app.logger.warning('Marked %d unfinished dual-probe jobs as interrupted', len(recovered))

    def api_errors(function):
        @wraps(function)
        def decorated(*args, **kwargs):
            try:
                return function(*args, **kwargs)
            except InputError as error:
                return jsonify(error=str(error), field=error.field), error.status
            except Exception:
                app.logger.exception('Dual-probe web workflow failed')
                return jsonify(error='设计未完成。请检查参考文件与比对环境后重试；已填写内容仍保留。'), 500
        return decorated

    def links(identifier):
        return dict(status_url=url_for('split_status', identifier=identifier),
                    result_url=url_for('split_result', identifier=identifier),
                    alignment_url=url_for('split_alignment', identifier=identifier),
                    download_url=url_for('split_download', identifier=identifier),
                    run_url=url_for('split_run', identifier=identifier))

    def status_payload(service, identifier):
        status = service.status(identifier)
        status.update(links(identifier), queue_position=runner.position(identifier))
        return status

    def page(mode='single'):
        genomes = available_genomes()
        default_reference = next((g['filename'] for g in genomes if g['filename'] == 'bol_bju_cdna.fasta'), '')
        config = dict(mode=mode, defaults=DEFAULTS, bridges=bridges,
                      single_url=url_for('split'), batch_url=url_for('split_batch'),
                      prepare_url=url_for('split_prepare'), jobs_url=url_for('jobs'),
                      default_reference=default_reference)
        return render_template('split.html', config=config, available_genomes=genomes)

    def legacy(mode):
        service = workflow()
        with STORAGE_LOCK:
            prepared = service.prepare(request.form, request.files, mode=mode)
            service.queue(prepared['id'], request.form.to_dict())
        service.execute(prepared['id'])
        return send_file(service.job_dir(prepared['id']) / 'results.zip', as_attachment=True,
                         download_name=f"split_{mode}_results_{prepared['id']}.zip")

    @app.route('/split', methods=['GET', 'POST'])
    @api_errors
    def split():
        return legacy('single') if request.method == 'POST' else page()

    @app.route('/split/batch', methods=['GET', 'POST'])
    @api_errors
    def split_batch():
        return legacy('batch') if request.method == 'POST' else page('batch')

    @app.post('/split/prepare')
    @api_errors
    def split_prepare():
        with STORAGE_LOCK:
            prepared = workflow().prepare(request.form, request.files)
        prepared['run_url'] = url_for('split_run', identifier=prepared['id'])
        return jsonify(prepared)

    @app.post('/split/jobs/<identifier>/run')
    @api_errors
    def split_run(identifier):
        """Queue the job and answer at once; the browser follows status_url."""
        service = workflow()
        with STORAGE_LOCK:
            queued = service.queue(identifier, request.get_json(silent=True))
        if queued:
            runner.submit(service, identifier)
            return jsonify(status_payload(service, identifier)), 202
        return jsonify(status_payload(service, identifier)), 200

    @app.get('/split/jobs/<identifier>/status')
    @api_errors
    @storage_locked
    def split_status(identifier):
        return jsonify(status_payload(workflow(), identifier))

    @app.get('/split/results/<identifier>')
    @api_errors
    @storage_locked
    def split_result(identifier):
        service = workflow()
        status = status_payload(service, identifier)
        if request.accept_mimetypes.best == 'application/json':
            return jsonify(dict(service.result(identifier), **links(identifier)))
        result = dict(service.result(identifier), **links(identifier)) if status['status'] == 'completed' else None
        if status['status'] == 'prepared':
            raise InputError('该任务尚未生成完整结果。', status=409)
        return render_template('split_result.html', status=status, result=result, jobs_url=url_for('jobs'))

    @app.get('/split/results/<identifier>/alignment')
    @api_errors
    def split_alignment(identifier):
        # The slow, read-only BLAST request must not hold the shared job/storage lock.
        index = request.args.get('task', '0')
        if not index.isdecimal() or len(index) > 6:
            raise InputError('请选择有效的任务编号。')
        return jsonify(workflow().homolog_alignment(identifier, int(index)))

    @app.get('/split/results/<identifier>/download')
    @api_errors
    @storage_locked
    def split_download(identifier):
        service = workflow()
        service.result(identifier)
        return send_file(service.job_dir(identifier) / 'results.zip', as_attachment=True,
                         download_name=f'split_results_{identifier}.zip')

    @app.get('/split/batch-template.csv')
    def split_batch_template():
        example = 'task_name,BP_ID,gene_id,max_selected,allowed_ids\nABCG13,BP_0079,BjuVA06G01970,5,\n'
        return send_file(io.BytesIO(example.encode('utf-8-sig')), mimetype='text/csv',
                         as_attachment=True, download_name='split_batch_template.csv')

    @app.get('/jobs')
    @api_errors
    @storage_locked
    def jobs():
        service = workflow()
        listing = [dict(job, **links(job['id']), queue_position=runner.position(job['id']))
                   for job in service.list_jobs()]
        if request.accept_mimetypes.best == 'application/json':
            return jsonify(jobs=listing)
        return render_template('jobs.html', jobs=listing)
