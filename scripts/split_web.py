"""HTTP adapters for the shared single/batch guided interface."""
from functools import wraps
import io

from flask import jsonify, render_template, request, send_file, url_for
from scripts.split_workflow import SplitWorkflow, InputError, DEFAULTS


def register_split_routes(app, bridges, cdna, cds, aliases, available_genomes):
    def workflow():
        return SplitWorkflow(app.config['UPLOAD_FOLDER'], bridges, cdna, cds, aliases)

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

    def result_links(result):
        return dict(result, download_url=url_for('split_download', identifier=result['id']),
                    result_url=url_for('split_result', identifier=result['id']))

    def page(mode='single', result=None):
        genomes = available_genomes()
        default_reference = next((g['filename'] for g in genomes if g['filename'] == 'bol_bju_cdna.fasta'), '')
        config = dict(mode=mode, defaults=DEFAULTS, bridges=bridges,
                      single_url=url_for('split'), batch_url=url_for('split_batch'),
                      prepare_url=url_for('split_prepare'),
                      default_reference=default_reference, result=result)
        return render_template('split.html', config=config, available_genomes=genomes)

    def legacy(mode):
        service = workflow()
        prepared = service.prepare(request.form, request.files, mode=mode)
        service.run(prepared['id'], request.form.to_dict())
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
        prepared = workflow().prepare(request.form, request.files)
        prepared['run_url'] = url_for('split_run', identifier=prepared['id'])
        return jsonify(prepared)

    @app.post('/split/jobs/<identifier>/run')
    @api_errors
    def split_run(identifier):
        return jsonify(result_links(workflow().run(identifier, request.get_json(silent=True))))

    @app.get('/split/results/<identifier>')
    @api_errors
    def split_result(identifier):
        result = result_links(workflow().result(identifier))
        return page(result['mode'], result)

    @app.get('/split/results/<identifier>/download')
    @api_errors
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
