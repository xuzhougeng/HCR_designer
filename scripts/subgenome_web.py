"""A stateless sequence tool: uploaded FASTA and generated CSV are never kept in uploads/."""
import json

from flask import jsonify, render_template, request

from src.common.subgenome_generator import generate_subgenomes


def register_subgenome_routes(app):
    from scripts.i18n import register_i18n
    register_i18n(app)
    @app.route('/subgenomes', methods=['GET', 'POST'])
    def subgenomes():
        if request.method == 'GET':
            return render_template('subgenomes.html')
        try:
            mode = request.form.get('mode', 'files')
            uploads = [upload for upload in request.files.getlist('fasta_files') if upload.filename]
            labels = json.loads(request.form['labels']) if 'labels' in request.form else None
            rules = json.loads(request.form['rules']) if 'rules' in request.form else None
            return jsonify(generate_subgenomes(
                uploads, mode, labels=labels, rules=rules,
                position=request.form.get('position', 'contains'),
                case_sensitive=request.form.get('case_sensitive') == 'on',
                rule_mode=request.form.get('rule_mode', 'feature')))
        except (ValueError, TypeError) as error:
            return jsonify(error=str(error)), 400
        except Exception:
            app.logger.exception('Subgenome annotation generation failed')
            return jsonify(error='生成未完成，请检查 FASTA 文件后重试。'), 500
