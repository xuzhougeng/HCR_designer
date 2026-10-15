"""Public reference uploads and configuration; deletion remains an administrator action."""
from functools import wraps
from pathlib import Path

from flask import jsonify, render_template, request

from scripts.reference_library import ReferenceLibrary


def register_reference_routes(app):
    from scripts.i18n import register_i18n
    register_i18n(app)
    app.config.setdefault('REFERENCE_CATALOG_FILE', str(Path(app.instance_path) / 'references.json'))
    app.config.setdefault('REFERENCE_PUBLIC_EDIT', True)
    library = ReferenceLibrary(app.config['UPLOAD_FOLDER'], app.config['REFERENCE_CATALOG_FILE'])
    app.extensions['reference_library'] = library
    auth = app.extensions['admin_auth']

    def can_edit():
        return bool(app.config['REFERENCE_PUBLIC_EDIT'] or auth['authenticated']())

    app.jinja_env.globals['reference_can_edit'] = can_edit
    app.extensions['reference_can_edit'] = can_edit

    def writable(function):
        @wraps(function)
        def decorated(*args, **kwargs):
            if not can_edit():
                return jsonify(error='请先登录管理员账号再管理参考库。'), 401
            auth['csrf']()
            try:
                return function(*args, **kwargs)
            except ValueError as error:
                return jsonify(error=str(error)), 400
        return decorated

    @app.get('/references')
    def references():
        listing = library.list()
        if request.accept_mimetypes.best == 'application/json':
            return jsonify(references=listing, can_edit=can_edit())
        return render_template('references.html', references=listing, can_edit=can_edit(),
                               can_delete=bool(auth['authenticated']()))

    @app.post('/references/upload')
    @writable
    def reference_upload():
        upload = request.files.get('reference')
        if not upload or not upload.filename:
            return jsonify(error='请选择参考 FASTA 文件。'), 400
        ploidy = request.form.get('polyploid', '')
        if ploidy not in ('', 'true', 'false'):
            raise ValueError('请选择有效的多倍体配置。')
        result = library.upload(upload, request.form.get('display_name', ''),
                                None if not ploidy else ploidy == 'true', request.files.get('subgenome_map'))
        return jsonify(result), 200 if result['duplicate'] else 201

    @app.post('/references/configure')
    @writable
    def reference_configure():
        data = request.get_json(silent=True)
        if request.mimetype == 'multipart/form-data':
            ploidy = request.form.get('polyploid', '')
            if ploidy not in ('', 'true', 'false'):
                raise ValueError('请选择有效的多倍体配置。')
            data = dict(request.form, polyploid=None if not ploidy else ploidy == 'true')
        if not isinstance(data, dict):
            raise ValueError('请提供参考库配置。')
        return jsonify(reference=library.configure(data.get('filename'), data.get('display_name'), data.get('polyploid'),
                                                   request.files.get('subgenome_map')))

    @app.post('/references/hash')
    @writable
    def reference_hash():
        data = request.get_json(silent=True)
        return jsonify(reference=library.index(data.get('filename') if isinstance(data, dict) else None))
