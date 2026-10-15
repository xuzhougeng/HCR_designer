"""Administrator authentication, storage UI and explicit maintenance commands."""
from datetime import timedelta
from functools import wraps
import hmac
import json
import os
from pathlib import Path
import secrets
import threading
import time

import click
from flask import abort, g, jsonify, redirect, render_template, request, session, url_for
from werkzeug.security import check_password_hash, generate_password_hash

from scripts.storage_maintenance import CleanupScheduler, STORAGE_LOCK, StorageMaintenance, three_months_ago, utcnow


def register_admin(app):
    from scripts.i18n import register_i18n
    register_i18n(app)
    app.config.setdefault('ADMIN_CREDENTIALS_FILE', str(Path(app.instance_path) / 'admin.json'))
    app.config.update(SESSION_COOKIE_HTTPONLY=True, SESSION_COOKIE_SAMESITE='Lax',
                      PERMANENT_SESSION_LIFETIME=timedelta(hours=8))
    app.config['SESSION_COOKIE_SECURE'] = os.environ.get('HCR_COOKIE_SECURE', '0') == '1'

    def credentials():
        path = Path(app.config['ADMIN_CREDENTIALS_FILE'])
        if not path.exists():
            return None
        value = json.loads(path.read_text(encoding='utf-8'))
        if not all(isinstance(value.get(key), str) and value[key] for key in
                   ('username', 'password_hash', 'secret_key', 'version')):
            raise ValueError('Invalid administrator configuration')
        return value

    configured = credentials()
    app.secret_key = (configured or {}).get('secret_key') or secrets.token_hex(32)
    storage = StorageMaintenance(app.config['UPLOAD_FOLDER'])
    app.extensions['storage_maintenance'] = storage
    attempts, attempts_lock = {}, threading.Lock()
    # PBKDF2 also works on the project's macOS Python builds without hashlib.scrypt.
    password_method = 'pbkdf2:sha256:1000000'
    dummy_hash = generate_password_hash(secrets.token_hex(32), method=password_method)

    def csrf_token():
        if 'admin_csrf' not in session:
            session['admin_csrf'] = secrets.token_urlsafe(32)
        return session['admin_csrf']

    app.jinja_env.globals['admin_csrf_token'] = csrf_token

    def check_csrf():
        expected = session.get('admin_csrf', '')
        supplied = request.headers.get('X-CSRF-Token') or request.form.get('csrf_token', '')
        if not expected or not isinstance(supplied, str) or not hmac.compare_digest(expected.encode(), supplied.encode()):
            abort(400, description='页面验证已失效，请刷新后重试。')

    def authenticated():
        current = credentials()
        return current and session.get('admin_version') == current['version']

    def admin_only(function):
        @wraps(function)
        def decorated(*args, **kwargs):
            if not authenticated():
                if request.method != 'GET' or request.accept_mimetypes.best == 'application/json':
                    return jsonify(error='请先登录管理员账号。'), 401
                return redirect(url_for('admin_login'))
            if request.method == 'POST':
                check_csrf()
            return function(*args, **kwargs)
        return decorated

    app.extensions['admin_auth'] = dict(authenticated=authenticated, required=admin_only, csrf=check_csrf)

    @app.after_request
    def private_admin_responses(response):
        if request.path == '/admin' or request.path.startswith('/admin/'):
            response.headers['Cache-Control'] = 'no-store'
            response.headers['X-Frame-Options'] = 'DENY'
        return response

    # Legacy triplet requests have no persisted status; track their whole lifetime.
    @app.before_request
    def protect_legacy_jobs():
        if request.method == 'POST' and request.endpoint in ('triplet', 'triplet_batch'):
            with STORAGE_LOCK:
                storage.legacy_active += 1
                g.storage_legacy_active = True

    @app.teardown_request
    def release_legacy_jobs(error):
        if g.pop('storage_legacy_active', False):
            with STORAGE_LOCK:
                storage.legacy_active -= 1

    @app.route('/admin/login', methods=['GET', 'POST'])
    def admin_login():
        destination = request.args.get('next', 'admin')
        if destination not in ('admin', 'references') or destination not in app.view_functions:
            destination = 'admin'
        current = credentials()
        if current is None:
            return render_template('admin_login.html', initialized=False), 503
        if authenticated():
            return redirect(url_for(destination))
        error, status = None, 200
        if request.method == 'POST':
            check_csrf()
            address, now = request.remote_addr or 'unknown', time.monotonic()
            with attempts_lock:
                for key in list(attempts):
                    attempts[key] = [stamp for stamp in attempts[key] if now - stamp < 300]
                    if not attempts[key]:
                        del attempts[key]
                limited = len(attempts.get(address, [])) >= 5
                if not limited:
                    attempts.setdefault(address, []).append(now)
            if limited:
                error, status = '登录尝试过多，请在 5 分钟后重试。', 429
            else:
                username = request.form.get('username', '')
                matches = hmac.compare_digest(username.encode(), current['username'].encode())
                valid = check_password_hash(current['password_hash'] if matches else dummy_hash,
                                            request.form.get('password', ''))
                if matches and valid:
                    with attempts_lock:
                        attempts.pop(address, None)
                    session.clear()
                    session.permanent = True
                    session['admin_version'] = current['version']
                    return redirect(url_for(destination))
                error, status = '账号或密码错误。', 401
        return render_template('admin_login.html', initialized=True, error=error), status

    @app.post('/admin/logout')
    @admin_only
    def admin_logout():
        session.clear()
        return redirect(url_for('admin_login'))

    @app.get('/admin')
    @admin_only
    def admin():
        entries = [{k: v for k, v in entry.items() if not k.startswith('_')} for entry in storage.inventory()]
        data = dict(entries=entries, total_bytes=sum(e['size'] for e in entries),
                    expired_bytes=sum(e['size'] for e in entries if e['expired'] and not e['blocked']),
                    cutoff=three_months_ago(utcnow()).isoformat(), last_cleanup=storage.last_cleanup,
                    automatic_cleanup='storage_cleanup_scheduler' in app.extensions)
        if request.accept_mimetypes.best == 'application/json':
            return jsonify(data)
        return render_template('admin.html', **data)

    @app.post('/admin/delete')
    @admin_only
    def admin_delete():
        body = request.get_json(silent=True)
        try:
            report = storage.delete(keys=body.get('keys') if isinstance(body, dict) else None)
        except ValueError as error:
            return jsonify(error=str(error)), 409
        app.logger.info('Administrator deleted storage entries: %s', report)
        return jsonify(report), 207 if report['errors'] else 200

    @app.post('/admin/cleanup')
    @admin_only
    def admin_cleanup():
        report = storage.cleanup()
        app.logger.info('Administrator ran result retention cleanup: %s', report)
        return jsonify(report), 207 if report['errors'] else 200

    @app.cli.command('init-admin')
    @click.option('--username', prompt='管理员账号')
    @click.option('--reset', is_flag=True, help='重设已有管理员账号及密码，使已有登录失效。')
    def init_admin(username, reset):
        """Initialize credentials locally; passwords are never accepted on the command line."""
        path = Path(app.config['ADMIN_CREDENTIALS_FILE'])
        current = credentials()
        if current and not reset:
            raise click.ClickException('管理员已初始化；需要重设时请使用 --reset。')
        username = username.strip()
        if not username or len(username) > 80:
            raise click.ClickException('账号需为 1–80 个字符。')
        password = click.prompt('管理员密码（至少 12 位）', hide_input=True, confirmation_prompt=True)
        if len(password) < 12:
            raise click.ClickException('密码至少需要 12 位。')
        value = dict(username=username, password_hash=generate_password_hash(password, method=password_method),
                     secret_key=(current or {}).get('secret_key') or secrets.token_hex(32),
                     version=secrets.token_hex(32))
        path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        temporary = path.with_name(path.name + '.' + secrets.token_hex(8) + '.tmp')
        try:
            fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            with os.fdopen(fd, 'w', encoding='utf-8') as handle:
                json.dump(value, handle)
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
        app.secret_key = value['secret_key']
        click.echo('管理员账号已保存，密码仅保存哈希。请重启 Web 服务后访问 /admin。')

    @app.cli.command('cleanup-results')
    @click.option('--dry-run', is_flag=True, help='只预览超过三个月的结果，不删除文件。')
    def cleanup_results(dry_run):
        """Remove expired analysis results; shared references are retained."""
        report = storage.cleanup(dry_run=dry_run)
        click.echo(json.dumps(report, ensure_ascii=False, indent=2))
        if report['errors']:
            raise click.ClickException('部分结果未能清理，请检查文件权限。')

    # Start with the first server request, never as a side effect of an import or CLI command.
    if os.environ.get('HCR_CLEANUP_ENABLED') == '1' and not app.testing:
        scheduler = CleanupScheduler(storage, app.logger)
        app.extensions['storage_cleanup_scheduler'] = scheduler
        scheduler_lock = threading.Lock()

        @app.before_request
        def start_storage_cleanup():
            with scheduler_lock:
                if not app.testing and scheduler.thread.ident is None:
                    scheduler.start()
