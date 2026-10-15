"""Request-local zh/en UI translations, shared by Jinja and browser scripts.

Chinese messages are the source keys. Only presentation strings are translated;
stored jobs, user input, sequence data and download formats remain unchanged.
"""
import json
import re
from pathlib import Path
from urllib.parse import urlencode

from flask import g, has_request_context, request


LOCALES = Path(__file__).resolve().parents[1] / 'locales'
LANGUAGES = {'zh': 'zh-CN', 'en': 'en'}
COOKIE = 'hcr_language'
PLACEHOLDER = re.compile(r'\{(\w+)\}')


def _catalog(filename):
    with (LOCALES / filename).open(encoding='utf-8') as handle:
        return json.load(handle)


ENGLISH = _catalog('en.json')
MESSAGES = _catalog('messages.en.json')


def _message_patterns():
    patterns = []
    for source, target in MESSAGES.items():
        keys = PLACEHOLDER.findall(source)
        if not keys:
            continue
        parts = PLACEHOLDER.split(source)
        pattern = ''.join(re.escape(part) if i % 2 == 0 else '(.+?)'
                          for i, part in enumerate(parts))
        specificity = sum(len(part) for part in parts[::2])
        translated_values = {'error'}
        if source.startswith('{label}需为 ') or source == '{label}不能大于{error}。':
            translated_values.add('label')
        patterns.append((specificity, re.compile(pattern, re.DOTALL), keys, target, translated_values))
    return sorted(patterns, key=lambda item: item[0], reverse=True)


MESSAGE_PATTERNS = _message_patterns()


def language():
    if not has_request_context():
        return 'zh'
    if 'hcr_language' not in g:
        explicit = request.args.get('lang')
        saved = request.cookies.get(COOKIE)
        g.hcr_language = explicit if explicit in LANGUAGES else saved if saved in LANGUAGES else 'zh'
    return g.hcr_language


def translate(source, **values):
    if source is None:
        return ''
    text = str(source)
    if language() == 'en':
        text = ENGLISH.get(text, MESSAGES.get(text, text))
    return PLACEHOLDER.sub(lambda match: str(values.get(match[1], match[0])), text)


def message(source, depth=0):
    """Localize existing validation messages at the HTTP/display boundary.

    Legacy services format errors before returning them. Match their catalogued
    templates here, preserving variable values such as filenames and transcript
    IDs. Nested error wrappers and parameter labels can themselves be translated.
    This is deliberately limited to error/status fields, never arbitrary data.
    """
    if source is None:
        return ''
    text = str(source)
    if language() != 'en' or depth > 5:
        return text
    if text in MESSAGES or text in ENGLISH:
        return translate(text)
    for delimiter, localized in (('\n', '\n'), ('；', '; ')):
        if delimiter in text:
            return localized.join(message(part, depth + 1) for part in text.split(delimiter))
    for _, pattern, keys, target, translated_values in MESSAGE_PATTERNS:
        matched = pattern.fullmatch(text)
        if matched:
            values = dict(zip(keys, matched.groups()))
            for key in translated_values:
                if key in values:
                    values[key] = message(values[key], depth + 1)
            return PLACEHOLDER.sub(lambda match: values[match[1]], target)
    return text


def language_url(selected):
    query = [(key, value) for key, values in request.args.lists() if key != 'lang' for value in values]
    return request.script_root + request.path + '?' + urlencode(query + [('lang', selected)])


def localize_status(value):
    """Return display copies; language choices must never modify persisted jobs."""
    if isinstance(value, list):
        return [localize_status(item) for item in value]
    if not isinstance(value, dict):
        return value
    item = dict(value)
    if isinstance(item.get('error'), str):
        item['error'] = message(item['error'])
    if isinstance(item.get('progress'), dict) and isinstance(item['progress'].get('label'), str):
        item['progress'] = dict(item['progress'], label=translate(item['progress']['label']))
    return item


def register_i18n(app):
    """Idempotent so individual route modules also work in isolated Flask apps."""
    if 'i18n' in app.extensions:
        return
    app.extensions['i18n'] = True
    app.jinja_env.globals.update(_=translate, language_url=language_url)
    app.jinja_env.filters['message'] = message
    app.jinja_env.filters['localized_status'] = localize_status

    @app.context_processor
    def i18n_context():
        selected = language()
        return dict(language=selected, html_language=LANGUAGES[selected],
                    i18n_catalog=dict(MESSAGES, **ENGLISH) if selected == 'en' else {})

    @app.after_request
    def localize_response(response):
        if request.endpoint == 'static':
            return response
        selected = language()
        response.headers['Content-Language'] = LANGUAGES[selected]
        response.vary.add('Cookie')
        if request.args.get('lang') in LANGUAGES:
            response.set_cookie(COOKIE, selected, max_age=365 * 24 * 60 * 60,
                                httponly=True, samesite='Lax',
                                secure=request.is_secure or app.config.get('SESSION_COOKIE_SECURE', False),
                                path=app.config.get('APPLICATION_ROOT') or '/')
        if selected == 'en' and response.is_json and not response.direct_passthrough:
            payload = response.get_json(silent=True)

            payload = localize_status(payload)
            if isinstance(payload, dict) and isinstance(payload.get('jobs'), list):
                payload['jobs'] = localize_status(payload['jobs'])
            if (request.endpoint == 'split_alignment' and isinstance(payload, dict)
                    and isinstance(payload.get('message'), str)):
                payload['message'] = message(payload['message'])
            if payload is not None:
                response.set_data(app.json.dumps(payload))
        return response
