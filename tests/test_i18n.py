import ast
import io
import json
from html.parser import HTMLParser
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from flask import Flask, jsonify, render_template, render_template_string
from jinja2 import Environment

from scripts.i18n import ENGLISH, MESSAGES, PLACEHOLDER, localize_status, message, register_i18n


ROOT = Path(__file__).resolve().parents[1]
HAN = re.compile(r'[\u3400-\u9fff]')
PAGES = ['/', '/split', '/split/batch', '/triplet', '/triplet/batch', '/bridge', '/bridge/list',
         '/bridge/query', '/bridge/sequence', '/bridge/analyze', '/bridge/validate', '/formatting',
         '/probe_recommend', '/references', '/subgenomes', '/help', '/jobs', '/orders']


class VisibleText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hidden = False
        self.text = []

    def handle_starttag(self, tag, attrs):
        if tag in ('script', 'style'):
            self.hidden = True
        for name, value in attrs:
            if name in ('aria-label', 'title', 'placeholder') and value:
                self.text.append(value)

    def handle_endtag(self, tag):
        if tag in ('script', 'style'):
            self.hidden = False

    def handle_data(self, value):
        if not self.hidden:
            self.text.append(value)


class LocaleTests(unittest.TestCase):
    def setUp(self):
        self.app = Flask(__name__)
        register_i18n(self.app)
        register_i18n(self.app)  # Route modules may register independently.
        self.original = {'error': '第 2 行：在所选序列库中找不到 帮助，请核对完整转录本编号。',
                         'progress': {'stage': 'candidates', 'label': '生成候选'},
                         'tasks': [{'name': '帮助', 'gene_id': '基因1', 'sequence': 'ATCG'}]}

        @self.app.get('/')
        def page():
            return render_template_string('<html lang="{{ html_language }}">{{ _("双探针设计") }}'
                                          '<a href="{{ language_url(\'en\') }}">English</a></html>')

        @self.app.get('/status')
        def status():
            return jsonify(self.original)

        @self.app.get('/jobs')
        def jobs():
            return jsonify(jobs=[self.original])

        @self.app.get('/alignment', endpoint='split_alignment')
        def alignment():
            return jsonify(status='no_match', subject_id='帮助',
                           message='归属文件中没有可用于比对的其他分组序列。')

        self.client = self.app.test_client()

    def test_default_and_cookie_persistence_in_both_directions(self):
        self.assertIn('<html lang="zh-CN">双探针设计', self.client.get('/').text)
        english = self.client.get('/?lang=en')
        self.assertIn('<html lang="en">Split-probe design', english.text)
        cookie = english.headers['Set-Cookie']
        self.assertIn('hcr_language=en', cookie)
        self.assertIn('HttpOnly', cookie)
        self.assertIn('SameSite=Lax', cookie)
        self.assertIn('Split-probe design', self.client.get('/').text)
        self.client.get('/?lang=zh')
        self.assertIn('双探针设计', self.client.get('/').text)

    def test_invalid_language_uses_saved_or_default_language(self):
        self.assertEqual(self.client.get('/?lang=fr').headers['Content-Language'], 'zh-CN')
        self.client.get('/?lang=en')
        response = self.client.get('/?lang=../../invalid')
        self.assertEqual(response.headers['Content-Language'], 'en')
        self.assertNotIn('Set-Cookie', response.headers)

    def test_language_is_request_local_and_marked_for_caches(self):
        self.client.get('/?lang=en')
        self.assertIn('双探针设计', self.app.test_client().get('/').text)
        response = self.client.get('/')
        self.assertEqual(response.headers['Content-Language'], 'en')
        self.assertIn('Cookie', response.vary)
        secure = self.app.test_client().get('/?lang=en', base_url='https://example.test')
        self.assertIn('Secure', secure.headers['Set-Cookie'])

    def test_language_links_keep_repeated_query_parameters(self):
        page = self.client.get('/?q=Bo&page=2&tag=A&tag=B&lang=zh').text
        self.assertIn('/?q=Bo&amp;page=2&amp;tag=A&amp;tag=B&amp;lang=en', page)

    def test_json_errors_and_progress_are_localized_without_changing_user_data(self):
        for path in ('/status?lang=en', '/jobs?lang=en'):
            response = self.client.get(path).json
            status = response['jobs'][0] if 'jobs' in response else response
            self.assertEqual(status['error'], 'Row 2: Could not find 帮助 in the selected sequence library. Check the full transcript ID.')
            self.assertEqual(status['progress']['label'], 'Generate candidates')
            self.assertEqual(status['progress']['stage'], 'candidates')
            self.assertEqual(status['tasks'], self.original['tasks'])
        self.assertEqual(self.original['progress']['label'], '生成候选')
        self.assertTrue(self.original['error'].startswith('第 2 行'))

    def test_formatted_and_multiple_validation_messages(self):
        with self.app.test_request_context('/?lang=en'):
            self.assertEqual(message('最小臂长需为 15–20 之间的整数。'),
                             'Minimum arm length must be an integer between 15 and 20.')
            self.assertEqual(message('最低 Tm不能大于最高 Tm。'), 'Minimum Tm cannot exceed Maximum Tm.')
            self.assertEqual(message('长度不符合要求: 3bp；GC含量不适合: 0.0%'),
                             'Length outside the allowed range: 3 bp; GC content outside the allowed range: 0.0%')
            self.assertEqual(message('第 1 行：CSV 第 1 行需为 3–5 列。\n第 2 行：目标序列至少需要 30 nt。'),
                             'Row 1: CSV row 1 must have 3–5 columns.\nRow 2: Target sequence must be at least 30 nt.')
            self.assertEqual(message('无法读取 中文参考.fa，请检查文件是否完整且为 UTF-8 FASTA。'),
                             'Could not read 中文参考.fa. Check file integrity and UTF-8 FASTA encoding.')
            self.assertEqual(message('亚基因组名称 帮助 重复，请为每个文件使用不同名称。'),
                             'Duplicate subgenome name 帮助. Use a different name for each file.')

    def test_status_display_copy_and_html_escaping(self):
        with self.app.test_request_context('/?lang=en'):
            display = localize_status(self.original)
            self.assertEqual(display['progress']['label'], 'Generate candidates')
            rendered = render_template_string('{{ error|message }}', error='桥接编号 <script>alert(1)</script> 不在当前探针库中。')
            self.assertIn('Bridge ID &lt;script&gt;', rendered)
            self.assertNotIn('<script>', rendered)
        self.assertEqual(self.original['progress']['label'], '生成候选')

    def test_homolog_alignment_messages_follow_the_selected_language(self):
        result = self.client.get('/alignment?lang=en').json
        self.assertEqual(result['status'], 'no_match')
        self.assertEqual(result['subject_id'], '帮助')
        self.assertEqual(result['message'], 'No sequences from other annotated groups are available for alignment.')
        self.assertEqual(self.client.get('/alignment?lang=zh').json['message'],
                         '归属文件中没有可用于比对的其他分组序列。')


class TranslatedPageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import app
        cls.web = app

    def setUp(self):
        self.client = self.web.app.test_client()

    def test_all_public_pages_render_english_including_accessibility_text(self):
        for path in PAGES:
            with self.subTest(path=path):
                response = self.client.get(path + '?lang=en')
                self.assertEqual(response.status_code, 200)
                self.assertIn('<html lang="en">', response.text)
                visible = VisibleText()
                visible.feed(response.text)
                self.assertEqual([s for s in visible.text if HAN.search(s) and s.strip() != '中文'], [])
                self.assertIn('/static/js/i18n.js', response.text)

    def test_existing_chinese_pages_and_user_reference_names_are_preserved(self):
        reference = dict(filename='reference.fa', display_name='我的参考库')
        with patch.object(self.web, 'get_available_genomes', return_value=[reference]):
            page = self.client.get('/?lang=en').text
        self.assertIn('我的参考库', page)
        chinese = self.client.get('/?lang=zh').text
        self.assertIn('从转录本到可复核的探针清单', chinese)

    def test_real_validation_endpoints_and_downloads(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(self.web.app.config, UPLOAD_FOLDER=folder):
            response = self.client.post('/split/prepare?lang=en', data={
                'input_mode': 'single', 'name': '测试任务', 'geneID': 'BjuTarget', 'bp_id': 'BP_9999',
                'sequenceType': 'custom', 'seq': 'ACGT' * 20, 'reference_mode': 'none'})
            self.assertEqual(response.status_code, 400)
            self.assertIn('Bridge ID BP_9999', response.json['error'])
            response = self.client.post('/triplet?lang=en', data={
                'name': 'test', 'bp_id': 'BP_9999', 'geneID': '', 'seq': 'ACGT'})
        self.assertIn('Bridge probe ID BP_9999', response.json['error'])
        for lang in ('zh', 'en'):
            file = self.client.get('/split/batch-template.csv?lang=' + lang)
            self.assertIn('task_name,BP_ID,gene_id,max_selected,allowed_ids', file.text)
            self.assertIn('attachment;', file.headers['Content-Disposition'])

    def test_legacy_validation_reasons_are_english(self):
        response = self.client.post('/bridge/validate?lang=en', data={
            'probe_file': (io.BytesIO(b'id,sequence\np1,AAA\n'), 'probes.csv')},
            content_type='multipart/form-data')
        self.assertEqual(response.status_code, 200)
        self.assertIn('Length outside the allowed range', response.text)

    def test_administrator_categories_and_protection_reasons_are_english(self):
        date = '2026-10-15T12:00:00+00:00'
        entry = dict(key='reference.fa', name='实验参考.fa', kind='参考基因组', size=1024,
                     modified_at=date, file_count=2, blocked='参考库正在被任务使用', expired=False)
        with self.web.app.test_request_context('/admin?lang=en'):
            page = render_template('admin.html', entries=[entry], total_bytes=1024, expired_bytes=0,
                                   automatic_cleanup=True, cutoff=date,
                                   last_cleanup=dict(ran_at=date, deleted=[], skipped=[], errors=[]))
        self.assertIn('Reference genome', page)
        self.assertIn('Reference is in use by a job', page)
        self.assertIn('实验参考.fa', page)
        self.assertIn('value="参考基因组"', page)  # Filtering keys stay stable across languages.

    def test_polyploid_results_keep_alignment_and_evidence_navigation_in_english(self):
        status = dict(tasks=[dict(name='测试任务')], id='a' * 32, status='completed',
                      mode='single', reference='我的参考库', download_url='/download')
        with self.web.app.test_request_context('/split/results/' + status['id'] + '?lang=en'):
            page = render_template('split_result.html', status=status, result={'tasks': status['tasks']}, jobs_url='/jobs')
        for text in ('Polyploid mode', 'Most similar sequence in another group', 'Previous pair', 'Next pair',
                     'id="polyploid-view"', 'id="homolog-panel"', '/static/js/result-evidence.js', '测试任务', '我的参考库'):
            self.assertIn(text, page)


class CatalogTests(unittest.TestCase):
    def test_all_catalog_entries_are_translated_with_matching_placeholders(self):
        for source, target in dict(MESSAGES, **ENGLISH).items():
            with self.subTest(source=source):
                self.assertTrue(target)
                self.assertIsNone(HAN.search(target))
                self.assertEqual(set(PLACEHOLDER.findall(source)), set(PLACEHOLDER.findall(target)))

    def test_all_explicit_template_and_javascript_messages_have_translations(self):
        catalog = dict(MESSAGES, **ENGLISH)
        for path in ROOT.glob('templates/**/*.html'):
            if path.name in ('hcr.html', 'snail.html'):  # Historical, unrouted templates.
                continue
            for _, kind, raw in Environment().lex(path.read_text()):
                if kind == 'string':
                    value = ast.literal_eval(raw)
                    if HAN.search(value):
                        self.assertIn(value, catalog, str(path))
        calls = re.compile(r'''HCR\.t\(("(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')''')
        for path in [*ROOT.glob('static/**/*.js'), *ROOT.glob('templates/**/*.html')]:
            for match in calls.finditer(path.read_text()):
                value = ast.literal_eval(match[1])
                self.assertIn(value, catalog, str(path))

    @unittest.skipUnless(shutil.which('node'), 'Node is required for browser helper verification')
    def test_browser_translation_and_saved_order_labels(self):
        # Run the real browser helpers without DOM rendering or external packages.
        script = '''
const fs = require('fs'), vm = require('vm'), assert = require('assert');
const messages = JSON.parse(fs.readFileSync('locales/en.json', 'utf8'));
const config = {language: 'en', messages};
global.window = {addEventListener() {}};
global.document = {getElementById() {return {textContent: JSON.stringify(config)}}, addEventListener() {}};
vm.runInThisContext(fs.readFileSync('static/js/i18n.js', 'utf8'));
global.HCR = window.HCR;
vm.runInThisContext(fs.readFileSync('static/js/app.js', 'utf8'));
assert.equal(HCR.locale, 'en-US');
assert.equal(HCR.t('请求未完成（HTTP {status}）。', {status: 500}), 'Request incomplete (HTTP 500).');
assert.equal(HCR.t('选择 {name}', {name: '$&<中文>'}), 'Select $&<中文>');
assert.equal(HCR.t('unknown message'), 'unknown message');
assert.equal(HCR.t(undefined), '');
assert.equal(HCR.orderStateLabel({state: 'keep', state_label: '保留待复核'}), 'Retained for review');
assert.equal(HCR.orderClosestLabel({closest: '未检出'}), 'Not detected');
config.language = 'zh'; config.messages = {};
vm.runInThisContext(fs.readFileSync('static/js/i18n.js', 'utf8'));
global.HCR = window.HCR;
vm.runInThisContext(fs.readFileSync('static/js/app.js', 'utf8'));
assert.equal(HCR.locale, 'zh-CN');
assert.equal(HCR.orderStateLabel({state: 'keep', state_label: 'Retained for review'}), '保留待复核');
assert.equal(HCR.orderClosestLabel({closest: 'Not detected'}), '未检出');
'''
        result = subprocess.run(['node', '-e', script], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
