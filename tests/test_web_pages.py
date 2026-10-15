import io
import json
import re
import unittest
from unittest.mock import patch

import app as web


PAGES = {
    '/': None,
    '/split': '双探针设计',
    '/split/batch': '双探针设计',
    '/triplet': '三探针设计',
    '/triplet/batch': '三探针设计',
    '/bridge/list': '浏览探针库',
    '/bridge/query': '冲突检查',
    '/bridge/sequence': '冲突检查',
    '/bridge/analyze': '冲突检查',
    '/bridge': '推荐新探针',
    '/bridge/validate': '合规验证',
    '/formatting': 'FASTA 格式化',
    '/probe_recommend': '探针数量推荐',
    '/filter': 'AI 筛选',
    '/help': '帮助',
}


class WorkbenchShellTests(unittest.TestCase):
    def setUp(self):
        self.client = web.app.test_client()

    def test_every_page_uses_shell_and_marks_its_nav_entry(self):
        for path, label in PAGES.items():
            with self.subTest(path=path):
                response = self.client.get(path)
                self.assertEqual(response.status_code, 200)
                page = response.get_data(as_text=True)
                self.assertIn('class="sidebar"', page)
                self.assertIn('/static/js/app.js', page)
                current = re.findall(r'class="side-link"[^>]*aria-current="page"[^>]*>.*?</svg>([^<]+)', page, re.S)
                self.assertEqual([c.strip() for c in current], [label] if label else [])
                # The intranet server may have no internet access: no third-party assets.
                self.assertNotRegex(page, r'<(link|script)[^>]+(href|src)="https?://')

    def test_home_reports_library_size(self):
        page = self.client.get('/').get_data(as_text=True)
        self.assertIn('{:,} 条'.format(len(web.bridge_seq_dict)), page)


class LegacyToolTests(unittest.TestCase):
    def setUp(self):
        self.client = web.app.test_client()

    def test_filter_upload_reaches_llm_handler(self):
        """POST /filter used to be swallowed by a duplicate route that ignored the upload."""
        filtered = {'report_title': '筛选报告', 'probe_sets': []}
        with patch.object(web, 'filter_probe_by_llm', return_value=(filtered, '推理过程示例')) as llm:
            response = self.client.post('/filter', content_type='multipart/form-data', data={
                'probe_file': (io.BytesIO(json.dumps({'probe_sets': []}).encode()), 'probes.json')})
        self.assertEqual(response.status_code, 200)
        llm.assert_called_once_with({'probe_sets': []})
        page = response.get_data(as_text=True)
        self.assertIn('筛选报告', page)
        self.assertIn('推理过程示例', page)

    def test_bridge_results_offer_download(self):
        bp_id, sequence = next(iter(web.bridge_seq_dict.items()))
        responses = [
            self.client.post('/bridge/query', data={'bp_id': bp_id, 'kmer': 5}),
            self.client.post('/bridge/sequence', data={'sequences': sequence, 'kmer': 9}),
            self.client.post('/bridge/analyze', data={'input_content': f'{bp_id}\n{sequence}', 'kmer': 9}),
        ]
        for response in responses:
            with self.subTest(url=response.request.path):
                self.assertEqual(response.status_code, 200)
                self.assertIn('id="download-button"', response.get_data(as_text=True))

    def test_triplet_errors_are_json(self):
        response = self.client.post('/triplet', data={'name': 't', 'bp_id': 'BP_9999', 'geneID': '', 'seq': 'ACGT'})
        self.assertEqual(response.status_code, 400)
        self.assertIn('BP_9999', response.get_json()['error'])
        response = self.client.post('/triplet/batch', data={})
        self.assertEqual(response.status_code, 400)
        self.assertIn('CSV', response.get_json()['error'])

    def test_triplet_batch_reports_rows_it_cannot_design(self):
        csv_text = 'task1,BP_9999,AT1G01010\ntask2,BP_1101,NOT_A_GENE\nshort,row\n'
        with patch.object(web, 'create_triplet_probe') as design:
            response = self.client.post('/triplet/batch', content_type='multipart/form-data', data={
                'batch_file': (io.BytesIO(csv_text.encode()), 'tasks.csv')})
        design.assert_not_called()
        self.assertEqual(response.status_code, 400)
        error = response.get_json()['error']
        for fragment in ['第 1 行', 'BP_9999', '第 2 行', 'NOT_A_GENE', '第 3 行']:
            self.assertIn(fragment, error)

    def test_formatting_returns_summary_header(self):
        fasta = '>AT1G01010.1 desc\nacgtacgt\n>AT1G01010.2\nACGTACGTACGT\n>bad\nACGTXZ\n'
        response = self.client.post('/formatting', content_type='multipart/form-data', data={
            'fasta_file': (io.BytesIO(fasta.encode()), 'genes.fasta'),
            'clean_headers': 'on', 'uppercase': 'on', 'filter_invalid': 'on', 'deduplicate': 'on'})
        self.assertEqual(response.status_code, 200)
        stats = json.loads(response.headers['X-Format-Stats'])
        self.assertEqual(set(stats), {'original_count', 'formatted_count', 'removed_count', 'headers_cleaned',
                                      'deduplicated_count', 'detected_type'})
        self.assertEqual(stats['original_count'], 3)
        self.assertIn('ACGTACGTACGT', response.get_data(as_text=True))


if __name__ == '__main__':
    unittest.main()
