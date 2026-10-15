import io
import json
import unittest
from unittest.mock import patch

import app as web


FONT_AWESOME = 'font-awesome/5.15.4/css/all.min.css'


class LegacyPageTests(unittest.TestCase):
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

    def test_filter_page_and_nav_link(self):
        self.assertEqual(self.client.get('/filter').status_code, 200)
        self.assertIn('href="/filter"', self.client.get('/').get_data(as_text=True))

    def test_bridge_pages_render_download_scripts(self):
        for path in ['/bridge', '/bridge/query', '/bridge/sequence', '/bridge/analyze']:
            with self.subTest(path=path):
                page = self.client.get(path).get_data(as_text=True)
                self.assertIn('function downloadResults', page)
        self.assertNotIn('toggleQueryForm', self.client.get('/bridge/query').get_data(as_text=True))

    def test_icon_library_loaded_on_pages_using_icons(self):
        for path in ['/', '/formatting', '/probe_recommend']:
            with self.subTest(path=path):
                page = self.client.get(path).get_data(as_text=True)
                self.assertEqual(page.count(FONT_AWESOME), 1)


if __name__ == '__main__':
    unittest.main()
