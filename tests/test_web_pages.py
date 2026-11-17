import io
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

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
    '/subgenomes': '亚基因组编号',
    '/references': '参考库管理',
    '/probe_recommend': '探针数量推荐',
    '/filter': 'AI 筛选',
    '/help': '帮助',
    '/jobs': '任务记录',
    '/orders': '下单清单',
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

    def test_home_does_not_list_every_reference(self):
        """The production server has 52 references; the tile must stay one line."""
        genomes = [{'filename': f'ref{i}.fa', 'display_name': f'ref{i:02d}'} for i in range(52)]
        with patch.object(web, 'get_available_genomes', return_value=genomes):
            page = self.client.get('/').get_data(as_text=True)
        self.assertIn('<strong>52 个</strong>', page)
        self.assertIn('ref00 · ref01 等', page)
        self.assertNotIn('ref02', page)


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

    def test_formatting_keeps_the_record_after_a_removed_one(self):
        text, stats = web.format_fasta_sequences('>bad\nACGTACGTAC!\n>kept1\nACGTACGT\n>kept2\nTTTTACGT\n',
                                                 deduplicate=False)
        self.assertEqual(text, '>kept1\nACGTACGT\n>kept2\nTTTTACGT')
        self.assertEqual((stats['original_count'], stats['removed_count'], stats['formatted_count']), (3, 1, 2))
        self.assertEqual(stats['removed_sequences'][0][0], 'bad')

    def test_recommendation_reports_unusable_input_instead_of_an_empty_table(self):
        def post(cds, counts):
            return self.client.post('/probe_recommend', content_type='multipart/form-data', data={
                'cds_file': (io.BytesIO(cds), 'cds.fa'), 'counts_file': (io.BytesIO(counts), 'counts.txt')})
        with tempfile.TemporaryDirectory() as directory, patch.dict(web.app.config, UPLOAD_FOLDER=directory):
            good = post(b'>g1\nATGCATGCAT\n>g2\nATGC\n', b'g1\t10\ng2\t3\n')
            self.assertEqual(good.status_code, 200)
            self.assertEqual(good.get_data(as_text=True).splitlines()[0], 'gene_id,tpm,probe_number')
            self.assertEqual(len(good.get_data(as_text=True).splitlines()), 3)
            good.close()
            for cds, counts in [(b'>g1\nATGC\n', b'other\t10\n'), (b'>g1\nATGC\n', b'g1\t0\n'),
                                (b'\x1f\x8b\xff\xfe', b'g1\t10\n')]:
                response = post(cds, counts)
                self.assertEqual(response.status_code, 400)
                self.assertTrue(response.get_json()['error'])
            self.assertEqual(len(list(Path(directory).iterdir())), 1)  # only the successful table was written


class TripletStorageTests(unittest.TestCase):
    """Legacy triplet requests keep one ZIP and never write outside their own results folder."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.uploads = Path(self.tmp.name) / 'site' / 'uploads'
        self.uploads.mkdir(parents=True)
        self.config = patch.dict(web.app.config, UPLOAD_FOLDER=str(self.uploads))
        self.config.start()
        self.client = web.app.test_client()
        self.sequence = web.cdna_dict['AT1G01010']

    def tearDown(self):
        self.config.stop()
        self.tmp.cleanup()

    def files(self):
        return sorted(str(path.relative_to(self.tmp.name)) for path in Path(self.tmp.name).rglob('*') if path.is_file())

    def test_single_rejects_path_like_names_and_leaves_no_working_directory(self):
        form = {'bp_id': 'BP_1101', 'geneID': '', 'seq': self.sequence, 'sequenceType': 'cdna'}
        for name in ['../../../escaped', 'a/b', '..', '']:
            with self.subTest(name=name):
                response = self.client.post('/triplet', data=dict(form, name=name))
                self.assertEqual(response.status_code, 400)
                self.assertIn('任务名称', response.get_json()['error'])
        self.assertEqual(self.files(), [])
        self.assertEqual(self.client.post('/triplet', data=dict(form, name='t', bp_id='BP_9999')).status_code, 400)
        self.assertEqual(list(self.uploads.iterdir()), [])

        response = self.client.post('/triplet', data=dict(form, name='pasted'))
        self.assertEqual(response.status_code, 200)
        with ZipFile(io.BytesIO(response.get_data())) as archive:
            names = archive.namelist()
        response.close()
        self.assertIn('target.fasta', names)  # a pasted sequence has no gene ID to name the file
        self.assertIn('pasted_filtered.json', names)
        results = sorted(path.name for path in (self.uploads / 'triplet_results').iterdir())
        self.assertEqual(len(results), 1)
        self.assertRegex(results[0], r'^triplet_results_[a-f0-9]{32}\.zip$')
        self.assertEqual(len(self.files()), 1)

    def test_targets_with_ambiguous_bases_are_refused_not_crashed(self):
        gene = next(k for k, v in web.cdna_dict.items() if set(v) - set('ACGT'))
        single = self.client.post('/triplet', data={'name': 't', 'bp_id': 'BP_1101', 'geneID': gene, 'seq': '', 'sequenceType': 'cdna'})
        self.assertEqual(single.status_code, 400)
        self.assertIn('A/T/C/G', single.get_json()['error'])
        batch = self.client.post('/triplet/batch', content_type='multipart/form-data', data={
            'batch_file': (io.BytesIO(f't,BP_1101,{gene}\n'.encode()), 'tasks.csv')})
        self.assertEqual(batch.status_code, 400)
        self.assertEqual(self.files(), [])

    def test_batch_skips_path_like_names_and_cleans_up_after_failures(self):
        def post(text):
            return self.client.post('/triplet/batch', content_type='multipart/form-data', data={
                'batch_file': (io.BytesIO(text.encode()), 'tasks.csv')})
        rejected = post('../../../escaped,BP_1101,AT1G01010\n')
        self.assertEqual(rejected.status_code, 400)
        self.assertIn('任务名称', rejected.get_json()['error'])
        self.assertEqual(self.files(), [])
        self.assertEqual(list((self.uploads / 'triplet_results').iterdir()), [])

        with patch.object(web, 'create_triplet_probe', side_effect=RuntimeError('engine failed')), \
             self.assertLogs(web.app.logger, 'ERROR'):
            self.assertEqual(post('ok,BP_1101,AT1G01010\n').status_code, 500)
        self.assertEqual(list((self.uploads / 'triplet_results').iterdir()), [])

        mixed = post('ok,BP_1101,AT1G01010\n../bad,BP_1101,AT1G01010\n')
        self.assertEqual(mixed.status_code, 200)
        self.assertEqual((mixed.headers['X-Tasks-Completed'], mixed.headers['X-Tasks-Skipped']), ('1', '1'))
        with ZipFile(io.BytesIO(mixed.get_data())) as archive:
            names = archive.namelist()
            self.assertIn('ok/AT1G01010.fasta', names)
            self.assertIn('skipped_tasks.tsv', names)
            self.assertIn('任务名称', archive.read('skipped_tasks.tsv').decode())
        mixed.close()
        self.assertTrue(all(name.startswith('ok/') or name == 'skipped_tasks.tsv' for name in names))
        self.assertEqual(len(self.files()), 1)


if __name__ == '__main__':
    unittest.main()
