import io
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

import openpyxl

from src.dual_probe import probe_selector as selector


L_ADAPTER = 'TCTCGTTTTCTGAACTTAGCCGGTATCAAG'
R_ADAPTER = 'CTGTTTAAGATTAGACTCTCGATCAGCCCA'


def write_gene(folder, name, bp_id='BP_0079', pairs=12, seed=7, scores=True):
    """A results folder shaped like the dual-probe pipeline output (BED 0-based, CSV 1-based)."""
    rng = random.Random(seed)
    folder = Path(folder)
    folder.mkdir(parents=True)
    bed, dual, arms, detailed = [], ['Probe ID,Sequence'], ['Set_ID,Primer_Type,Sequence,Start,End,Length,Tm,GC_Content,Has_GC_Junction(Left only)'], ['双探针设计结果报告', '=' * 80, '']
    position = 10
    for number in range(1, pairs + 1):
        position += rng.randrange(4, 60)
        left = ''.join(rng.choice('ACGT') for _ in range(rng.randrange(17, 21)))
        right = ''.join(rng.choice('ACGT') for _ in range(rng.randrange(17, 21)))
        l_start, l_end = position, position + len(left)
        r_start, r_end = l_end, l_end + len(right)
        bed += [f'{name}\t{l_start}\t{l_end}\tL-{number}\t0\t+', f'{name}\t{r_start}\t{r_end}\tR-{number}\t0\t+']
        dual += [f'{name}-{number}-L,{left}{L_ADAPTER}', f'{name}-{number}-R-{bp_id},{R_ADAPTER}{right}']
        arms += [f'{name}_{number},Left,{left},{l_start + 1},{l_end},{len(left)},{rng.uniform(47, 53):.1f},50.0,No',
                 f'{name}_{number},Right,{right},{r_start + 1},{r_end},{len(right)},{rng.uniform(47, 53):.1f},50.0,-']
        detailed += [f'探针组合 {number}:', '-' * 80, '', '探针组合评估:',
                     f'探针组合质量得分: {rng.choice([0, 5, 10, 15, 35, 60]):.1f} (分数越低越好)', '', '=' * 80, '']
        position = r_end
    (folder / 'dual_probe.bed').write_text('\n'.join(bed) + '\n')
    (folder / 'dual_probe_dual.csv').write_text('\n'.join(dual) + '\n')
    (folder / 'dual_probe.csv').write_text('\n'.join(arms) + '\n')
    if scores:
        (folder / 'dual_probe_detailed.txt').write_text('\n'.join(detailed) + '\n', encoding='utf-8')
    (folder / f'{name}_filtered.json').write_text('{"probe_sets": []}')


def sheet_rows(data):
    sheet = openpyxl.load_workbook(io.BytesIO(data) if isinstance(data, bytes) else data).active
    return sheet.title, [tuple(row) for row in sheet.iter_rows(values_only=True)]


# Produced by the original ProbeSelector_Share/probe_filter.py (with only the Set_ID underscore
# fix) for write_gene(<dir>/'ABC_1', 'ABC_1') and TPM 12 → 7 pairs. Any change here is a change
# in the algorithm, which must stay identical to the reference script.
GOLDEN_ABC_1_TPM_12 = [
    (None, None, 'ABC_1-1-L', 'TAAAGACAATTACATAACTCTCGTTTTCTGAACTTAGCCGGTATCAAG', 48, 35, None),
    (None, None, 'ABC_1-1-R-BP_0079', 'CTGTTTAAGATTAGACTCTCGATCAGCCCATACACGTCAGCACGAAA', 47, 35, None),
    (None, None, 'ABC_1-3-L', 'ACGCCTTTACTTGCTGTGTCTCTCGTTTTCTGAACTTAGCCGGTATCAAG', 50, 60, None),
    (None, None, 'ABC_1-3-R-BP_0079', 'CTGTTTAAGATTAGACTCTCGATCAGCCCAACCCCATCGGACTGGCAT', 48, 60, None),
    (None, None, 'ABC_1-4-L', 'TTTATTACACTCAGAAACAGTCTCGTTTTCTGAACTTAGCCGGTATCAAG', 50, 10, None),
    ('ABC_1', None, 'ABC_1-4-R-BP_0079', 'CTGTTTAAGATTAGACTCTCGATCAGCCCAACTCGGGTAATTTTGAC', 47, 10, None),
    ('GENE1', 'ABC_1_7X_BP79', 'ABC_1-6-L', 'TGAATCTCTGATTTACCTCTCGTTTTCTGAACTTAGCCGGTATCAAG', 47, 35, None),
    (12, None, 'ABC_1-6-R-BP_0079', 'CTGTTTAAGATTAGACTCTCGATCAGCCCAACTCTGCCAAACTCCAGC', 48, 35, None),
    (None, None, 'ABC_1-8-L', 'TACGACGCGCTCATTCTCTCGTTTTCTGAACTTAGCCGGTATCAAG', 46, 5, 'L->R C'),
    (None, None, 'ABC_1-8-R-BP_0079', 'CTGTTTAAGATTAGACTCTCGATCAGCCCACTTGTCGGAGAGTTATGGC', 49, 5, None),
    (None, None, 'ABC_1-10-L', 'TCGGCGTCGGAGAAACTCTTCTCGTTTTCTGAACTTAGCCGGTATCAAG', 49, 15, 'R->L T'),
    (None, None, 'ABC_1-10-R-BP_0079', 'CTGTTTAAGATTAGACTCTCGATCAGCCCATTTGCCGCCTGACAAG', 46, 15, None),
    (None, None, 'ABC_1-12-L', 'TGTCGCATCACAAACGATTCTCGTTTTCTGAACTTAGCCGGTATCAAG', 48, 10, None),
    (None, None, 'ABC_1-12-R-BP_0079', 'CTGTTTAAGATTAGACTCTCGATCAGCCCAAACTGATAAATGAGCCCTTT', 50, 10, None),
]


class AlgorithmTests(unittest.TestCase):
    def test_tpm_tiers(self):
        cases = {0: 10, 0.99: 10, 1: 9, 4.9: 9, 5: 8, 9.9: 8, 10: 7, 29.9: 7, 30: 6, 49.9: 6, 50: 5, 99.9: 5, 100: 4, 1e6: 4}
        for tpm, pairs in cases.items():
            self.assertEqual(selector.get_target_pairs(tpm), pairs, tpm)
        self.assertEqual(selector.tpm_rule_rows()[0], ('< 1', 10))
        self.assertEqual(selector.tpm_rule_rows()[-1], ('≥ 100', 4))

    def test_corrections_follow_the_reference_rules(self):
        # Full synthesis oligos, as the reference script passes them.
        l_arm, r_arm = 'CATGCATTAGCATGCA', 'TTGACTAGCATGCATA'
        new_l, new_r, note = selector.correct_probe_sequence(l_arm + L_ADAPTER, R_ADAPTER + r_arm)
        self.assertEqual(note, 'R->L A')
        self.assertEqual(new_l, 'A' + l_arm + L_ADAPTER)
        self.assertEqual(new_r, R_ADAPTER + r_arm[:-1])
        self.assertEqual(selector.correct_probe_sequence('ATGC', 'GGCC'), ('ATGC', 'GGCC', ''))
        self.assertEqual(selector.correct_probe_sequence('GCGATT', 'ATTGCG')[2], '跳过(L探针5端3个以上C/G且R端CG)')
        self.assertEqual(selector.correct_probe_sequence('GCGATTACGATTAGCA', 'ATTGCATTAGCAAGCT')[2], '跳过(TM>3)')
        self.assertEqual(selector.correct_probe_sequence('GCNATTACGATTAGCA', 'ATTGCATTAGCAAGCA')[2], '跳过(非A/T)')

    def test_set_ids_with_underscores(self):
        with tempfile.TemporaryDirectory() as folder:
            write_gene(Path(folder) / 'ABC_1', 'ABC_1', pairs=3)
            parsed = selector.parse_dual_probe_csv(str(Path(folder) / 'ABC_1' / 'dual_probe.csv'))
        self.assertEqual(sorted(parsed), [1, 2, 3])
        self.assertTrue(parsed[2]['L']['tm'] and parsed[2]['R']['sequence'])

    def test_workbook_matches_reference_output(self):
        with tempfile.TemporaryDirectory() as folder:
            write_gene(Path(folder) / 'ABC_1', 'ABC_1')
            table = Path(folder) / 'genes.csv'
            table.write_text('Name,GeneID,BP_ID,TPM\nABC_1,GENE1,BP_0079,12\n')
            results, warnings = selector.run(str(table), folder, str(Path(folder) / 'out.xlsx'))
            title, rows = sheet_rows(str(Path(folder) / 'out.xlsx'))
        self.assertEqual(title, 'Selected Probes')
        self.assertEqual(rows[0], ('GeneInfo', 'ProbeName', 'ProbeID', 'Sequence', 'Length', 'Score', 'Correction'))
        self.assertEqual(rows[1:], GOLDEN_ABC_1_TPM_12)
        self.assertEqual((results[0]['TargetCount'], results[0]['SelectedPairs']), (7, 7))
        self.assertEqual(warnings, [])

    def test_table_problems_are_reported_not_hidden(self):
        with self.assertRaisesRegex(selector.ProbeSelectorError, 'GeneID, BP_ID'):
            selector.parse_tpm_table('Name,TPM\nA,1\n')
        data, warnings = selector.parse_tpm_table('﻿Name,GeneID,BP_ID,TPM\nA,G1,BP_0001,x\nA,G1,BP_0001,3\nB,G1,BP_0002,\n')
        self.assertEqual(data['A']['TPM'], 3.0)   # last row wins, as in the reference script
        self.assertEqual(data['B']['TPM'], 0.0)   # unparsable TPM counts as 0
        codes = [item['code'] for item in warnings]
        self.assertEqual(codes.count('tpm_not_number'), 2)
        self.assertIn('duplicate_name', codes)
        self.assertIn('duplicate_gene_id', codes)
        with tempfile.TemporaryDirectory() as folder:
            write_gene(Path(folder) / 'A', 'A', pairs=2)
            write_gene(Path(folder) / 'Other', 'Other', pairs=2)
            results, warnings = selector.collect_gene_results(selector.discover_gene_dirs(folder), data)
        self.assertEqual([r['GeneName'] for r in results], ['A'])
        # The table says BP_0001 but the folder was designed with BP_0079: R probes cannot be found.
        self.assertEqual(sorted(w['code'] for w in warnings),
                         ['fewer_pairs', 'folder_missing', 'folder_not_in_table', 'missing_sequences'])
        self.assertEqual([p['ProbeID'] for p in results[0]['Probes']], ['A-1-L', 'A-2-L'])


class WebTests(unittest.TestCase):
    def setUp(self):
        import app as web
        self.web = web
        self.tmp = tempfile.TemporaryDirectory()
        self.uploads = Path(self.tmp.name) / 'uploads'
        self.config = patch.dict(web.app.config, UPLOAD_FOLDER=str(self.uploads))
        self.config.start()
        self.client = web.app.test_client()
        self.job = 'a' * 32
        job = self.uploads / 'split_results' / self.job
        write_gene(job / 'output' / 'ABC_1', 'ABC_1')
        write_gene(job / 'output' / 'XYZ', 'XYZ', pairs=3, seed=3, scores=False)
        tasks = [dict(name='ABC_1', gene_id='GENE1', bridge_id='BP_0079', sequence_length=900),
                 dict(name='XYZ', gene_id='GENE2', bridge_id='BP_0079', sequence_length=300)]
        (job / 'request.json').write_text(json.dumps(dict(
            id=self.job, status='completed', mode='batch', created_at='2026-11-01T00:00:00+00:00',
            queued_at='2026-11-01T00:00:00+00:00', completed_at='2026-11-01T00:01:00+00:00',
            reference=None, summary=dict(candidates=15, selected=8), tasks=tasks)))
        self.table = 'Name,GeneID,BP_ID,TPM\nABC_1,GENE1,BP_0079,12\nXYZ,GENE2,BP_0079,200\n'

    def tearDown(self):
        self.config.stop()
        self.tmp.cleanup()

    def post(self, **data):
        return self.client.post('/split/selector?lang=zh', data=data, content_type='multipart/form-data')

    def test_page_lists_completed_jobs_and_preselects_one(self):
        page = self.client.get('/split/selector?job=' + self.job).text
        config = json.loads(page.split('id="selector-config">')[1].split('</script>')[0])
        self.assertEqual(config['selected'], self.job)
        self.assertEqual([t['name'] for t in config['jobs'][0]['tasks']], ['ABC_1', 'XYZ'])
        self.assertIn('≥ 100', page)

    def test_job_source_downloads_the_reference_workbook(self):
        response = self.post(source='job', job_id=self.job, gene_text=self.table)
        self.assertEqual(response.status_code, 200)
        self.assertIn('selected_probes_v3.xlsx', response.headers['Content-Disposition'])
        title, rows = sheet_rows(response.data)
        self.assertEqual(rows[1:15], GOLDEN_ABC_1_TPM_12)
        self.assertEqual([r[2] for r in rows[15:]], ['XYZ-1-L', 'XYZ-1-R-BP_0079', 'XYZ-2-L', 'XYZ-2-R-BP_0079',
                                                     'XYZ-3-L', 'XYZ-3-R-BP_0079'])
        self.assertEqual({r[5] for r in rows[15:]}, {'N/A'})   # no dual_probe_detailed.txt
        summary = json.loads(response.headers['X-Selection-Summary'])
        self.assertEqual((summary['genes'], summary['pairs']), (2, 10))
        self.assertEqual(summary['rows'][1][:5], ['XYZ', 'GENE2', 200.0, 4, 3])
        self.assertEqual(summary['warnings'], ['XYZ：目标 4 对，只有 3 对可选。'])

    def test_zip_source_accepts_nested_and_single_target_layouts(self):
        archive = io.BytesIO()
        with ZipFile(archive, 'w') as bundle:
            root = self.uploads / 'split_results' / self.job / 'output'
            for path in (root / 'ABC_1').iterdir():
                bundle.write(path, 'batch/ABC_1/' + path.name)
                bundle.write(path, path.name)  # single-target ZIPs keep files at the root
            bundle.writestr('../escape/dual_probe.bed', 'x')
            bundle.writestr('__MACOSX/batch/ABC_1/._dual_probe.bed', 'x')
        archive.seek(0)
        response = self.post(source='zip', results_zip=(archive, 'results.zip'),
                             gene_text='Name,GeneID,BP_ID,TPM\nABC_1,GENE1,BP_0079,12\n')
        self.assertEqual(response.status_code, 200, response.data[:300])
        self.assertEqual(sheet_rows(response.data)[1][1:], GOLDEN_ABC_1_TPM_12)
        summary = json.loads(response.headers['X-Selection-Summary'])
        self.assertEqual(summary['warnings'], ['ZIP 中有多个名为 ABC_1 的结果文件夹，只使用 /。'])
        self.assertFalse((Path(self.tmp.name) / 'escape').exists())

    def test_errors_are_explained_in_both_languages(self):
        bad = self.post(source='job', job_id=self.job, gene_text='Name,TPM\nABC_1,1\n')
        self.assertEqual(bad.status_code, 400)
        self.assertIn('GeneID, BP_ID', bad.json['error'])
        english = self.client.post('/split/selector?lang=en', data=dict(source='job', job_id=self.job,
                                                                        gene_text='Name,GeneID,BP_ID,TPM\nNope,,BP_1,1\n'))
        self.assertEqual(english.status_code, 400)
        self.assertEqual(english.json['error'].splitlines()[:2], [
            "No genes to output: the gene table's Name must match the results folder name (task name).",
            'Results folder ABC_1 is not in the gene table; skipped.'])
        self.assertEqual(self.post(source='job', job_id='b' * 32, gene_text=self.table).status_code, 404)
        self.assertEqual(self.post(source='zip', results_zip=(io.BytesIO(b'not a zip'), 'r.zip'),
                                   gene_text=self.table).json['error'], '结果文件需为 ZIP 压缩包。')
        english = self.client.post('/split/selector?lang=en', data=dict(source='job', job_id=self.job, gene_text=self.table + 'Q,,BP_1,x\n'))
        notes = json.loads(english.headers['X-Selection-Summary'])['warnings']
        self.assertIn('Q: TPM "x" is not a number and counts as 0 (10 pairs).', notes)

    def test_probe_number_recommendation_uses_the_same_rule(self):
        self.assertEqual([self.web.get_recommended_probe_number(t) for t in (0.5, 3, 7, 20, 40, 70, 150)],
                         [10, 9, 8, 7, 6, 5, 4])
        page = self.client.get('/probe_recommend').text
        self.assertIn('≥ 100', page)
        self.assertNotIn('≥ 150', page)


if __name__ == '__main__':
    unittest.main()
