import csv
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest

from flask import Flask
from werkzeug.datastructures import FileStorage

from scripts.subgenome_web import register_subgenome_routes
from src.common.reference_annotation import SubgenomeAnnotation, parse_subgenome_map
from src.common.subgenome_generator import generate_subgenomes


def upload(name, content, compressed=False):
    data = content.encode('utf-8')
    return FileStorage(stream=io.BytesIO(gzip.compress(data) if compressed else data), filename=name)


class SubgenomeGeneratorTests(unittest.TestCase):
    def test_multiple_files_extract_original_ids_and_merge_one_compatible_csv(self):
        data = generate_subgenomes([
            upload('A.fasta', '>BjuVA01G1 description with spaces\nACGT\n>BjuVA02G2\nACGT\n'),
            upload('B.fna.gz', '>BjuVB01G3 note\nACGT\n', compressed=True)], 'files')
        self.assertEqual(data['summary']['subgenome_counts'], {'A': 2, 'B': 1})
        rows = list(csv.reader(io.StringIO(data['csv'])))
        self.assertEqual(rows, [['transcript_id', 'subgenome'], ['BjuVA01G1', 'A'], ['BjuVA02G2', 'A'], ['BjuVB01G3', 'B']])
        annotation = SubgenomeAnnotation(parse_subgenome_map(data['csv']))
        self.assertEqual(annotation.of('bjuvb01g3'), 'B')
        self.assertIsNone(annotation.of('description'))

    def test_explicit_labels_and_duplicates_within_same_subgenome(self):
        data = generate_subgenomes([upload('A1.fa', '>t1\nA\n>t1\nA\n'), upload('A2.fa', '>T1\nA\n>t2\nA\n')],
                                  'files', labels=['A', 'A'])
        self.assertEqual(data['summary']['duplicates'], 2)
        self.assertEqual(data['summary']['total'], 2)
        self.assertEqual(data['summary']['headers'], 4)
        self.assertEqual(data['summary']['subgenome_counts'], {'A': 2})

    def test_same_identifier_in_different_subgenomes_blocks_csv(self):
        data = generate_subgenomes([upload('A.fa', '>t1\nA\n'), upload('B.fa', '>T1\nA\n')], 'files')
        self.assertEqual(data['summary']['conflicts'], 1)
        self.assertIsNone(data['csv'])
        self.assertEqual(data['issues'][0], dict(identifier='t1', issue='conflict', labels=['A', 'B']))
        self.assertIn('t1,conflict,A|B', data['issues_csv'])

    def test_numeric_features_match_anywhere_and_ignore_header_descriptions(self):
        fasta = '>BjuVA01G1 note\nACGT\n>BjuVA12G2\nACGT\n>BjuVB09G3\nACGT\n>other contains A01\nACGT\n'
        data = generate_subgenomes([upload('combined.fa', fasta)], 'rules', rules=[
            dict(pattern='Axx', label='A'), dict(pattern='B0X', label='B')])
        self.assertEqual(data['summary']['subgenome_counts'], {'A': 2, 'B': 1})
        self.assertEqual(data['summary']['unmatched'], 1)
        self.assertEqual(data['issues'][0]['identifier'], 'other')
        self.assertNotIn('other', data['csv'])

    def test_prefix_suffix_literal_and_case_matching(self):
        fasta = '>A01_gene\nACGT\n>prefixA01\nACGT\n>a01_lower\nACGT\n>Axx_literal\nACGT\n'
        def run(**options):
            return generate_subgenomes([upload('combined.fa', fasta)], 'rules', **options)
        rule = [dict(pattern='A0X', label='A')]
        self.assertEqual(run(rules=rule, position='prefix')['summary']['assigned'], 2)
        self.assertEqual(run(rules=rule, position='suffix')['summary']['assigned'], 1)
        self.assertEqual(run(rules=rule, position='prefix', case_sensitive=True)['summary']['assigned'], 1)
        self.assertEqual(run(rules=[dict(pattern='Axx', label='A')], rule_mode='literal')['preview'],
                         [dict(identifier='Axx_literal', subgenome='A')])

    def test_overlapping_rules_conflict_only_when_labels_differ(self):
        rules = [dict(pattern='Axx', label='A'), dict(pattern='A0X', label='A')]
        data = generate_subgenomes([upload('combined.fa', '>A01\nACGT\n')], 'rules', rules=rules)
        self.assertEqual(data['summary']['assigned'], 1)
        rules[1]['label'] = 'B'
        data = generate_subgenomes([upload('combined.fa', '>A01\nACGT\n')], 'rules', rules=rules)
        self.assertEqual(data['summary']['conflicts'], 1)
        self.assertIsNone(data['csv'])

    def test_all_unmatched_and_previews_do_not_truncate_export(self):
        fasta = ''.join('>A%03d\nACGT\n' % i for i in range(125))
        data = generate_subgenomes([upload('A.fa', fasta)], 'files')
        self.assertEqual(len(data['preview']), 100)
        self.assertEqual(len(list(csv.reader(io.StringIO(data['csv'])))), 126)
        data = generate_subgenomes([upload('A.fa', fasta)], 'rules', rules=[dict(pattern='Bxx', label='B')])
        self.assertIsNone(data['csv'])
        self.assertEqual(data['summary']['unmatched'], 125)
        self.assertEqual(len(data['issues']), 100)
        self.assertEqual(len(list(csv.reader(io.StringIO(data['issues_csv'])))), 126)

    def test_invalid_or_truncated_fasta_is_rejected_without_partial_result(self):
        for name, data in [('bad.fa', b'ACGT\n>id\nA\n'), ('empty.fa', b''), ('header.fa', b'>  \nA\n'),
                           ('invalid.fa', b'>id\n\xff\n'), ('bad.csv', b'>id\nA\n'),
                           ('separator.fa', b'>id,bad\nA\n'), ('broken.fa.gz', gzip.compress(b'>id\nA\n')[:-4])]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                generate_subgenomes([FileStorage(stream=io.BytesIO(data), filename=name)], 'files')

    def test_invalid_rules_labels_modes_and_file_counts_are_rejected(self):
        options = [dict(mode='invalid'), dict(mode='files', labels=[]), dict(mode='files', labels=['A,B']),
                   dict(mode='rules', rules=[]), dict(mode='rules', rules=[dict(pattern='', label='A')]),
                   dict(mode='rules', rules=[dict(pattern='Axx', label='')]),
                   dict(mode='rules', rules=[dict(pattern='Axx', label='A')], position='invalid')]
        for values in options:
            with self.subTest(values=values), self.assertRaises(ValueError):
                generate_subgenomes([upload('A.fa', '>id\nA\n')], **values)
        with self.assertRaises(ValueError):
            generate_subgenomes([], 'files')
        with self.assertRaises(ValueError):
            generate_subgenomes([upload('A.fa', '>id\nA\n'), upload('B.fa', '>id2\nA\n')], 'rules')


class SubgenomeHttpTests(unittest.TestCase):
    def setUp(self):
        self.app = Flask(__name__)
        self.app.config['TESTING'] = True
        register_subgenome_routes(self.app)
        self.client = self.app.test_client()

    def test_file_and_rule_uploads_return_csv_without_creating_saved_results(self):
        with tempfile.TemporaryDirectory() as directory:
            self.app.config['UPLOAD_FOLDER'] = directory
            files = self.client.post('/subgenomes', content_type='multipart/form-data', data=dict(
                mode='files', labels=json.dumps(['A', 'B']), fasta_files=[
                    (io.BytesIO(b'>t1 description\nACGT\n'), 'A.fa'),
                    (io.BytesIO(gzip.compress(b'>t2\nACGT\n')), 'B.fa.gz')]))
            self.assertEqual(files.status_code, 200, files.json)
            self.assertIn('t2,B', files.json['csv'])
            rules = self.client.post('/subgenomes', content_type='multipart/form-data', data=dict(
                mode='rules', rules=json.dumps([dict(pattern='Axx', label='A')]),
                fasta_files=(io.BytesIO(b'>prefixA12\nACGT\n'), 'all.fa')))
            self.assertEqual(rules.status_code, 200, rules.json)
            self.assertEqual(rules.json['summary']['assigned'], 1)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_missing_files_and_bad_rule_json_are_reported(self):
        for data in (dict(mode='files'), dict(mode='rules', rules='{bad'), dict(mode='rules', rules='null')):
            response = self.client.post('/subgenomes', data=data)
            self.assertEqual(response.status_code, 400)
            self.assertTrue(response.json['error'])


if __name__ == '__main__':
    unittest.main()
