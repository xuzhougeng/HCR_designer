import hashlib
import json
from pathlib import Path
import random
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from flask import Flask

from scripts.split_web import register_split_routes
from src.common.blast_runtime import blast_executable
from src.common.homolog_alignment import AlignmentUnavailable, find_homolog, parse_alignments


class HomologAlignmentTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        rng = random.Random(34)
        self.query = ''.join(rng.choice('ACGT') for _ in range(240))
        mutated = list(self.query)
        for index in (10, 80, 160):
            mutated[index] = 'A' if mutated[index] != 'A' else 'C'
        self.reference = self.root / 'reference.fa'
        self.reference.write_text('>BjuTarget\n' + self.query + '\n>BjuOther\n' + self.query
                                  + '\n>BoOther\n' + ''.join(mutated) + '\n>Unmapped\n' + self.query + '\n')
        self.digest = hashlib.sha256(self.reference.read_bytes()).hexdigest()
        self.annotation = self.root / 'subgenomes.tsv'
        self.annotation.write_text('BjuTarget\tA\nBjuOther\tB\nBoOther\tC\n')

    def tearDown(self):
        self.tmp.cleanup()

    def find(self):
        return find_homolog(self.reference, self.annotation, self.digest, 'BjuTarget', self.query)

    def test_real_blast_prefers_annotated_cross_species_over_perfect_same_species(self):
        try:
            blast_executable('blastn')
        except FileNotFoundError:
            self.skipTest('BLAST+ unavailable')
        result = self.find()
        self.assertEqual(result['subject_id'], 'BoOther')
        self.assertEqual(result['target_group'], 'A')
        self.assertEqual(result['subject_group'], 'C')
        self.assertEqual(result['searched_sequences'], 2)  # target and unannotated exact copy excluded
        self.assertEqual(result['query_coverage'], 100)
        self.assertEqual(result['identity'], 98.75)
        with patch('src.common.homolog_alignment.subprocess.run') as blast:
            self.assertEqual(self.find(), result)
            blast.assert_not_called()

    def test_same_group_and_unannotated_sequences_are_not_comparison_candidates(self):
        self.annotation.write_text('BjuTarget\tA\nBjuOther\tA\n')
        with patch('src.common.homolog_alignment.subprocess.run') as blast:
            result = self.find()
        self.assertEqual(result['status'], 'no_match')
        self.assertEqual(result['searched_sequences'], 0)
        blast.assert_not_called()

    def test_other_subgenome_fallback_uses_uploaded_labels_and_no_hit_is_explicit(self):
        self.annotation.write_text('BjuTarget\tcustom_target\nBjuOther\tcustom_other\n')
        output = 's0\t1\t240\t1\t240\t400\t1e-80\t100\t' + self.query + '\t' + self.query
        with patch('src.common.homolog_alignment.blast_executable', return_value='blastn'), \
                patch('src.common.homolog_alignment.subprocess.run', return_value=subprocess.CompletedProcess([], 0, output)):
            result = self.find()
        self.assertEqual(result['relation'], 'other_subgenome')
        self.assertEqual(result['subject_id'], 'BjuOther')
        self.assertEqual(result['subject_group'], 'custom_other')
        self.annotation.write_text('BjuTarget\tA\nBoOther\tC\n')
        with patch('src.common.homolog_alignment.blast_executable', return_value='blastn'), \
                patch('src.common.homolog_alignment.subprocess.run', return_value=subprocess.CompletedProcess([], 0, '')):
            result = self.find()
        self.assertEqual(result['status'], 'no_match')
        self.assertNotIn('subject_id', result)

    def test_reference_content_change_is_rejected_even_with_same_size_and_mtime(self):
        import os
        # Prime the verified content cache without needing BLAST.
        self.annotation.write_text('BjuTarget\tA\n')
        self.find()
        stat = self.reference.stat()
        self.reference.write_bytes(self.reference.read_bytes().replace(b'AC', b'TG'))
        os.utime(self.reference, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        with self.assertRaisesRegex(AlignmentUnavailable, '不一致'):
            self.find()

    def test_unknown_target_annotation_and_timeout_have_explicit_errors(self):
        self.annotation.write_text('BoOther\tC\n')
        with self.assertRaisesRegex(AlignmentUnavailable, '未标注目标'):
            self.find()
        self.annotation.write_text('BjuTarget\tA\nBoOther\tC\n')
        with patch('src.common.homolog_alignment.blast_executable', return_value='blastn'), \
                patch('src.common.homolog_alignment.subprocess.run', side_effect=subprocess.TimeoutExpired('blastn', 90)), \
                self.assertRaisesRegex(AlignmentUnavailable, '超时'):
            self.find()

    def test_parser_preserves_insertions_deletions_and_reverse_subject_coordinates(self):
        subjects = {'s0': dict(subject_id='BoOther', subject_group='B', subject_length=20, relation='cross_species')}
        result = parse_alignments('s0\t2\t7\t12\t7\t42\t1e-8\t71.4\tCG-TACG\tCGTT-CG', subjects, 'ACGTACGT')[0]
        self.assertEqual(result['strand'], '-')
        self.assertEqual(result['query_coverage'], 75)
        self.assertEqual(result['query_aligned'], 'CG-TACG')
        self.assertEqual((result['subject_start'], result['subject_end']), (12, 7))
        with self.assertRaises(ValueError):
            parse_alignments('s0\t1\t7\t12\t7\t42\t1e-8\t71.4\tCG-TACG\tCGTT-CG', subjects, 'ACGTACGT')

    def test_http_uses_frozen_job_annotation_and_validates_mode_and_task(self):
        identifier = 'a' * 32
        directory = self.root / 'split_results' / identifier
        directory.mkdir(parents=True)
        frozen = directory / 'subgenomes.tsv'
        frozen.write_bytes(self.annotation.read_bytes())
        manifest = dict(id=identifier, status='completed', process_multiple_sets=True,
                        parameters={'target_mode': 'subgenome'}, reference_path=str(self.reference),
                        reference={'sha256': self.digest}, subgenome_map='/previous/server/path.tsv',
                        tasks=[{'gene_id': 'BjuTarget', 'sequence': self.query}])
        path = directory / 'request.json'
        path.write_text(json.dumps(manifest))
        app = Flask(__name__)
        app.config.update(TESTING=True, UPLOAD_FOLDER=str(self.root), SPLIT_JOBS_INLINE=True)
        register_split_routes(app, {}, {}, {}, {}, lambda: [])
        client = app.test_client()
        url = '/split/results/' + identifier + '/alignment'
        with patch('scripts.split_workflow.find_homolog', return_value={'status': 'no_match'}) as find:
            self.assertEqual(client.get(url).status_code, 200)
            self.assertEqual(find.call_args.args[1], frozen.resolve())
            self.assertEqual(find.call_args.args[2], self.digest)
            self.assertEqual(client.get(url + '?task=-1').status_code, 400)
            self.assertEqual(client.get(url + '?task=1').status_code, 404)
            self.assertEqual(client.get(url + '?task=x').status_code, 400)
        manifest.update(process_multiple_sets=False, parameters={'target_mode': 'species'})
        path.write_text(json.dumps(manifest))
        self.assertEqual(client.get(url).status_code, 409)
        self.assertEqual(client.get('/split/results/' + 'b' * 32 + '/alignment').status_code, 404)


if __name__ == '__main__':
    unittest.main()
