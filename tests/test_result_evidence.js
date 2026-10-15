const test = require('node:test');
const assert = require('node:assert/strict');
const {closest, compare, alignmentColumns, neighbor} = require('../static/js/result-evidence.js');
const same = {subject_id: 'BjuB', relation: 'same_species_other', subject_subgenome: 'B', left_mismatches: 0, right_mismatches: 0};
const cross = {subject_id: 'BoC', relation: 'cross_species', subject_subgenome: 'C', left_mismatches: 1, right_mismatches: 2};

test('polyploid view prioritizes cross-species while ordinary view retains mismatch ordering', () => {
  const probe = {specificity: {paired_matches: [same, cross]}};
  assert.equal(closest(probe, true, 'A'), cross);
  assert.equal(closest(probe, false, 'A'), same);
  assert.ok(compare(cross, same, true, 'A') < 0);
});

test('within a tier the better worst-arm mismatch count wins; targets are never the nearest off-target', () => {
  const better = {...cross, subject_id: 'BoBetter', left_mismatches: 1, right_mismatches: 1};
  const target = {...same, relation: 'target'};
  assert.equal(closest({specificity: {paired_matches: [cross, better, target]}}, true, 'A'), better);
  assert.equal(closest({specificity: {paired_matches: [target]}}, true, 'A'), null);
});

test('gapped alignment preserves target coordinates, reverse subject coordinates, and uncovered flanks', () => {
  const columns = alignmentColumns('ACGTACGT', {query_start: 2, query_aligned: 'CG-TACG', subject_aligned: 'CGTT-CG', subject_start: 12, strand: '-'});
  assert.equal(columns.map(c => c.query).join(''), 'ACG-TACGT');
  assert.equal(columns.map(c => c.subject).join(''), '·CGTT-CG·');
  assert.deepEqual(columns.map(c => c.pos), [1, 2, 3, null, 4, 5, 6, 7, 8]);
  assert.deepEqual(columns.map(c => c.subjectPos), [null, 12, 11, 10, 9, null, 8, 7, null]);
  assert.equal(columns[3].anchor, 3.5);
});

test('evidence navigation follows the filtered/sorted rows and stops at either end', () => {
  const a = {id: 7}, b = {id: 2}, c = {id: 3};
  assert.equal(neighbor([b, a, c], a, -1), b);
  assert.equal(neighbor([b, a, c], a, 1), c);
  assert.equal(neighbor([a, c], b, 1), null);
  assert.equal(neighbor([a, c], a, -1), null);
  assert.equal(neighbor([a, c], c, 1), null);
});
