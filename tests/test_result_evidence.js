const test = require('node:test');
const assert = require('node:assert/strict');
const {closest, compare, targetMismatchPositions, alignmentColumns, neighbor} = require('../static/js/result-evidence.js');
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

test('mismatch bars agree with sense-strand sequence evidence for unequal L/R arms', () => {
  // Lj_AGL62-28: asymmetric mismatches expose a reversal within either arm.
  const arms = [
    {target: 'CTCTCTCTTCCCGTCCT', matched: 'TTCTCTCTTCCCTCCCT', positions: [4, 5, 17]},
    {target: 'CCTTCTCGCTTCTAGTTC', matched: 'CCTTCTCCCTTCTCCTTC', positions: [4, 5, 11]}
  ];
  const bar = [], sequence = [];
  arms.forEach(({target, matched, positions}) => {
    const converted = targetMismatchPositions(target.length, positions);
    const differences = [...matched].flatMap((base, i) => base === target[i] ? [] : [i + 1]);
    assert.deepEqual(converted, differences);
    bar.push(...[...target].map((_, i) => converted.includes(i + 1)));
    sequence.push(...[...matched].map((base, i) => base !== target[i]));
  });
  assert.deepEqual(bar, sequence);
  assert.equal(bar.filter(Boolean).length, 6);
});

test('sense-strand conversion handles end bases and empty hits without changing stored positions', () => {
  const positions = Object.freeze([1, 3, 20]);
  assert.deepEqual(targetMismatchPositions(20, positions), [1, 18, 20]);
  assert.deepEqual(positions, [1, 3, 20]);
  assert.deepEqual(targetMismatchPositions(17, []), []);
  assert.deepEqual(targetMismatchPositions(18), []);
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
