/* Pure result-view helpers, shared with the Node regression tests. */
((root) => {
  'use strict';
  const risk = hit => [Math.max(hit.left_mismatches, hit.right_mismatches), hit.left_mismatches + hit.right_mismatches];
  const byRisk = (a, b) => risk(a)[0] - risk(b)[0] || risk(a)[1] - risk(b)[1];
  function priority(hit, group) {
    if (hit.relation === 'cross_species') return 0;
    if (group && hit.subject_subgenome && hit.subject_subgenome !== group) return 1;
    return hit.relation === 'allowed' ? 3 : 2;
  }
  function compare(a, b, polyploid, group) {
    return (polyploid ? priority(a, group) - priority(b, group) : 0) || byRisk(a, b)
      || a.subject_id.localeCompare(b.subject_id) || a.left_start - b.left_start;
  }
  function closest(probe, polyploid, group) {
    return [...(probe.specificity?.paired_matches || [])].filter(hit => hit.relation !== 'target')
      .sort((a, b) => compare(a, b, polyploid, group))[0] || null;
  }
  function targetMismatchPositions(length, probePositions = []) {
    // Each targeting arm is reverse-complementary to the target sense strand.
    // Reverse positions within each arm, keeping the L/R order on the target.
    return probePositions.map(pos => length - pos + 1).sort((a, b) => a - b);
  }
  function alignmentColumns(sequence, alignment) {
    const columns = [];
    const unaligned = pos => ({pos, anchor: pos, query: sequence[pos - 1], subject: '·', subjectPos: null, state: 'unaligned'});
    for (let pos = 1; pos < alignment.query_start; pos++) columns.push(unaligned(pos));
    let pos = alignment.query_start, subjectPos = alignment.subject_start;
    const step = alignment.strand === '-' ? -1 : 1;
    for (let i = 0; i < alignment.query_aligned.length; i++) {
      const query = alignment.query_aligned[i], subject = alignment.subject_aligned[i];
      columns.push({pos: query === '-' ? null : pos, anchor: query === '-' ? pos - .5 : pos,
        query, subject, subjectPos: subject === '-' ? null : subjectPos,
        state: query === '-' || subject === '-' ? 'gap' : query === subject ? 'match' : 'mismatch'});
      if (query !== '-') pos++;
      if (subject !== '-') subjectPos += step;
    }
    for (; pos <= sequence.length; pos++) columns.push(unaligned(pos));
    return columns;
  }
  function neighbor(rows, selected, direction) {
    const index = rows.indexOf(selected);
    return index < 0 ? null : rows[index + direction] || null;
  }
  const api = {compare, closest, targetMismatchPositions, alignmentColumns, neighbor};
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.HCRResults = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
