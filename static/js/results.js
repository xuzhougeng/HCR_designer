/* Dual-probe job page: follows a queued/running job, then shows the results workspace. */
(() => {
  'use strict';
  const byId = id => document.getElementById(id);
  const status = JSON.parse(byId('job-status').textContent);
  const STAGES = ['queued', 'reference', 'candidates', 'specificity', 'selection', 'packaging'];
  const STATE = {queued: [HCR.t("排队中"), 'info'], running: [HCR.t("计算中"), 'info'], completed: [HCR.t("已完成"), 'ok'],
                 failed: [HCR.t("失败"), 'bad'], interrupted: [HCR.t("已中断"), 'warn']};
  const REASONS = {
    cross_species_exact_pair: HCR.t("另一物种存在双臂完全配对"),
    cross_species_near_pair: HCR.t("另一物种存在阈值内的近似配对"),
    unclassified_exact_pair: HCR.t("未知物种的非目标转录本存在完全配对"),
    unclassified_near_pair: HCR.t("未知物种的非目标转录本存在阈值内的近似配对"),
    same_species_exact_pair: HCR.t("同物种非允许转录本存在完全配对"),
    same_species_near_pair: HCR.t("同物种非允许转录本存在阈值内的近似配对"),
    specificity_not_evaluated: HCR.t("配对特异性证据缺失或不完整"),
    no_exact_on_target_pair: HCR.t("未找到目标转录本上的双臂完全配对"),
    gc_junction: HCR.t("L 臂连接位点未满足当前筛选规则"),
    too_close: HCR.t("与其他已选组合间距不足"),
    selection_limit_or_rank: HCR.t("排序后未进入本次保留数量"),
    selected: HCR.t("按当前规则保留，仍需进一步复核")
  };
  const RELATIONS = {target: HCR.t("目标转录本"), allowed: HCR.t("允许共同检测"), same_species_other: HCR.t("同物种其他转录本"),
                     cross_species: HCR.t("跨物种"), unknown_species: HCR.t("物种未知")};
  const MODES = {species: HCR.t("物种模式"), transcript: HCR.t("严格模式"), allowlist: HCR.t("允许集合"), subgenome: HCR.t("亚基因组模式")};

  function el(tag, text, className) {
    const node = document.createElement(tag);
    if (text !== undefined && text !== null) node.textContent = text;
    if (className) node.className = className;
    return node;
  }
  const svgEl = (tag, attrs) => {
    const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
    Object.entries(attrs || {}).forEach(([key, value]) => node.setAttribute(key, value));
    return node;
  };
  const clock = seconds => {
    seconds = Math.max(0, Math.round(seconds));
    const m = Math.floor(seconds / 60), s = seconds % 60;
    return (m ? m + HCR.t(" 分 ") : '') + s + HCR.t(" 秒");
  };

  /* ---------------- progress ---------------- */
  function followProgress() {
    let current = status, timer = null;
    const started = () => Date.parse(current.started_at || current.queued_at || current.created_at);

    function paint() {
      const state = current.status;
      const chip = byId('status-chip');
      chip.textContent = (STATE[state] || [state])[0];
      chip.className = 'chip ' + ((STATE[state] || [])[1] || '');
      const progress = current.progress || {};
      const stage = state === 'queued' ? 'queued' : progress.stage;
      const index = STAGES.indexOf(stage);
      document.querySelectorAll('#stage-track li').forEach((item, i) => {
        item.classList.toggle('done', index > i || state === 'completed');
        item.classList.toggle('current', index === i && (state === 'running' || state === 'queued'));
      });
      const total = progress.task_total || current.tasks.length;
      const doneTasks = state === 'completed' ? total : Math.max(0, (progress.task_index || 1) - 1);
      byId('progress-bar').style.width = (total ? Math.round(100 * doneTasks / total) : 0) + '%';
      byId('progress-count').textContent = HCR.t("任务 ") + Math.min(progress.task_index || 0, total) + ' / ' + total;
      let detail;
      if (state === 'queued') detail = current.queue_position ? HCR.t("前面还有 ") + current.queue_position + HCR.t(" 个任务") : HCR.t("等待后台开始");
      else if (state === 'running') detail = (progress.task_name ? progress.task_name + ' · ' : '') + (HCR.t(progress.label) || HCR.t("计算中"));
      else if (state === 'failed') detail = HCR.t("任务失败");
      else if (state === 'interrupted') detail = HCR.t("任务被中断");
      else detail = HCR.t("已完成，正在打开结果……");
      byId('progress-detail').textContent = detail;
      const live = state === 'queued' || state === 'running';
      byId('elapsed').textContent = live ? clock((Date.now() - started()) / 1000) : '';
      const error = byId('progress-error');
      error.hidden = !current.error;
      error.textContent = current.error || '';
      byId('rerun-button').hidden = !(state === 'failed' || state === 'interrupted') || !current.parameters;
    }

    async function poll() {
      try {
        const response = await fetch(status.status_url, {headers: {Accept: 'application/json'}});
        if (response.ok) current = Object.assign({}, current, await response.json());
      } catch (error) { /* keep polling through brief network errors */ }
      paint();
      if (current.status === 'completed') { window.location.reload(); return; }
      if (current.status === 'queued' || current.status === 'running') timer = setTimeout(poll, 2000);
    }

    setInterval(() => {
      if (current.status === 'queued' || current.status === 'running') byId('elapsed').textContent = clock((Date.now() - started()) / 1000);
    }, 1000);
    byId('rerun-button').addEventListener('click', async event => {
      event.currentTarget.disabled = true;
      const response = await fetch(status.run_url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(current.parameters)});
      const data = await response.json().catch(() => ({}));
      if (!response.ok) { current.error = data.error || HCR.t("无法重新运行。"); paint(); event.currentTarget.disabled = false; return; }
      current = Object.assign({}, current, data);
      clearTimeout(timer);
      poll();
    });
    paint();
    poll();
  }

  /* ---------------- results ---------------- */
  const assessed = p => p.specificity?.status === 'completed';
  const policyExcluded = p => /_(exact|near)_pair$/.test(p.reason_code || '');
  const crossExact = p => (p.specificity?.counts?.cross_species?.[0] || 0) > 0;
  function stateOf(p) {
    if (!assessed(p)) return ['unknown', p.keep ? HCR.t("未评估 · 按质量保留") : HCR.t("特异性未评估")];
    if (p.keep) return ['keep', HCR.t("保留待复核")];
    return policyExcluded(p) ? ['excluded', HCR.t("按目标策略排除")] : ['other', HCR.t("未保留")];
  }
  const span = p => [Math.min(p.left_probe.position.start, p.right_probe.position.start),
                     Math.max(p.left_probe.position.end, p.right_probe.position.end)];
  const risk = m => [Math.max(m.left_mismatches, m.right_mismatches), m.left_mismatches + m.right_mismatches];
  const byRisk = (a, b) => risk(a)[0] - risk(b)[0] || risk(a)[1] - risk(b)[1];

  function showResults(result) {
    let task = result.tasks[0], selected = null, range = null;
    const polyploidAvailable = result.process_multiple_sets || result.parameters.target_mode === 'subgenome';
    let polyploid = Boolean(polyploidAvailable), alignment = null, alignmentColumns = [], windowStart = 1;
    let alignmentRequest = null, alignmentGeneration = 0;
    const alignmentCache = new Map();
    const compareMatches = (a, b) => HCRResults.compare(a, b, polyploid, task.subgenome);
    const closest = probe => HCRResults.closest(probe, polyploid, task.subgenome);
    byId('polyploid-controls').hidden = !polyploidAvailable;

    // summary
    const stats = byId('result-stats');
    [[HCR.t("原始候选"), 'candidates'], [HCR.t("跨物种完全配对"), 'cross_exact'], [HCR.t("按目标策略排除"), 'excluded'],
     [HCR.t("保留待复核"), 'selected'], [HCR.t("特异性未评估"), 'unassessed']].forEach(([label, key]) => {
      const tile = el('div', undefined, 'stat');
      tile.append(el('span', label), el('strong', result.tasks.reduce((sum, t) => sum + (t[key] || 0), 0).toLocaleString(HCR.locale) + HCR.t(" 对")));
      stats.append(tile);
    });
    byId('result-limit').textContent = !result.reference
      ? HCR.t("本次未指定参考库，未进行特异性评估。保留候选仅反映质量筛选，不能用于判断物种特异性。")
      : result.search_exhaustive
        ? HCR.t("双臂配对已在“正向转录本、无缺口、每臂最多 3 个替换”的范围内完整检出；插入缺失未复核。候选保留和“未检出”均不等于特异性合格。")
        : HCR.t("本次有任务未完成完整配对扫描，相关候选标记为未评估；未进行插入缺失复核。");
    const p = result.parameters, settings = byId('run-settings');
    const line = text => settings.append(el('div', text));
    if (result.reference) {
      line(result.reference.name + ' · ' + result.reference.sequence_count.toLocaleString(HCR.locale) + HCR.t(" 条序列"));
      const sha = el('div'); sha.append(HCR.t("参考 SHA256："), el('code', result.reference.sha256)); settings.append(sha);
    }
    line(HCR.t("臂长 ") + p.min_length + '–' + p.max_length + ' nt · Tm ' + p.min_tm + '–' + p.max_tm + ' °C · GC ' + p.min_gc + '–' + p.max_gc + '%');
    line(HCR.t("默认最多 ") + p.max_selected + HCR.t(" 对/基因 · 组合间距 ") + p.min_gap + HCR.t(" nt · 臂内间隔 0 nt · 连续同碱基阈值 ") + p.poly_n);
    line(HCR.t("目标判定 ") + (MODES[p.target_mode] || p.target_mode) + HCR.t(" · 非允许配对排除阈值 每臂 ≤ ") + p.exclusion_max_mismatches + HCR.t(" 错配"));
    line(HCR.t("提交 ") + new Date(result.created_at).toLocaleString(HCR.locale) + HCR.t(" · 完成 ") + new Date(result.completed_at).toLocaleString(HCR.locale));

    const picker = byId('result-task');
    result.tasks.forEach((t, i) => { const o = el('option', t.name + ' · ' + t.gene_id + ' · ' + t.bridge_id); o.value = i; picker.append(o); });
    picker.addEventListener('change', () => { task = result.tasks[Number(picker.value)]; selected = null; range = null; alignment = null; alignmentColumns = []; windowStart = 1; byId('evidence').hidden = true; render(); loadAlignment(); });
    byId('result-filter').addEventListener('change', renderTable);
    byId('result-sort').addEventListener('change', renderTable);
    byId('close-evidence').addEventListener('click', () => { selected = null; byId('evidence').hidden = true; renderTrack(); renderTable(); });
    document.addEventListener('hcr:orders', renderTable);
    byId('polyploid-view').addEventListener('change', event => {
      polyploid = event.target.checked;
      render();
      if (selected) renderEvidence(selected, false);
      loadAlignment();
    });
    byId('homolog-retry').addEventListener('click', () => loadAlignment());
    byId('homolog-back').addEventListener('click', () => moveAlignmentWindow(-80));
    byId('homolog-forward').addEventListener('click', () => moveAlignmentWindow(80));
    byId('evidence-previous').addEventListener('click', () => moveEvidence(-1));
    byId('evidence-next').addEventListener('click', () => moveEvidence(1));
    document.addEventListener('keydown', event => {
      if (byId('evidence').hidden || event.defaultPrevented || event.altKey || event.ctrlKey || event.metaKey || event.shiftKey
          || event.target.closest('input, textarea, select, [contenteditable]:not([contenteditable="false"])')) return;
      const direction = {ArrowUp: -1, ArrowLeft: -1, ArrowDown: 1, ArrowRight: 1}[event.key];
      if (direction) { event.preventDefault(); moveEvidence(direction); }
    });

    function orderKey(probe) { return status.id + ':' + task.name + ':' + probe.id; }
    function orderItem(probe) {
      const [state, label] = stateOf(probe);
      const near = closest(probe);
      return {key: orderKey(probe), job_id: status.id, result_url: status.result_url, task: task.name,
              gene_id: task.gene_id, label: task.name + '-' + probe.id, probe_id: probe.id,
              bridge_id: task.bridge_id, bridge_sequence: task.bridge_sequence, state, state_label: label, quality_retained: !!probe.keep,
              closest_status: near ? null : (assessed(probe) ? 'not_detected' : 'not_assessed'),
              left: {sequence: probe.left_probe.sequence, start: probe.left_probe.position.start, end: probe.left_probe.position.end},
              right: {sequence: probe.right_probe.sequence, start: probe.right_probe.position.start, end: probe.right_probe.position.end},
              synthesis: probe.synthesis || null,
              closest: near ? near.subject_id + ' L' + near.left_mismatches + '/R' + near.right_mismatches : (assessed(probe) ? HCR.t("未检出") : HCR.t("未评估"))};
    }

    function render() {
      byId('track-caption').textContent = task.gene_id + ' · ' + task.sequence_length.toLocaleString(HCR.locale) + HCR.t(" nt · 桥接 ") + task.bridge_id + HCR.t(" · 点击候选或密度柱查看");
      byId('closest-caption').textContent = (polyploid ? HCR.t("多倍体模式：优先展示跨物种，其次跨亚基因组的非目标配对；同类按错配数排序。") : HCR.t("错配条按探针 5′ → 3′ 显示与最接近的非目标配对；")) + HCR.t("灰色为一致碱基，橙色为错配碱基");
      renderTrack(); renderTable();
    }

    async function loadAlignment() {
      const generation = ++alignmentGeneration;
      alignmentRequest?.abort();
      byId('homolog-panel').hidden = !polyploid;
      byId('homolog-detail').hidden = true;
      byId('homolog-retry').hidden = true;
      if (!polyploid) return;
      const index = result.tasks.indexOf(task);
      byId('homolog-summary').textContent = HCR.t("正在从其他已标注分组中查找最相似序列……");
      try {
        let data = alignmentCache.get(index);
        if (!data) {
          alignmentRequest = new AbortController();
          const url = new URL(status.alignment_url, location.href); url.searchParams.set('task', index);
          const response = await fetch(url, {headers: {Accept: 'application/json'}, signal: alignmentRequest.signal});
          data = await response.json();
          if (!response.ok) throw new Error(data.error || HCR.t("同源序列比对未完成。"));
          alignmentCache.set(index, data);
        }
        if (generation !== alignmentGeneration) return;
        alignment = data.status === 'completed' ? data : null;
        alignmentColumns = alignment ? HCRResults.alignmentColumns(task.sequence, alignment) : [];
        byId('homolog-summary').textContent = alignment
          ? HCR.t("{v0}（{v1}） ↔ {v2}（{v3}） · {v4} · 相似度 {v5}% · 目标覆盖 {v6}% · 方向 {v7}", {v0: task.gene_id, v1: data.target_group, v2: data.subject_id, v3: data.subject_group, v4: data.relation === 'cross_species' ? HCR.t('跨物种') : HCR.t('跨亚基因组'), v5: data.identity.toFixed(1), v6: data.query_coverage.toFixed(1), v7: data.strand})
          : data.message;
        byId('homolog-detail').hidden = !alignment;
        if (selected) windowStart = Math.max(1, span(selected)[0] - 12);
        renderTrack(); renderHomolog();
      } catch (error) {
        if (generation !== alignmentGeneration || error.name === 'AbortError') return;
        alignment = null; alignmentColumns = [];
        byId('homolog-summary').textContent = error.message;
        byId('homolog-retry').hidden = false;
        renderTrack();
      }
    }

    function moveAlignmentWindow(delta) {
      windowStart = Math.min(Math.max(1, task.sequence_length - 79), Math.max(1, windowStart + delta));
      renderHomolog(); renderTrack();
    }

    function renderHomolog() {
      if (!polyploid || !alignment) return;
      windowStart = Math.min(Math.max(1, task.sequence_length - 79), Math.max(1, windowStart));
      const end = Math.min(task.sequence_length, windowStart + 79);
      const columns = alignmentColumns.filter(column => column.anchor >= windowStart && column.anchor <= end);
      const subjectPositions = columns.map(column => column.subjectPos).filter(pos => pos !== null);
      byId('homolog-window-label').textContent = HCR.t("目标 {v0}–{v1} / {v2} nt", {v0: windowStart, v1: end, v2: task.sequence_length}) + (subjectPositions.length ? HCR.t(" · 比对序列 {v0}–{v1}", {v0: subjectPositions[0], v1: subjectPositions[subjectPositions.length - 1]}) : HCR.t(" · 此区域未覆盖"));
      byId('homolog-back').disabled = windowStart === 1;
      byId('homolog-forward').disabled = end === task.sequence_length;
      const view = byId('homolog-sequence'); view.replaceChildren();
      const grid = el('div', undefined, 'seq-grid');
      const ruler = el('div', undefined, 'seq-row ruler'), target = el('div', undefined, 'seq-row'), subject = el('div', undefined, 'seq-row');
      columns.forEach(column => {
        const side = !selected || !column.pos ? '' : column.pos >= selected.left_probe.position.start && column.pos <= selected.left_probe.position.end ? 'l'
          : column.pos >= selected.right_probe.position.start && column.pos <= selected.right_probe.position.end ? 'r' : '';
        ruler.append(el('span', column.pos && column.pos % 10 === 0 ? '|' : ''));
        target.append(el('span', column.query, side || 'f'));
        subject.append(el('span', column.subject, column.state === 'mismatch' ? 'x' : column.state === 'gap' ? 'homolog-gap' : column.state === 'match' ? 'same' : 'f'));
      });
      grid.append(el('span', HCR.t("目标坐标 →"), 'name'), ruler, el('span', task.gene_id, 'name'), target,
        el('span', alignment.subject_id + ' · ' + alignment.subject_group, 'name'), subject);
      view.append(grid);
    }

    /* coordinate track: kept pairs to scale, candidate density by status below */
    function renderTrack() {
      const svg = byId('track'); svg.replaceChildren();
      const width = 1000, left = 70, right = 16, inner = width - left - right, length = task.sequence_length;
      const x = position => left + (position - 1) / Math.max(1, length - 1) * inner;
      const kept = task.probe_sets.filter(q => q.keep);
      const bins = Math.min(120, Math.max(20, Math.round(length / 20)));
      const binSize = length / bins;
      const counts = Array.from({length: bins}, () => ({keep: 0, excluded: 0, other: 0, unknown: 0}));
      task.probe_sets.forEach(q => {
        const [a, b] = span(q); const mid = (a + b) / 2;
        counts[Math.min(bins - 1, Math.floor((mid - 1) / binSize))][stateOf(q)[0]]++;
      });
      const maxCount = Math.max(1, ...counts.map(c => c.keep + c.excluded + c.other + c.unknown));
      const keptY = 28, densityTop = 64, densityHeight = 60, axisY = densityTop + densityHeight + 8;
      const hasAlignment = polyploid && alignment, height = axisY + (hasAlignment ? 62 : 26);
      svg.setAttribute('viewBox', '0 0 ' + width + ' ' + height);
      if (range) svg.append(svgEl('rect', {class: 'range', x: x(range[0]), y: 4, width: Math.max(2, x(range[1]) - x(range[0])), height: axisY - 4}));
      if (selected) {
        const [a, b] = span(selected);
        svg.append(svgEl('rect', {class: 'focus', x: x(a) - 2, y: 4, width: Math.max(4, x(b) - x(a) + 4), height: axisY - 4}));
      }
      const label1 = svgEl('text', {class: 'lane-label', x: 0, y: keptY + 10}); label1.textContent = HCR.t("保留 ") + kept.length; svg.append(label1);
      const label2 = svgEl('text', {class: 'lane-label', x: 0, y: densityTop + 14}); label2.textContent = HCR.t("候选密度"); svg.append(label2);
      const label3 = svgEl('text', {class: 'lane-label', x: 0, y: densityTop + 30}); label3.textContent = HCR.t("共 ") + task.probe_sets.length; svg.append(label3);
      svg.append(svgEl('line', {class: 'axis', x1: left, x2: left + inner, y1: keptY + 7, y2: keptY + 7, 'stroke-dasharray': '2 3'}));
      kept.forEach(q => {
        const g = svgEl('g', {class: 'pair', tabindex: 0, role: 'button'});
        const title = svgEl('title'); title.textContent = task.name + '-' + q.id + ' · ' + span(q).join('–'); g.append(title);
        g.append(svgEl('rect', {class: 'arm-l', x: x(q.left_probe.position.start), y: keptY, rx: 1.5, height: 14, width: Math.max(3, x(q.left_probe.position.end) - x(q.left_probe.position.start))}));
        g.append(svgEl('rect', {class: 'arm-r', x: x(q.right_probe.position.start), y: keptY, rx: 1.5, height: 14, width: Math.max(3, x(q.right_probe.position.end) - x(q.right_probe.position.start))}));
        const text = svgEl('text', {class: 'kept-label', x: x(span(q)[0]), y: keptY - 5}); text.textContent = q.id; g.append(text);
        const open = () => select(q);
        g.addEventListener('click', open);
        g.addEventListener('keydown', event => { if (event.key === 'Enter') open(); });
        svg.append(g);
      });
      const barWidth = inner / bins;
      counts.forEach((c, i) => {
        const total = c.keep + c.excluded + c.other + c.unknown;
        if (!total) return;
        const g = svgEl('g', {class: 'bin'});
        const from = Math.floor(i * binSize) + 1, to = Math.min(length, Math.floor((i + 1) * binSize));
        const title = svgEl('title'); title.textContent = from + '–' + to + HCR.t(" nt：保留 ") + c.keep + HCR.t(" · 排除 ") + c.excluded + HCR.t(" · 其他 ") + c.other + HCR.t(" · 未评估 ") + c.unknown; g.append(title);
        let y = densityTop + densityHeight;
        ['keep', 'excluded', 'other', 'unknown'].forEach(key => {
          if (!c[key]) return;
          const h = c[key] / maxCount * densityHeight;
          y -= h;
          g.append(svgEl('rect', {class: 'st-' + key, x: left + i * barWidth + .5, y, width: Math.max(1, barWidth - 1), height: h}));
        });
        g.addEventListener('click', () => { range = [from, to]; byId('result-filter').value = 'all'; renderTrack(); renderTable(); });
        svg.append(g);
      });
      svg.append(svgEl('line', {class: 'axis', x1: left, x2: left + inner, y1: axisY, y2: axisY}));
      const step = [50, 100, 200, 250, 500, 1000, 2000, 5000].find(s => length / s <= 10) || 10000;
      for (let t = 0; t <= length - step * 0.4; t += step) {  // keep clear of the end label
        const pos = Math.max(1, t);
        svg.append(svgEl('line', {class: 'axis', x1: x(pos), x2: x(pos), y1: axisY, y2: axisY + 4}));
        const label = svgEl('text', {x: x(pos), y: axisY + 16, 'text-anchor': t === 0 ? 'start' : 'middle'}); label.textContent = pos; svg.append(label);
      }
      const end = svgEl('text', {x: left + inner, y: axisY + 16, 'text-anchor': 'end'}); end.textContent = length + ' nt'; svg.append(end);
      if (hasAlignment) {
        const y = axisY + 32, lane = svgEl('g', {class: 'homolog-lane'});
        const label = svgEl('text', {class: 'lane-label', x: 0, y: y + 11}); label.textContent = HCR.t("跨组比对"); svg.append(label);
        lane.append(svgEl('rect', {class: 'homolog-background', x: left, y, width: inner, height: 14}));
        alignmentColumns.forEach(column => {
          if (column.state === 'unaligned') return;
          const position = x(column.anchor);
          if (!column.pos) lane.append(svgEl('line', {class: 'homolog-insertion', x1: position, x2: position, y1: y - 2, y2: y + 16}));
          else lane.append(svgEl('rect', {class: 'homolog-' + column.state, x: position, y, width: Math.max(.5, inner / length), height: 14}));
        });
        lane.append(svgEl('rect', {class: 'homolog-window', x: x(windowStart), y: y - 2, width: Math.max(2, x(Math.min(length, windowStart + 79)) - x(windowStart)), height: 18}));
        lane.addEventListener('click', event => {
          const point = svg.createSVGPoint(); point.x = event.clientX; point.y = event.clientY;
          const position = point.matrixTransform(svg.getScreenCTM().inverse()).x;
          windowStart = Math.round(1 + (position - left) / inner * (length - 1)) - 40;
          moveAlignmentWindow(0);
        });
        svg.append(lane);
      }
    }

    function strip(probe) {
      const cell = el('td');
      if (!assessed(probe)) { cell.append(el('span', HCR.t("未评估"), 'hint')); return cell; }
      const near = closest(probe);
      if (!near) { cell.append(el('span', HCR.t("当前检索未检出"), 'hint')); return cell; }
      const bar = el('span', undefined, 'strip');
      bar.setAttribute('aria-hidden', 'true');
      const paint = (length, positions, side) => {
        for (let i = 1; i <= length; i++) bar.append(el('i', undefined, positions.includes(i) ? 'x' : side));
      };
      paint(probe.left_probe.sequence.length, near.left_mismatch_positions_probe_5to3 || [], 'l');
      bar.append(el('span', undefined, 'gap'));
      paint(probe.right_probe.sequence.length, near.right_mismatch_positions_probe_5to3 || [], 'r');
      cell.append(bar, el('span', near.subject_id + ' · ' + (RELATIONS[near.relation] || near.relation) + ' · L ' + near.left_mismatches + ' / R ' + near.right_mismatches, 'strip-label'));
      return cell;
    }

    function visible() {
      const filter = byId('result-filter').value;
      let rows = task.probe_sets.filter(q => filter === 'all' || (filter === 'keep' && q.keep) ||
        (filter === 'cross' && crossExact(q)) || (filter === 'excluded' && policyExcluded(q)) || (filter === 'unassessed' && !assessed(q)));
      if (range) rows = rows.filter(q => { const [a, b] = span(q); const mid = (a + b) / 2; return mid >= range[0] && mid <= range[1]; });
      if (byId('result-sort').value === 'risk') {
        rows = rows.slice().sort((a, b) => {
          const ca = closest(a), cb = closest(b);
          if (!ca || !cb) return (ca ? 0 : 1) - (cb ? 0 : 1);
          return byRisk(ca, cb);
        });
      } else rows = rows.slice().sort((a, b) => span(a)[0] - span(b)[0]);
      return rows;
    }

    function renderTable() {
      const rows = visible(), body = byId('candidate-rows');
      body.replaceChildren();
      const count = byId('visible-count'); count.replaceChildren(HCR.t("显示 ") + rows.length + ' / ' + task.probe_sets.length + HCR.t(" 对"));
      if (range) {
        const chip = el('span', undefined, 'chip info range-chip');
        chip.append(HCR.t("位置 ") + range[0] + '–' + range[1]);
        const clear = el('button', '×'); clear.type = 'button'; clear.setAttribute('aria-label', HCR.t("清除位置范围"));
        clear.addEventListener('click', () => { range = null; renderTrack(); renderTable(); });
        chip.append(clear); count.prepend(chip, ' ');
      }
      const empty = byId('empty-results');
      empty.hidden = rows.length > 0;
      empty.textContent = task.probe_sets.length ? HCR.t("这个范围内没有候选。可切换显示范围或清除位置筛选。") : HCR.t("当前参数未生成候选。可调整臂长、Tm 或 GC 范围后重新设计。");
      rows.forEach(q => {
        const row = el('tr'); if (selected === q) row.className = 'is-selected';
        const id = el('td'); id.append(el('code', task.name + '-' + q.id));
        const [state, label] = stateOf(q);
        const st = el('td'); st.append(el('span', label, 'status-pill ' + (state === 'keep' ? '' : state)));
        const [a, b] = span(q);
        const actions = el('td'); const box = el('div', undefined, 'row-actions');
        const view = el('button', HCR.t("证据"), 'btn small'); view.type = 'button';
        view.addEventListener('click', () => select(q));
        const inOrder = HCR.orders.has(orderKey(q));
        const add = el('button', inOrder ? HCR.t("已加入") : HCR.t("加入清单"), 'btn small' + (inOrder ? ' in-order' : '')); add.type = 'button';
        add.addEventListener('click', () => HCR.orders.toggle(orderItem(q)));
        box.append(view, add); actions.append(box);
        row.append(id, el('td', a + '–' + b, 'num'), st, strip(q), el('td', REASONS[q.reason_code] || q.reason || '—'), actions);
        body.append(row);
      });
      updateEvidenceNavigation(rows);
    }

    function select(probe) {
      selected = probe;
      if (alignment) { windowStart = Math.max(1, span(probe)[0] - 12); renderHomolog(); }
      renderTrack(); renderTable(); renderEvidence(probe);
    }

    function updateEvidenceNavigation(rows = visible()) {
      const index = rows.indexOf(selected);
      byId('evidence-previous').disabled = index <= 0;
      byId('evidence-next').disabled = index < 0 || index >= rows.length - 1;
      byId('evidence-position').textContent = index < 0 ? HCR.t("所选探针对不在当前筛选范围") : HCR.t("第 {v0} / {v1} 对", {v0: index + 1, v1: rows.length});
    }

    function moveEvidence(direction) {
      const next = HCRResults.neighbor(visible(), selected, direction);
      if (next) select(next);
    }

    function renderEvidence(probe, focus = true) {
      const panel = byId('evidence');
      byId('evidence-title').textContent = task.name + '-' + probe.id + HCR.t(" · 配对证据");
      const spec = probe.specificity || {};
      const matches = [...(spec.paired_matches || [])].sort((a, b) => (Number(b.relation === 'target') - Number(a.relation === 'target')) || compareMatches(a, b));
      byId('evidence-summary').textContent = assessed(probe)
        ? (spec.on_target_exact ? HCR.t("已找到目标转录本上的双臂完全配对") : HCR.t("未找到目标完全配对")) + HCR.t(" · 共 ") + matches.length + HCR.t(" 个配对位点 · ") + (REASONS[probe.reason_code] || '')
        : HCR.t("本对探针未完成配对特异性评估。");
      const orderButton = byId('evidence-order');
      const refresh = () => { const has = HCR.orders.has(orderKey(probe)); orderButton.textContent = has ? HCR.t("已加入清单") : HCR.t("加入下单清单"); orderButton.classList.toggle('in-order', has); };
      orderButton.onclick = () => { HCR.orders.toggle(orderItem(probe)); refresh(); };
      refresh();

      const cards = byId('arm-cards'); cards.replaceChildren();
      [['L', 'left_probe', 'l'], ['R', 'right_probe', 'r']].forEach(([name, key, cls]) => {
        const arm = probe[key], card = el('div', undefined, 'arm-card ' + cls);
        card.append(el('strong', name + HCR.t(" 靶向臂 · 5′ → 3′")), el('code', arm.sequence, 'seq'),
          el('span', HCR.t("目标坐标 ") + arm.position.start + '–' + arm.position.end + ' · ' + arm.sequence.length + ' nt · Tm ' + arm.tm + ' °C · GC ' + arm.gc_content + '%', 'hint'));
        const synthesis = probe.synthesis && probe.synthesis[name];
        if (synthesis) {
          const box = el('div', undefined, 'synth');
          const copy = el('button', HCR.t("复制"), 'text-button'); copy.type = 'button'; copy.dataset.copy = synthesis.sequence;
          box.append(HCR.t("合成序列 "), el('b', synthesis.id), ' ', copy, el('code', synthesis.sequence));
          card.append(box);
        }
        cards.append(card);
      });

      renderSequence(probe, matches);

      const body = byId('evidence-matches'); body.replaceChildren();
      matches.forEach(hit => {
        const row = el('tr');
        [hit.subject_id, (RELATIONS[hit.relation] || hit.relation) + (hit.subject_subgenome ? ' · ' + hit.subject_subgenome : ''),
         hit.left_start + '–' + hit.left_end + ' / ' + hit.left_mismatches, hit.right_start + '–' + hit.right_end + ' / ' + hit.right_mismatches,
         hit.strand + ' / ' + hit.gap + ' nt'].forEach(v => row.append(el('td', v)));
        body.append(row);
      });
      panel.hidden = false;
      if (focus) {
        panel.focus({preventScroll: true});
        panel.scrollIntoView({block: 'start'});
      }
    }

    /* sense-strand view: target footprint with flanks, off-target footprints aligned underneath */
    function renderSequence(probe, matches) {
      const view = byId('seqview'); view.replaceChildren();
      const sequence = task.sequence;
      const L = probe.left_probe.position, R = probe.right_probe.position;
      if (!sequence) { view.append(el('p', HCR.t("该结果未保存目标序列，无法显示序列视图。"), 'hint')); return; }
      const from = Math.max(1, Math.min(L.start, R.start) - 12), to = Math.min(sequence.length, Math.max(L.end, R.end) + 12);
      const arm = pos => (pos >= L.start && pos <= L.end ? 'l' : pos >= R.start && pos <= R.end ? 'r' : '');
      const grid = el('div', undefined, 'seq-grid');
      const ruler = el('div', undefined, 'seq-row ruler');
      for (let pos = from; pos <= to; pos++) ruler.append(el('span', pos % 10 === 0 ? '|' : pos % 5 === 0 ? '·' : ''));
      grid.append(el('span', from + ' →', 'name'), ruler);
      const target = el('div', undefined, 'seq-row');
      for (let pos = from; pos <= to; pos++) target.append(el('span', sequence[pos - 1], arm(pos) || 'f'));
      const targetName = el('span', undefined, 'name'); targetName.append(HCR.t("目标 "), el('b', task.gene_id));
      grid.append(targetName, target);
      matches.filter(m => m.relation !== 'target').slice(0, 10).forEach(m => {
        const row = el('div', undefined, 'seq-row');
        for (let pos = from; pos <= to; pos++) {
          const side = arm(pos);
          if (!side) { row.append(el('span', '·', 'f')); continue; }
          const segment = side === 'l' ? m.left_matched_target : m.right_matched_target;
          const start = side === 'l' ? L.start : R.start;
          const base = segment ? segment[pos - start] : '?';
          row.append(el('span', base, base === sequence[pos - 1] ? side + ' same' : 'x'));
        }
        const name = el('span', undefined, 'name');
        name.append(el('b', m.subject_id), ' ' + (RELATIONS[m.relation] || m.relation) + ' · L' + m.left_mismatches + '/R' + m.right_mismatches);
        grid.append(name, row);
      });
      const others = matches.filter(m => m.relation !== 'target').length;
      view.append(grid);
      if (others > 10) view.append(el('p', HCR.t("另有 ") + (others - 10) + HCR.t(" 个非目标配对见下表。"), 'hint'));
      if (!others && assessed(probe)) view.append(el('p', HCR.t("当前检索范围内没有非目标配对。"), 'hint'));
    }

    render(); loadAlignment();
  }

  if (status.status !== 'completed') { followProgress(); return; }
  fetch(status.result_url, {headers: {Accept: 'application/json'}})
    .then(response => response.ok ? response.json() : Promise.reject(new Error('HTTP ' + response.status)))
    .then(showResults)
    .catch(error => {
      const box = el('div', HCR.t("结果读取失败：") + error.message + HCR.t("。可以下载完整 ZIP 查看。"), 'alert error');
      byId('result-stats').before(box);
    });
})();
