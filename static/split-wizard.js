(() => {
  'use strict';
  const config = JSON.parse(document.getElementById('split-config').textContent);
  const root = document.getElementById('split-workspace');
  const form = document.getElementById('split-form');
  const byId = id => document.getElementById(id);
  const field = name => form.elements.namedItem(name);
  let mode = config.mode, step = 1, prepared = null, result = config.result, busy = false;
  const reasons = {
    cross_species_exact_pair: '另一物种存在双臂完全配对',
    cross_species_near_pair: '另一物种存在阈值内的近似配对',
    unclassified_exact_pair: '未知物种的非目标转录本存在完全配对',
    unclassified_near_pair: '未知物种的非目标转录本存在阈值内的近似配对',
    same_species_exact_pair: '同物种非允许转录本存在完全配对',
    same_species_near_pair: '同物种非允许转录本存在阈值内的近似配对',
    specificity_not_evaluated: '配对特异性证据缺失或不完整',
    no_exact_on_target_pair: '未找到目标转录本上的双臂完全配对',
    gc_junction: 'L 臂连接位点未满足当前筛选规则',
    too_close: '与其他已选组合间距不足',
    selection_limit_or_rank: '排序后未进入本次保留数量',
    selected: '按当前规则保留，仍需进一步复核'
  };
  const relations = {target:'目标转录本',allowed:'允许共同检测',same_species_other:'同物种其他转录本',cross_species:'跨物种',unknown_species:'物种未知'};
  const modes = {species:'物种模式',transcript:'严格模式',allowlist:'允许集合',subgenome:'亚基因组模式'};
  function element(tag, text, className) {
    const node = document.createElement(tag);
    if (text !== undefined) node.textContent = text;
    if (className) node.className = className;
    return node;
  }
  function showMessage(text, kind = '') {
    const box = byId('workflow-message');
    box.hidden = !text; box.textContent = text;
    box.className = 'wizard-alert ' + kind;
    box.setAttribute('role', kind === 'error' ? 'alert' : 'status');
  }
  function setStep(next) {
    step = next;
    for (let n = 1; n <= 3; n++) byId('step-' + n).hidden = n !== step;
    root.querySelectorAll('[data-step]').forEach(button => {
      const n = Number(button.dataset.step);
      if (n === step) button.setAttribute('aria-current', 'step'); else button.removeAttribute('aria-current');
      button.disabled = busy || (n === 2 && !prepared) || (n === 3 && !result);
    });
  }
  function syncInput() {
    field('input_mode').value = mode;
    root.querySelectorAll('[data-mode]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.mode === mode)));
    byId('single-fields').hidden = mode !== 'single';
    byId('batch-fields').hidden = mode !== 'batch';
    byId('paste-field').hidden = field('sequenceType').value !== 'custom';
    const refMode = field('reference_mode').value;
    byId('existing-reference-field').hidden = refMode !== 'existing';
    byId('upload-reference-field').hidden = refMode !== 'upload';
    byId('subgenome-files-field').hidden = refMode !== 'subgenome_files';
    byId('subgenome-map-field').hidden = refMode !== 'existing' && refMode !== 'upload';
    byId('subgenome-note').hidden = refMode === 'none';
    field('existing_ref_genome').disabled = busy || refMode !== 'existing';
    field('ref_genome').disabled = busy || refMode !== 'upload';
    field('subgenome_files').disabled = busy || refMode !== 'subgenome_files';
    field('subgenome_map').disabled = busy || (refMode !== 'existing' && refMode !== 'upload');
    const bridge = field('bp_id').value.trim().toUpperCase();
    byId('bridge-sequence').textContent = config.bridges[bridge] || '编号不在当前桥接探针库中';
    byId('reference-note').textContent = refMode === 'none'
      ? '本次只做质量筛选，结果将标记“特异性未评估”。用于物种来源判定时请选择参考库。'
      : '参考库用于查找目标及非目标转录本上的双臂配对。Bo/Bju cDNA 是转录本参考库。';
  }
  function invalidate() {
    if (result) history.replaceState(null, '', mode === 'batch' ? config.batch_url : config.single_url);
    prepared = null;
    result = null;
    setStep(1); syncInput(); showMessage('');
  }
  function setBusy(value, message = '') {
    busy = value;
    form.querySelectorAll('input,select,textarea,button').forEach(e => e.disabled = value);
    root.querySelectorAll('[data-mode]').forEach(e => e.disabled = value);
    setStep(step); syncInput();
    if (value) showMessage(message, 'busy');
    root.setAttribute('aria-busy', String(value));
  }
  async function requestJSON(url, options) {
    const response = await fetch(url, options);
    const data = await response.json().catch(() => ({error:'服务器未返回完整结果，请重试。'}));
    if (!response.ok) throw Object.assign(new Error(data.error || '请求未完成。'), {field:data.field});
    return data;
  }
  function referenceText(reference) {
    if (!reference) return '未指定参考库 · 特异性未评估';
    const species = Object.entries(reference.species_counts).map(([key, value]) => key + ' ' + value.toLocaleString()).join(' / ');
    let text = reference.name + ' · ' + reference.sequence_count.toLocaleString() + ' 条序列\n' + species;
    const sub = reference.subgenome;
    if (sub) {
      const counts = Object.entries(sub.subgenome_counts).map(([key, value]) => key + ' ' + value.toLocaleString()).join(' / ');
      text += '\n亚基因组（' + sub.source + '）：' + counts + ' · 未标注 ' + (sub.unannotated || 0).toLocaleString() + ' 条';
    }
    return text;
  }
  function renderPrepared() {
    byId('reference-summary').textContent = referenceText(prepared.reference);
    const body = byId('prepared-tasks'); body.replaceChildren();
    prepared.tasks.forEach(task => {
      const row = element('tr');
      [task.name, task.gene_id || '未指定', task.sequence_length + ' nt', task.bridge_id,
        task.max_selected === null ? '默认 ' + field('max_selected').value + ' 对' : task.max_selected + ' 对',
        task.subgenome || '未标注',
        (task.allowed_ids || []).length ? task.allowed_ids.join('; ') : '无（按目标判定模式）'].forEach(value => row.append(element('td', value)));
      body.append(row);
    });
  }
  async function prepare() {
    if (busy) return;
    if (prepared) {renderPrepared();setStep(2);return;}
    const data = new FormData(form);
    if (mode === 'single') {data.delete('batch_file');data.delete('batch_text');}
    else ['name','geneID','bp_id','seq','sequenceType'].forEach(k => data.delete(k));
    const refMode = field('reference_mode').value;
    if (refMode !== 'existing') data.delete('existing_ref_genome');
    if (refMode !== 'upload') data.delete('ref_genome');
    if (refMode !== 'subgenome_files') data.delete('subgenome_files');
    if (refMode !== 'existing' && refMode !== 'upload') data.delete('subgenome_map');
    setBusy(true, '正在核对目标、桥接探针和参考序列……大文件首次检查可能需要一些时间。');
    let error;
    try {
      prepared = await requestJSON(config.prepare_url, {method:'POST',body:data});
      renderPrepared(); showMessage(''); setStep(2);
    } catch (caught) {error = caught;showMessage(error.message, 'error');}
    finally {setBusy(false);}
    if (error) byId('workflow-message').scrollIntoView({block:'nearest'});
    else byId('parameter-heading').scrollIntoView({block:'nearest'});
  }
  function readParameters() {
    const values = {};
    for (const input of form.querySelectorAll('[data-parameter]')) {
      if (!input.value.trim() || !input.checkValidity()) {
        input.closest('details')?.setAttribute('open', '');
        input.focus(); input.reportValidity();
        throw new Error('请检查参数范围，整数参数不能填写小数。');
      }
      values[input.name] = input.value;
    }
    for (const [low, high, label] of [['min_length','max_length','臂长'],['min_tm','max_tm','Tm'],['min_gc','max_gc','GC']]) {
      if (Number(values[low]) > Number(values[high])) throw new Error(label + ' 的最小值不能大于最大值。');
    }
    return values;
  }
  async function run() {
    if (busy || !prepared) return;
    let values;
    try {values = readParameters();} catch (error) {showMessage(error.message, 'error');return;}
    setBusy(true, '正在生成候选、检查双臂配对并整理结果……完成前请保持本页打开，请勿重复提交。');
    try {
      result = await requestJSON(prepared.run_url, {method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(values)});
      renderResult();showMessage('');setStep(3);
      history.replaceState(null, '', result.result_url);
      prepared = null;
    } catch (error) {showMessage(error.message, 'error');}
    finally {setBusy(false);}
    if (result) byId('results-heading').scrollIntoView({block:'nearest'});
    else byId('workflow-message').scrollIntoView({block:'nearest'});
  }
  const assessed = p => p.specificity?.status === 'completed';
  const crossExact = p => (p.specificity?.counts?.cross_species?.[0] || 0) > 0;
  const policyExcluded = p => /_(exact|near)_pair$/.test(p.reason_code || '');
  function status(p) {
    if (!assessed(p)) return [p.keep ? '未评估 · 按质量保留' : '特异性未评估', 'unknown'];
    if (p.keep) return ['保留待复核', ''];
    return [policyExcluded(p) ? '按目标策略排除' : '未保留', 'excluded'];
  }
  function crossText(p) {
    if (!assessed(p)) return '未评估';
    const pairs = p.specificity.paired_matches.filter(h => h.relation === 'cross_species');
    if (!pairs.length) return '当前检索未检出';
    pairs.sort((a,b) => Math.max(a.left_mismatches,a.right_mismatches)-Math.max(b.left_mismatches,b.right_mismatches) || a.left_mismatches+a.right_mismatches-b.left_mismatches-b.right_mismatches);
    return 'L/R：' + pairs[0].left_mismatches + ' / ' + pairs[0].right_mismatches + ' 错配';
  }
  function renderResult() {
    byId('result-subtitle').textContent = result.tasks.length + ' 个任务完成 · ' + new Date(result.completed_at).toLocaleString('zh-CN');
    byId('download-results').href = result.download_url;
    const stats = byId('result-stats');stats.replaceChildren();
    [['原始候选','candidates'],['跨物种完全配对','cross_exact'],['按目标策略排除','excluded'],['保留待复核','selected'],['特异性未评估','unassessed']].forEach(([label,key]) => {
      const card = element('div',undefined,'result-stat');
      card.append(element('span',label),element('strong',result.tasks.reduce((sum,t) => sum+t[key],0)+' 对'));stats.append(card);
    });
    byId('result-limit').textContent = !result.reference
      ? '本次未指定参考库，未进行特异性评估。保留候选仅反映质量筛选，不能用于判断物种特异性。'
      : result.search_exhaustive
        ? '双臂配对已在“正向转录本、无缺口、每臂最多 3 个替换”的范围内完整检出；插入缺失未复核。候选保留和“未检出”均不等于特异性合格。'
        : '本次有任务未完成完整配对扫描，相关候选标记为未评估；未进行插入缺失复核。候选保留和“未检出”均不等于特异性合格。';
    const settings = byId('run-settings-content');settings.replaceChildren(element('div',referenceText(result.reference)));
    if (result.reference) settings.append(element('code','参考 SHA256：'+result.reference.sha256));
    const p = result.parameters;
    settings.append(element('div','臂长 '+p.min_length+'–'+p.max_length+' nt · Tm '+p.min_tm+'–'+p.max_tm+' °C · GC '+p.min_gc+'–'+p.max_gc+'%'));
    settings.append(element('div','默认最多 '+p.max_selected+' 对/基因 · 组合间距 '+p.min_gap+' nt · 臂内间隔 0 nt'));
    settings.append(element('div','连续同碱基阈值 '+p.poly_n+' · 互补长度 '+p.min_complementary_length+' nt · 重复片段 '+p.kmer_size+' nt / '+p.min_kmer_count+' 次 · 连接位点惩罚 '+p.gc_junction_penalty));
    settings.append(element('div','目标判定模式 '+(modes[p.target_mode] || p.target_mode)+' · 非允许配对排除阈值 每臂 ≤ '+p.exclusion_max_mismatches+' 错配'));
    if (p.target_mode === 'subgenome') settings.append(element('div','各任务允许的亚基因组：'+result.tasks.map(t => t.name+' → '+(t.target_policy?.target_subgenome || t.subgenome || '未标注')).join('；')));
    settings.append(element('div',result.tasks.map(t => t.name+'：'+(modes[t.target_policy?.mode] || '')+((t.target_policy?.allowed_ids || []).length ? '（'+t.target_policy.allowed_ids.join('; ')+'）' : '')+'，最多 '+t.max_selected+' 对，实际 '+t.selected+' 对'+(t.search?.seconds !== undefined ? '，扫描 '+t.search.seconds+' s' : '')).join('；')));
    const tasks = byId('result-task');tasks.replaceChildren();
    result.tasks.forEach((task,i) => {const o = element('option',task.name+' · '+task.gene_id);o.value = i;tasks.append(o);});
    renderCandidates();
  }
  function renderCandidates() {
    if (!result) return;
    const task = result.tasks[Number(byId('result-task').value) || 0];
    const filter = byId('result-filter').value;
    const probes = task.probe_sets.filter(p => filter === 'all' || (filter === 'keep' && p.keep) || (filter === 'cross' && crossExact(p)) || (filter === 'excluded' && policyExcluded(p)) || (filter === 'unassessed' && !assessed(p)));
    const body = byId('candidate-rows');body.replaceChildren();byId('pair-evidence').hidden = true;
    byId('visible-count').textContent = '显示 '+probes.length+' / '+task.candidates+' 对';
    byId('empty-results').hidden = probes.length > 0;
    byId('empty-results').textContent = task.candidates ? '这个范围内没有候选。可切换显示范围。' : '当前参数未生成候选。可调整臂长、Tm 或 GC 范围后重新设计。';
    for (const probe of probes) {
      const row = element('tr');
      const id = element('td');id.append(element('code',task.name+'-'+probe.id));
      const state = element('td');const [label,kind] = status(probe);state.append(element('span',label,'status-pill '+kind));
      const action = element('td'),button = element('button','查看配对','text-button');button.type = 'button';
      button.setAttribute('aria-label','查看 '+task.name+'-'+probe.id+' 的配对证据');
      button.addEventListener('click',() => {body.querySelectorAll('tr').forEach(r => r.classList.toggle('selected-row',r===row));renderEvidence(task,probe);});action.append(button);
      row.append(id,state,element('td',crossText(probe)),element('td',reasons[probe.reason_code] || '请查看完整结果中的筛选说明'),action);body.append(row);
    }
  }
  function renderEvidence(task, probe) {
    byId('evidence-heading').textContent = task.name+'-'+probe.id+' · 配对证据';
    const arms = byId('evidence-arms');arms.replaceChildren();
    [['L','left_probe'],['R','right_probe']].forEach(([label,key]) => {
      const arm = probe[key],box = element('div',undefined,'evidence-arm');
      box.append(element('strong',label+' 靶向臂'),element('code',arm.sequence),element('span','输入序列坐标 '+arm.position.start+'–'+arm.position.end+' · Tm '+arm.tm+' °C · GC '+arm.gc_content+'%'));arms.append(box);
    });
    const body = byId('evidence-matches');body.replaceChildren();
    const spec = probe.specificity || {};
    const matches = [...(spec.paired_matches || [])].sort((a,b) => (a.relation === 'cross_species' ? 0 : 1)-(b.relation === 'cross_species' ? 0 : 1));
    byId('evidence-summary').textContent = assessed(probe)
      ? (spec.on_target_exact ? '已找到目标转录本上的双臂完全配对。' : '未找到目标完全配对。')+' 以下按同一转录本、方向及设计间隔组合，共 '+matches.length+' 个配对位点。'
      : '本对探针未完成配对特异性评估，没有可用的配对结论。';
    matches.forEach(hit => {
      const row = element('tr');
      [hit.subject_id,(relations[hit.relation] || hit.relation) + (hit.subject_subgenome ? ' · ' + hit.subject_subgenome : ''),
        hit.left_start+'–'+hit.left_end+' / '+hit.left_mismatches+' 错配',
        hit.right_start+'–'+hit.right_end+' / '+hit.right_mismatches+' 错配',
        hit.strand+' / '+hit.gap+' nt'].forEach(v => row.append(element('td',v)));body.append(row);
    });
    byId('pair-evidence').hidden = false;byId('pair-evidence').scrollIntoView({block:'nearest'});
  }
  form.addEventListener('submit',e => e.preventDefault());
  byId('step-1').addEventListener('input',() => {if (!busy) invalidate();});
  byId('step-1').addEventListener('change',() => {if (!busy) invalidate();});
  field('max_selected').addEventListener('input',() => {if (prepared) renderPrepared();});
  root.querySelectorAll('[data-mode]').forEach(b => b.addEventListener('click',() => {mode = b.dataset.mode;invalidate();}));
  root.querySelectorAll('[data-step]').forEach(b => b.addEventListener('click',() => {if (!busy) {showMessage('');setStep(Number(b.dataset.step));}}));
  root.querySelectorAll('[data-back]').forEach(b => b.addEventListener('click',() => {showMessage('');setStep(1);}));
  byId('fill-example').addEventListener('click',() => {
    field('name').value = 'ABCG13_demo';field('geneID').value = 'BjuVA06G01970';field('bp_id').value = 'BP_0079';field('sequenceType').value = 'reference';
    field('batch_file').value = '';field('batch_text').value = 'task_name,BP_ID,gene_id,max_selected,allowed_ids\nABCG13,BP_0079,BjuVA06G01970,5,';field('allowed_ids').value = '';field('subgenome_map').value = '';field('subgenome_files').value = '';
    if (config.default_reference) {field('reference_mode').value = 'existing';field('existing_ref_genome').value = config.default_reference;}
    invalidate();
  });
  byId('prepare-button').addEventListener('click',prepare);
  byId('run-button').addEventListener('click',run);
  byId('result-task').addEventListener('change',renderCandidates);
  byId('result-filter').addEventListener('change',renderCandidates);
  byId('close-evidence').addEventListener('click',() => byId('pair-evidence').hidden = true);
  field('existing_ref_genome').value = config.default_reference;
  if (!config.default_reference && field('existing_ref_genome').options.length === 1) field('reference_mode').value = 'upload';
  syncInput();setStep(1);
  if (result) {renderResult();setStep(3);}
})();
