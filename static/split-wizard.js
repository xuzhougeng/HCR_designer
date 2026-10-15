(() => {
  'use strict';
  const config = JSON.parse(document.getElementById('split-config').textContent);
  const root = document.getElementById('split-workspace');
  const form = document.getElementById('split-form');
  const byId = id => document.getElementById(id);
  const field = name => form.elements.namedItem(name);
  let mode = config.mode, step = 1, prepared = null, busy = false;
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
    for (let n = 1; n <= 2; n++) byId('step-' + n).hidden = n !== step;
    root.querySelectorAll('[data-step]').forEach(button => {
      const n = Number(button.dataset.step);
      if (n === step) button.setAttribute('aria-current', 'step'); else button.removeAttribute('aria-current');
      button.disabled = busy || (n === 2 && !prepared) || n === 3;
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
    const selectedReference = field('existing_ref_genome').selectedOptions[0];
    const polyploid = refMode === 'existing' && selectedReference?.dataset.ploidy === 'true';
    const hasAnnotation = polyploid && selectedReference?.dataset.subgenomeReady === 'true';
    byId('multiple-sets-field').hidden = !polyploid;
    if (!hasAnnotation) field('process_multiple_sets').checked = false;
    byId('multiple-sets-hint').textContent = hasAnnotation
      ? '使用参考库已保存的亚基因组归属；结果页提供跨组同源比对与跨物种优先的配对证据，无需再次上传。'
      : '此参考库尚无有效的亚基因组归属文件，请先到参考库管理页上传，刷新本页后即可勾选。';
    byId('subgenome-note').hidden = refMode === 'none';
    field('existing_ref_genome').disabled = busy || refMode !== 'existing';
    field('ref_genome').disabled = busy || refMode !== 'upload';
    field('subgenome_files').disabled = busy || refMode !== 'subgenome_files';
    field('process_multiple_sets').disabled = busy || !hasAnnotation;
    const useAnnotation = hasAnnotation && field('process_multiple_sets').checked;
    const subgenomeOption = field('target_mode').querySelector('option[value="subgenome"]');
    subgenomeOption.disabled = !useAnnotation && refMode !== 'subgenome_files';
    if (useAnnotation) field('target_mode').value = 'subgenome';
    else if (subgenomeOption.disabled && field('target_mode').value === 'subgenome') field('target_mode').value = 'species';
    field('target_mode').disabled = busy || useAnnotation;
    const bridge = field('bp_id').value.trim().toUpperCase();
    byId('bridge-sequence').textContent = config.bridges[bridge] || '编号不在当前桥接探针库中';
    byId('reference-note').textContent = refMode === 'none'
      ? '本次只做质量筛选，结果将标记“特异性未评估”。用于物种来源判定时请选择参考库。'
      : '参考库用于查找目标及非目标转录本上的双臂配对。Bo/Bju cDNA 是转录本参考库。';
    if (refMode === 'existing' && field('existing_ref_genome').selectedOptions[0]?.dataset.ploidy === 'true') {
      byId('reference-note').textContent += ' 此参考已配置为多倍体，可为本次任务选择是否考虑多套染色体。';
    }
  }
  function invalidate() {
    prepared = null;
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
    if (reference.polyploid !== undefined && reference.polyploid !== null) text += '\n参考配置：' + (reference.polyploid ? '多倍体' : '非多倍体');
    const sub = reference.subgenome;
    if (sub) {
      const counts = Object.entries(sub.subgenome_counts).map(([key, value]) => key + ' ' + value.toLocaleString()).join(' / ');
      text += '\n亚基因组（' + sub.source + '）：' + counts + ' · 未标注 ' + (sub.unannotated || 0).toLocaleString() + ' 条';
    }
    return text;
  }
  function renderPrepared() {
    byId('reference-summary').textContent = referenceText(prepared.reference) + (prepared.reference?.polyploid ? '\n使用参考库的亚基因组归属：' + (prepared.process_multiple_sets ? '是' : '否') : '');
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
    if (refMode === 'none') data.delete('process_multiple_sets');
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
    setBusy(true, '正在提交到后台队列……');
    try {
      const job = await requestJSON(prepared.run_url, {method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(values)});
      window.location.href = job.result_url;
      return;
    } catch (error) {showMessage(error.message, 'error');}
    setBusy(false);
    byId('workflow-message').scrollIntoView({block:'nearest'});
  }
  form.addEventListener('submit',e => e.preventDefault());
  byId('step-1').addEventListener('input',() => {if (!busy) invalidate();});
  byId('step-1').addEventListener('change',() => {if (!busy) invalidate();});
  field('max_selected').addEventListener('input',() => {if (prepared) renderPrepared();});
  function defaultMultipleSets() {
    field('process_multiple_sets').checked = false;
  }
  field('existing_ref_genome').addEventListener('change', defaultMultipleSets);
  root.querySelectorAll('[data-mode]').forEach(b => b.addEventListener('click',() => {mode = b.dataset.mode;invalidate();}));
  root.querySelectorAll('[data-step]').forEach(b => b.addEventListener('click',() => {if (!busy) {showMessage('');setStep(Number(b.dataset.step));}}));
  root.querySelectorAll('[data-back]').forEach(b => b.addEventListener('click',() => {showMessage('');setStep(1);}));
  byId('fill-example').addEventListener('click',() => {
    field('name').value = 'ABCG13_demo';field('geneID').value = 'BjuVA06G01970';field('bp_id').value = 'BP_0079';field('sequenceType').value = 'reference';
    field('batch_file').value = '';field('batch_text').value = 'task_name,BP_ID,gene_id,max_selected,allowed_ids\nABCG13,BP_0079,BjuVA06G01970,5,';field('allowed_ids').value = '';field('subgenome_files').value = '';
    if (config.default_reference) {field('reference_mode').value = 'existing';field('existing_ref_genome').value = config.default_reference;}
    defaultMultipleSets();
    invalidate();
  });
  byId('prepare-button').addEventListener('click',prepare);
  byId('run-button').addEventListener('click',run);
  field('existing_ref_genome').value = config.default_reference;
  defaultMultipleSets();
  if (!config.default_reference && field('existing_ref_genome').options.length === 1) field('reference_mode').value = 'upload';
  syncInput();setStep(1);
})();
