(() => {
  'use strict';
  const form = document.getElementById('subgenome-form');
  const fileInput = document.getElementById('group-fastas');
  const ruleInput = document.getElementById('rule-fasta');
  const rulesBox = document.getElementById('subgenome-rules');
  const resultBox = document.getElementById('subgenome-result');
  const status = document.getElementById('subgenome-status');
  const submit = document.getElementById('generate-subgenomes');
  const download = document.getElementById('download-subgenomes');
  const accept = document.getElementById('accept-unmatched');
  let mode = 'files', result = null, version = 0;
  const node = (tag, text, cls) => {
    const element = document.createElement(tag);
    if (text != null) element.textContent = text;
    if (cls) element.className = cls;
    return element;
  };
  function invalidate() {
    version++; result = null; resultBox.hidden = true; download.disabled = true;
    HCR.setAlert(status, '');
  }
  function addRule(pattern = '', label = '') {
    const row = document.getElementById('subgenome-rule-template').content.cloneNode(true);
    row.querySelector('.rule-pattern').value = pattern;
    row.querySelector('.rule-label').value = label;
    row.querySelector('.remove-rule').addEventListener('click', event => { event.target.closest('tr').remove(); invalidate(); });
    rulesBox.append(row);
  }
  addRule('Axx', 'A'); addRule('Bxx', 'B');
  document.getElementById('add-subgenome-rule').addEventListener('click', () => {
    if (rulesBox.children.length >= 50) { HCR.setAlert(status, '最多支持 50 条规则。', 'error'); return; }
    addRule(); invalidate();
  });
  document.querySelectorAll('[data-mode]').forEach(button => button.addEventListener('click', () => {
    mode = button.dataset.mode;
    document.querySelectorAll('[data-mode]').forEach(b => b.setAttribute('aria-pressed', String(b === button)));
    document.getElementById('file-mode').hidden = mode !== 'files';
    document.getElementById('rule-mode').hidden = mode !== 'rules';
    invalidate();
  }));
  fileInput.addEventListener('change', () => {
    const body = document.getElementById('file-labels'); body.replaceChildren();
    [...fileInput.files].forEach(file => {
      const row = node('tr');
      const label = document.createElement('input'); label.className = 'file-subgenome-label';
      label.maxLength = 64; label.setAttribute('aria-label', file.name + ' 的亚基因组名称');
      label.value = file.name.replace(/\.gz$/i, '').replace(/\.(fa|fasta|fna)$/i, '').trim().replace(/\s+/g, '_');
      const cell = node('td'); cell.append(label); row.append(node('td', file.name), cell); body.append(row);
    });
    document.getElementById('file-labels-box').hidden = !fileInput.files.length;
  });
  form.addEventListener('input', invalidate);
  form.addEventListener('change', invalidate);
  function enableDownload() {
    download.disabled = !result?.csv || (result.summary.unmatched > 0 && !accept.checked);
  }
  accept.addEventListener('change', enableDownload);
  function render(data) {
    result = data; accept.checked = false;
    const s = data.summary;
    document.getElementById('subgenome-summary').textContent = `${s.files} 个文件，${s.headers.toLocaleString()} 条 FASTA 标题；重复编号合并 ${s.duplicates.toLocaleString()} 条。`;
    const tiles = [['独立编号', s.total], ['已归属', s.assigned], ['未匹配', s.unmatched], ['归属冲突', s.conflicts]];
    document.getElementById('subgenome-stats').replaceChildren(...tiles.map(([label, value]) => {
      const tile = node('div', null, 'stat'); tile.append(node('span', label), node('strong', value.toLocaleString())); return tile;
    }));
    document.getElementById('subgenome-counts').textContent = Object.entries(s.subgenome_counts).map(([label, count]) => `${label}：${count.toLocaleString()} 条`).join(' · ') || '尚无可导出的归属记录。';
    document.getElementById('subgenome-preview').replaceChildren(...data.preview.map(pair => {
      const row = node('tr'); row.append(node('td', pair.identifier), node('td', pair.subgenome)); return row;
    }));
    document.getElementById('subgenome-issues').replaceChildren(...data.issues.map(issue => {
      const row = node('tr'); row.append(node('td', issue.identifier), node('td', issue.issue === 'conflict' ? '多个亚基因组冲突' : '未匹配'), node('td', issue.labels.join('、') || '—')); return row;
    }));
    document.getElementById('subgenome-issues-box').hidden = !data.issues.length;
    document.getElementById('download-subgenome-issues').hidden = !data.issues_csv;
    document.getElementById('accept-unmatched-box').hidden = !s.unmatched || !data.csv;
    HCR.setAlert(document.getElementById('subgenome-notice'), s.conflicts
      ? `${s.conflicts} 个编号同时归入多个亚基因组，请调整文件归属或匹配规则后重新生成。`
      : s.unmatched ? `${s.unmatched} 个编号没有匹配。可调整规则后重新生成，或确认后仅导出已匹配编号。` : '', s.conflicts ? 'error' : 'notice');
    resultBox.hidden = false; enableDownload();
  }
  form.addEventListener('submit', async event => {
    event.preventDefault();
    const files = [...(mode === 'files' ? fileInput : ruleInput).files];
    if (!files.length) { HCR.setAlert(status, '请先选择 FASTA 文件。', 'error'); return; }
    invalidate();
    const current = version, body = new FormData();
    body.append('mode', mode); files.forEach(file => body.append('fasta_files', file));
    if (mode === 'files') {
      body.append('labels', JSON.stringify([...document.querySelectorAll('.file-subgenome-label')].map(input => input.value)));
    } else {
      body.append('rules', JSON.stringify([...rulesBox.children].map(row => ({pattern: row.querySelector('.rule-pattern').value, label: row.querySelector('.rule-label').value}))));
      body.append('position', document.getElementById('rule-position').value);
      body.append('rule_mode', document.getElementById('rule-syntax').value);
      if (document.getElementById('rule-case').checked) body.append('case_sensitive', 'on');
    }
    submit.disabled = true; HCR.setAlert(status, '正在提取编号并匹配归属……', 'busy');
    try {
      const response = await fetch(form.dataset.url, {method: 'POST', body});
      const data = await response.json();
      if (version !== current) return;
      if (!response.ok) throw new Error(data.error || '生成未完成。');
      render(data); HCR.setAlert(status, '');
    } catch (error) {
      if (version === current) HCR.setAlert(status, error.message, 'error');
    } finally { submit.disabled = false; }
  });
  const saveCSV = (text, filename) => HCR.saveBlob(new Blob(['\ufeff', text], {type: 'text/csv;charset=utf-8'}), filename);
  download.addEventListener('click', () => { if (!download.disabled && result?.csv) saveCSV(result.csv, 'subgenomes.csv'); });
  document.getElementById('download-subgenome-issues').addEventListener('click', () => { if (result?.issues_csv) saveCSV(result.issues_csv, 'subgenome_issues.csv'); });
})();
