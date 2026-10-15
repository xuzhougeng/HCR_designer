(() => {
  'use strict';
  const config = JSON.parse(document.getElementById('admin-config').textContent);
  const formatBytes = value => {
    const units = ['B', 'KB', 'MB', 'GB', 'TB'];
    let size = Number(value), i = 0;
    while (size >= 1024 && i < units.length - 1) { size /= 1024; i++; }
    return size.toLocaleString('zh-CN', {maximumFractionDigits: i ? 2 : 0}) + ' ' + units[i];
  };
  document.querySelectorAll('[data-bytes]').forEach(el => { el.textContent = formatBytes(el.dataset.bytes); });
  document.querySelectorAll('[data-date]').forEach(el => { el.textContent = new Date(el.dataset.date).toLocaleString('zh-CN'); });
  const rows = [...document.querySelectorAll('tr[data-key]')];
  const selectAll = document.getElementById('select-visible');
  const deleteButton = document.getElementById('delete-selected');
  const cleanupButton = document.getElementById('cleanup-results');
  const message = document.getElementById('admin-message');
  const selected = () => rows.filter(row => row.querySelector('input').checked);
  let busy = false;
  function update() {
    const items = selected();
    deleteButton.disabled = busy || !items.length;
    cleanupButton.disabled = busy;
    document.getElementById('selected-summary').textContent = `已选 ${items.length} 项 · ${formatBytes(items.reduce((sum, r) => sum + Number(r.dataset.size), 0))}`;
    const visible = rows.filter(r => !r.hidden && !r.querySelector('input').disabled);
    const count = visible.filter(r => r.querySelector('input').checked).length;
    selectAll.checked = visible.length > 0 && count === visible.length;
    selectAll.indeterminate = count > 0 && count < visible.length;
  }
  function filter() {
    const kind = document.getElementById('storage-kind').value;
    const query = document.getElementById('storage-search').value.trim().toLowerCase();
    rows.forEach(row => {
      row.hidden = Boolean((kind && row.dataset.kind !== kind) || !row.dataset.name.toLowerCase().includes(query));
      if (row.hidden) row.querySelector('input').checked = false;
    });
    update();
  }
  document.getElementById('storage-kind').addEventListener('change', filter);
  document.getElementById('storage-search').addEventListener('input', filter);
  document.querySelectorAll('.storage-select').forEach(input => input.addEventListener('change', update));
  selectAll.addEventListener('change', () => {
    rows.filter(r => !r.hidden && !r.querySelector('input').disabled).forEach(r => { r.querySelector('input').checked = selectAll.checked; });
    update();
  });
  async function remove(url, data) {
    busy = true; update();
    try {
      const response = await fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json', 'X-CSRF-Token': config.csrf}, body: JSON.stringify(data)});
      if (response.status === 401) { location.href = config.login_url; return; }
      const report = await response.json().catch(() => ({error: '页面验证已失效，请刷新后重试。'}));
      if (!response.ok) throw new Error(report.error || '删除未完成，请刷新后重试。');
      sessionStorage.setItem('admin-storage-message', `已删除 ${report.deleted.length} 项，释放 ${formatBytes(report.bytes)}；跳过 ${report.skipped.length} 项，失败 ${report.errors.length} 项。` + (report.errors.length ? ' 请检查文件权限并刷新列表。' : ''));
      location.reload();
    } catch (error) {
      message.textContent = error.message; message.className = 'alert error'; message.hidden = false;
    } finally { busy = false; update(); }
  }
  deleteButton.addEventListener('click', () => {
    const items = selected();
    if (!items.length || busy) return;
    const names = items.slice(0, 5).map(r => r.dataset.name).join('\n') + (items.length > 5 ? '\n…' : '');
    if (confirm(`确定永久删除所选 ${items.length} 项及其配套文件？\n${names}\n\n删除后无法恢复。`)) remove(config.delete_url, {keys: items.map(r => r.dataset.key)});
  });
  cleanupButton.addEventListener('click', () => {
    if (!busy && confirm('确定删除超过 3 个月的分析结果？进行中的任务会跳过，删除后无法恢复。')) remove(config.cleanup_url, {});
  });
  const saved = sessionStorage.getItem('admin-storage-message');
  if (saved) { message.textContent = saved; message.hidden = false; sessionStorage.removeItem('admin-storage-message'); }
  update();
})();
