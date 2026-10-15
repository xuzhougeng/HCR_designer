(() => {
  'use strict';
  const config = JSON.parse(document.getElementById('reference-config').textContent);
  const message = document.getElementById('reference-message');
  document.querySelectorAll('[data-bytes]').forEach(el => {
    let value = Number(el.dataset.bytes), i = 0;
    const units = ['B', 'KB', 'MB', 'GB', 'TB'];
    while (value >= 1024 && i < units.length - 1) { value /= 1024; i++; }
    el.textContent = value.toLocaleString('zh-CN', {maximumFractionDigits:2}) + ' ' + units[i];
  });
  const saved = sessionStorage.getItem('reference-library-message');
  if (saved) { HCR.setAlert(message, saved); sessionStorage.removeItem('reference-library-message'); }
  document.getElementById('reference-search').addEventListener('input', event => {
    const query = event.target.value.trim().toLowerCase();
    document.querySelectorAll('.reference-entry').forEach(row => { row.hidden = !(row.dataset.filename + ' ' + row.dataset.display).toLowerCase().includes(query); });
  });
  async function write(url, data, success) {
    document.querySelectorAll('form button').forEach(button => { button.disabled = true; });
    HCR.setAlert(message, '正在处理参考文件，请稍候……', 'busy');
    try {
      const isForm = data instanceof FormData;
      const response = await fetch(url, {method:'POST', headers:{'X-CSRF-Token':config.csrf, ...(isForm ? {} : {'Content-Type':'application/json'})}, body:isForm ? data : JSON.stringify(data)});
      if (response.status === 401) { location.href = config.login_url; return; }
      const result = await response.json();
      if (!response.ok) throw new Error(result.error || '参考库操作未完成。');
      const text = result.errors?.length ? '部分文件未能删除，请检查权限后刷新。' : success(result);
      sessionStorage.setItem('reference-library-message', text); location.reload();
    } catch (error) { HCR.setAlert(message, error.message, 'error'); }
    finally { document.querySelectorAll('form button').forEach(button => { button.disabled = false; }); }
  }
  const uploadForm = document.getElementById('reference-upload-form');
  if (uploadForm) {
    document.getElementById('library-file').addEventListener('change', event => {
      const file = event.target.files[0];
      if (file) document.getElementById('library-name').value = file.name.replace(/\.gz$/i, '').replace(/\.(fa|fasta|fna)$/i, '');
    });
    uploadForm.addEventListener('submit', event => {
      event.preventDefault();
      if (!uploadForm.reportValidity()) return;
      const data = new FormData();
      data.append('reference', document.getElementById('library-file').files[0]);
      data.append('display_name', document.getElementById('library-name').value);
      data.append('polyploid', document.getElementById('library-ploidy').value);
      const annotation = document.getElementById('library-subgenome-map').files[0];
      if (annotation) data.append('subgenome_map', annotation);
      write(config.upload_url, data, result => result.duplicate ? `内容重复，已复用参考库“${result.reference.display_name}”（${result.reference.filename}）。已有配置与归属已保留，如需更新请在该参考库下保存配置。` : `参考库“${result.reference.display_name}”已上传，SHA-256 已计算。`);
    });
  }
  document.querySelectorAll('.reference-entry').forEach(row => {
    const form = row.querySelector('form');
    if (!form) return;
    form.addEventListener('submit', event => {
      event.preventDefault(); if (!form.reportValidity()) return;
      const data = new FormData(form);
      data.append('filename', row.dataset.filename);
      if (!form.elements.subgenome_map.files.length) data.delete('subgenome_map');
      write(config.configure_url, data, () => '参考库配置已保存。');
    });
    row.querySelector('.compute-hash').addEventListener('click', () => write(config.hash_url, {filename:row.dataset.filename}, () => 'SHA-256 与序列数已更新。'));
    row.querySelector('.delete-reference')?.addEventListener('click', () => {
      if (confirm(`确定永久删除“${row.dataset.display}”及其全部 BLAST 索引？\n删除后无法恢复，在用参考会受到保护。`)) write(config.delete_url, {keys:[row.dataset.filename]}, () => '参考库及配套索引已删除。');
    });
  });
})();
