(() => {
  'use strict';
  const config = JSON.parse(document.getElementById('reference-config').textContent);
  const message = document.getElementById('reference-message');
  document.querySelectorAll('[data-bytes]').forEach(el => {
    let value = Number(el.dataset.bytes), i = 0;
    const units = ['B', 'KB', 'MB', 'GB', 'TB'];
    while (value >= 1024 && i < units.length - 1) { value /= 1024; i++; }
    el.textContent = value.toLocaleString(HCR.locale, {maximumFractionDigits:2}) + ' ' + units[i];
  });
  const saved = sessionStorage.getItem('reference-library-message');
  if (saved) { HCR.setAlert(message, saved); sessionStorage.removeItem('reference-library-message'); }
  const entries = [...document.querySelectorAll('.reference-entry')];
  const search = document.getElementById('reference-search');
  const pageSize = document.getElementById('reference-page-size');
  const pagers = [...document.querySelectorAll('.reference-pagination')];
  const view = new URL(location.href);
  search.value = view.searchParams.get('q') || '';
  pageSize.value = ['10', '20', '50'].includes(view.searchParams.get('page_size')) ? view.searchParams.get('page_size') : '10';
  let currentPage = Number(view.searchParams.get('page')) || 1;
  if (!Number.isSafeInteger(currentPage)) currentPage = 1;

  function renderPage() {
    const query = search.value.trim().toLowerCase();
    const matches = entries.filter(row => (row.dataset.filename + ' ' + row.dataset.display).toLowerCase().includes(query));
    const size = Number(pageSize.value);
    const pages = Math.max(1, Math.ceil(matches.length / size));
    currentPage = Math.min(Math.max(1, currentPage), pages);
    const start = (currentPage - 1) * size;
    const visible = new Set(matches.slice(start, start + size));
    entries.forEach(row => { row.hidden = !visible.has(row); });
    document.getElementById('reference-no-results').hidden = !entries.length || !!matches.length;
    const count = query ? HCR.t("匹配 {v0} / {v1} 个", {v0: matches.length, v1: entries.length}) : HCR.t("共 {v0} 个", {v0: entries.length});
    const summary = matches.length ? HCR.t("显示 {v0}–{v1} 个 · {v2} · 第 {v3} / {v4} 页", {v0: start + 1, v1: Math.min(start + size, matches.length), v2: count, v3: currentPage, v4: pages}) : count;
    pagers.forEach((pager, index) => {
      pager.hidden = !entries.length || (index > 0 && pages === 1);
      pager.querySelector('[data-page-summary]').textContent = summary;
      pager.querySelector('[data-page-step="-1"]').disabled = currentPage === 1;
      pager.querySelector('[data-page-step="1"]').disabled = currentPage === pages;
      const numbers = pager.querySelector('[data-page-numbers]');
      numbers.replaceChildren();
      let previous = 0;
      for (let page = 1; page <= pages; page++) {
        if (pages > 7 && page !== 1 && page !== pages && Math.abs(page - currentPage) > 1) continue;
        if (page > previous + 1) {
          const gap = document.createElement('span');
          gap.textContent = '…';
          numbers.append(gap);
        }
        const button = document.createElement('button');
        button.type = 'button';
        button.className = 'btn small';
        button.dataset.page = page;
        button.textContent = page;
        button.setAttribute('aria-label', HCR.t("第 {v0} 页", {v0: page}));
        if (page === currentPage) button.setAttribute('aria-current', 'page');
        numbers.append(button);
        previous = page;
      }
    });
    // Keep the current view when saving settings, computing hashes or deleting a reference reloads the page.
    const url = new URL(location.href);
    [['q', search.value.trim()], ['page_size', size === 10 ? '' : size], ['page', currentPage === 1 ? '' : currentPage]].forEach(([key, value]) => {
      if (value) url.searchParams.set(key, value); else url.searchParams.delete(key);
    });
    history.replaceState(null, '', url);
  }
  search.addEventListener('input', () => { currentPage = 1; renderPage(); });
  pageSize.addEventListener('change', () => { currentPage = 1; renderPage(); });
  pagers.forEach(pager => pager.addEventListener('click', event => {
    const button = event.target.closest('button');
    if (!button || button.disabled) return;
    currentPage = button.dataset.page ? Number(button.dataset.page) : currentPage + Number(button.dataset.pageStep);
    renderPage();
    pagers[0].querySelector('[aria-current="page"]').focus({preventScroll: true});
    document.getElementById('reference-list').scrollIntoView({block: 'start'});
  }));
  renderPage();
  async function write(url, data, success) {
    document.querySelectorAll('form button').forEach(button => { button.disabled = true; });
    HCR.setAlert(message, HCR.t("正在处理参考文件，请稍候……"), 'busy');
    try {
      const isForm = data instanceof FormData;
      const response = await fetch(url, {method:'POST', headers:{'X-CSRF-Token':config.csrf, ...(isForm ? {} : {'Content-Type':'application/json'})}, body:isForm ? data : JSON.stringify(data)});
      if (response.status === 401) { location.href = config.login_url; return; }
      const result = await response.json();
      if (!response.ok) throw new Error(result.error || HCR.t("参考库操作未完成。"));
      const text = result.errors?.length ? HCR.t("部分文件未能删除，请检查权限后刷新。") : success(result);
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
      write(config.upload_url, data, result => result.duplicate ? HCR.t("内容重复，已复用参考库“{v0}”（{v1}）。已有配置与归属已保留，如需更新请在该参考库下保存配置。", {v0: result.reference.display_name, v1: result.reference.filename}) : HCR.t("参考库“{v0}”已上传，SHA-256 已计算。", {v0: result.reference.display_name}));
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
      write(config.configure_url, data, () => HCR.t("参考库配置已保存。"));
    });
    row.querySelector('.compute-hash').addEventListener('click', () => write(config.hash_url, {filename:row.dataset.filename}, () => HCR.t("SHA-256 与序列数已更新。")));
    row.querySelector('.delete-reference')?.addEventListener('click', () => {
      if (confirm(HCR.t("确定永久删除“{v0}”及其全部 BLAST 索引？\n删除后无法恢复，在用参考会受到保护。", {v0: row.dataset.display}))) write(config.delete_url, {keys:[row.dataset.filename]}, () => HCR.t("参考库及配套索引已删除。"));
    });
  });
})();
