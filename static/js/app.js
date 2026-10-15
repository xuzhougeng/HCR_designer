/* Shared behaviour for every workbench page: mobile sidebar, file drop zones,
   segmented toggles and helpers for tools that answer with a file download. */
(() => {
  'use strict';

  function saveBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.append(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 60000);
    return url;
  }

  function filenameFrom(response, fallback) {
    const disposition = response.headers.get('Content-Disposition') || '';
    const encoded = /filename\*=UTF-8''([^;]+)/i.exec(disposition);
    if (encoded) return decodeURIComponent(encoded[1]);
    const plain = /filename="?([^";]+)"?/i.exec(disposition);
    return plain ? plain[1] : fallback;
  }

  /* POST a form and return {blob, filename, response}; errors come back as JSON or HTML. */
  async function postForFile(url, body, fallbackName) {
    const response = await fetch(url, {method: 'POST', body, headers: {'X-Requested-With': 'fetch'}});
    const type = response.headers.get('Content-Type') || '';
    if (!response.ok || type.includes('application/json') || type.includes('text/html')) {
      let message = HCR.t('请求未完成（HTTP {status}）。', {status: response.status});
      if (type.includes('application/json')) {
        const data = await response.json().catch(() => null);
        if (data && data.error) message = data.error;
      }
      throw new Error(message);
    }
    return {blob: await response.blob(), filename: filenameFrom(response, fallbackName), response};
  }

  function toDelimited(rows, delimiter) {
    const quote = value => {
      const text = String(value ?? '');
      return delimiter === ',' && /[",\n]/.test(text) ? '"' + text.replace(/"/g, '""') + '"' : text;
    };
    return rows.map(row => row.map(quote).join(delimiter)).join('\n') + '\n';
  }

  function downloadRows(rows, filename) {
    const delimiter = filename.endsWith('.csv') ? ',' : '\t';
    const bom = delimiter === ',' ? '﻿' : '';
    saveBlob(new Blob([bom + toDelimited(rows, delimiter)], {type: 'text/plain;charset=utf-8'}), filename);
  }

  function setAlert(box, text, kind) {
    if (!box) return;
    box.hidden = !text;
    box.textContent = text || '';
    box.className = 'alert' + (kind ? ' ' + kind : '');
    box.setAttribute('role', kind === 'error' ? 'alert' : 'status');
  }

  async function copyText(text, button) {
    try {
      await navigator.clipboard.writeText(text);
      if (button) {
        const label = button.textContent;
        button.textContent = HCR.t('已复制');
        setTimeout(() => { button.textContent = label; }, 1200);
      }
    } catch (error) {
      window.prompt(HCR.t('复制以下内容'), text);
    }
  }

  function initShell() {
    const shell = document.querySelector('.shell');
    const toggle = document.getElementById('nav-toggle');
    if (!shell || !toggle) return;
    toggle.addEventListener('click', () => {
      const open = shell.classList.toggle('nav-open');
      toggle.setAttribute('aria-expanded', String(open));
    });
    shell.addEventListener('click', event => {
      if (shell.classList.contains('nav-open') && !event.target.closest('.sidebar, #nav-toggle')) {
        shell.classList.remove('nav-open');
        toggle.setAttribute('aria-expanded', 'false');
      }
    });
  }

  function describeFile(file) {
    const size = file.size > 1048576 ? (file.size / 1048576).toFixed(1) + ' MB' : Math.max(1, Math.round(file.size / 1024)) + ' KB';
    return file.name + ' · ' + size;
  }

  function initDropzones() {
    document.querySelectorAll('[data-dropzone]').forEach(zone => {
      const input = zone.querySelector('input[type=file]');
      const label = zone.querySelector('[data-file-label]');
      const initial = label ? label.textContent : '';
      const refresh = () => {
        const files = [...input.files];
        zone.classList.toggle('has-file', files.length > 0);
        if (label) label.textContent = files.length ? files.map(describeFile).join('；') : initial;
      };
      // A <label> zone already opens the picker natively; only other elements need help.
      if (zone.tagName !== 'LABEL') zone.addEventListener('click', event => { if (event.target !== input) input.click(); });
      zone.addEventListener('keydown', event => {
        if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); input.click(); }
      });
      zone.addEventListener('dragover', event => { event.preventDefault(); zone.classList.add('is-over'); });
      zone.addEventListener('dragleave', () => zone.classList.remove('is-over'));
      zone.addEventListener('drop', event => {
        event.preventDefault();
        zone.classList.remove('is-over');
        if (event.dataTransfer.files.length) { input.files = event.dataTransfer.files; input.dispatchEvent(new Event('change', {bubbles: true})); }
      });
      input.addEventListener('change', refresh);
      zone.resetFile = () => { input.value = ''; refresh(); };
    });
  }

  function initCopyButtons() {
    document.addEventListener('click', event => {
      const button = event.target.closest('[data-copy]');
      if (button) copyText(button.dataset.copy, button);
    });
  }

  /* Order list: probe pairs picked across jobs, kept in this browser only. */
  const ORDER_KEY = 'hcr.orderList.v1';
  function orderStateLabel(item) {
    if (item.state === 'unknown' && (item.quality_retained || ['未评估 · 按质量保留', 'Not assessed · retained by quality'].includes(item.state_label))) {
      return HCR.t('未评估 · 按质量保留');
    }
    return {keep: HCR.t('保留待复核'), excluded: HCR.t('按目标策略排除'), other: HCR.t('未保留'), unknown: HCR.t('特异性未评估')}[item.state] || item.state_label;
  }

  function orderClosestLabel(item) {
    if (item.closest_status === 'not_detected' || ['未检出', 'Not detected'].includes(item.closest)) return HCR.t('未检出');
    if (item.closest_status === 'not_assessed' || ['未评估', 'Not assessed'].includes(item.closest)) return HCR.t('未评估');
    return item.closest;
  }

  const orders = {
    read() {
      try { return JSON.parse(localStorage.getItem(ORDER_KEY) || '[]'); } catch (error) { return []; }
    },
    write(items) {
      try { localStorage.setItem(ORDER_KEY, JSON.stringify(items)); } catch (error) { /* storage blocked */ }
      updateOrderBadge(items);
      document.dispatchEvent(new CustomEvent('hcr:orders', {detail: items}));
    },
    has(key) { return orders.read().some(item => item.key === key); },
    toggle(item) {
      const items = orders.read();
      const index = items.findIndex(existing => existing.key === item.key);
      if (index >= 0) items.splice(index, 1); else items.push(Object.assign({added_at: new Date().toISOString()}, item));
      orders.write(items);
      return index < 0;
    },
    remove(key) { orders.write(orders.read().filter(item => item.key !== key)); },
    clear() { orders.write([]); },
  };

  function updateOrderBadge(items) {
    const badge = document.querySelector('[data-order-count]');
    if (!badge) return;
    const count = (items || orders.read()).length;
    badge.hidden = count === 0;
    badge.textContent = count;
  }

  Object.assign(window.HCR, {saveBlob, postForFile, downloadRows, setAlert, copyText, orders, orderStateLabel, orderClosestLabel});
  window.addEventListener('storage', event => { if (event.key === ORDER_KEY) updateOrderBadge(); });
  document.addEventListener('DOMContentLoaded', () => { initShell(); initDropzones(); initCopyButtons(); updateOrderBadge(); });
})();
