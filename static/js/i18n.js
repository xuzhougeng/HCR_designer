/* Shared with Jinja through the catalog embedded by base.html. */
(() => {
  'use strict';
  const config = JSON.parse(document.getElementById('i18n-config').textContent);
  const t = (source, values = {}) => {
    if (source === null || source === undefined) return '';
    const text = Object.prototype.hasOwnProperty.call(config.messages, source) ? config.messages[source] : source;
    return String(text).replace(/\{(\w+)\}/g, (match, key) =>
      Object.prototype.hasOwnProperty.call(values, key) ? String(values[key]) : match);
  };
  window.HCR = {t, language: config.language, locale: config.language === 'en' ? 'en-US' : 'zh-CN'};
})();
