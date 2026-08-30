(function () {
    'use strict';

    const $ = (selector) => document.querySelector(selector);
    const form = $('#ai-settings-form');

    function showMessage(text, kind) {
        const node = $('#ai-settings-message');
        if (!node) return;
        node.className = 'alert mt-3 mb-0 alert-' + (kind || 'info');
        node.textContent = text || '';
    }

    function setBusy(button, value, label) {
        if (!button) return;
        button.disabled = value;
        if (value) {
            button.dataset.originalLabel = button.textContent;
            button.textContent = label || '处理中…';
        } else if (button.dataset.originalLabel) {
            button.textContent = button.dataset.originalLabel;
        }
    }

    function render(data) {
        if (!data) return;
        $('#ai-provider').value = data.provider || 'auto';
        $('#ai-model').value = data.model || '';
        $('#ai-api-url').value = data.api_url || '';
        $('#summary-provider').textContent = data.provider === 'auto' ? '自动识别' : (data.provider || '—');
        $('#summary-model').textContent = data.model || '—';
        $('#summary-url').textContent = data.api_url || '—';
        $('#summary-key').textContent = data.api_key_masked || '未配置';
        const state = $('#ai-settings-state');
        state.textContent = data.configured ? 'AI 已配置' : '待配置 API Key';
        state.className = 'badge rounded-pill ' + (data.configured ? 'text-bg-success' : 'text-bg-warning');
    }

    async function load() {
        try {
            const response = await fetch('/api/settings/ai');
            const data = await response.json();
            if (!response.ok || data.error) throw new Error(data.error || '读取失败');
            render(data);
        } catch (error) {
            showMessage('读取 AI 配置失败：' + error.message, 'danger');
            $('#ai-settings-state').textContent = '读取失败';
            $('#ai-settings-state').className = 'badge rounded-pill text-bg-danger';
        }
    }

    function payload() {
        return {
            provider: $('#ai-provider').value,
            model: $('#ai-model').value.trim(),
            api_url: $('#ai-api-url').value.trim(),
            api_key: $('#ai-api-key').value,
            clear_api_key: $('#ai-clear-key').checked,
        };
    }

    form.addEventListener('submit', async (event) => {
        event.preventDefault();
        const button = $('#ai-save');
        setBusy(button, true, '保存中…');
        showMessage('', 'info');
        try {
            const response = await fetch('/api/settings/ai', {
                method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(payload()),
            });
            const data = await response.json();
            if (!response.ok || data.error) throw new Error(data.error || '保存失败');
            $('#ai-api-key').value = '';
            $('#ai-clear-key').checked = false;
            render(data.settings);
            showMessage('配置已保存并立即应用到 AI 主工作台。', 'success');
        } catch (error) {
            showMessage(error.message, 'danger');
        } finally { setBusy(button, false); }
    });

    $('#ai-test').addEventListener('click', async () => {
        const button = $('#ai-test');
        setBusy(button, true, '测试中…');
        try {
            const response = await fetch('/api/settings/ai/test', {
                method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(payload()),
            });
            const data = await response.json();
            if (!response.ok || !data.ok) throw new Error(data.error || '连接测试失败');
            showMessage('连接成功（HTTP ' + data.status_code + '），可以保存此配置。', 'success');
        } catch (error) { showMessage(error.message, 'danger'); }
        finally { setBusy(button, false); }
    });

    $('#ai-reset').addEventListener('click', async () => {
        if (!window.confirm('确定移除平台页面保存的覆盖配置，恢复环境变量默认值吗？')) return;
        const button = $('#ai-reset');
        setBusy(button, true, '恢复中…');
        try {
            const response = await fetch('/api/settings/ai/reset', {method: 'POST'});
            const data = await response.json();
            if (!response.ok || data.error) throw new Error(data.error || '恢复失败');
            $('#ai-api-key').value = '';
            $('#ai-clear-key').checked = false;
            render(data.settings);
            showMessage('已恢复环境变量中的默认配置。', 'success');
        } catch (error) { showMessage(error.message, 'danger'); }
        finally { setBusy(button, false); }
    });

    load();
}());
