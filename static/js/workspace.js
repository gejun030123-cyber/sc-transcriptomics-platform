(function () {
    'use strict';

    const cfg = window.WORKSPACE_CONFIG || {};
    const projectId = cfg.projectId || '';
    let busy = false;
    let pendingTools = [];
    let activeAssay = 'scrna';

    const $ = (selector) => document.querySelector(selector);
    const $$ = (selector) => Array.from(document.querySelectorAll(selector));

    function escapeHtml(value) {
        const node = document.createElement('div');
        node.textContent = value == null ? '' : String(value);
        return node.innerHTML.replace(/\n/g, '<br>');
    }

    function formatToolArgs(args) {
        try {
            return JSON.stringify(args || {}, null, 2);
        } catch (error) {
            return String(args || '');
        }
    }

    function scrollMessages() {
        const box = $('#workspace-messages');
        if (box) box.scrollTop = box.scrollHeight;
    }

    function addAttachmentCards(container, attachments) {
        if (!container || !Array.isArray(attachments) || !attachments.length) return;
        const gallery = document.createElement('div');
        gallery.className = 'workspace-message-attachments';
        attachments.forEach((attachment) => {
            if (!attachment || typeof attachment.url !== 'string' || !attachment.url.startsWith('/api/projects/')) return;
            const card = document.createElement('div');
            card.className = 'workspace-attachment-card';
            const link = document.createElement('a');
            link.href = attachment.url;
            link.target = '_blank';
            link.rel = 'noopener';
            link.title = '打开原图';
            const image = document.createElement('img');
            image.src = attachment.url;
            image.alt = attachment.label || '分析结果图';
            image.loading = 'lazy';
            link.appendChild(image);
            const caption = document.createElement('div');
            caption.className = 'workspace-attachment-caption';
            caption.textContent = attachment.label || attachment.category || '分析结果图';
            const download = document.createElement('a');
            download.href = attachment.url;
            download.download = '';
            download.className = 'workspace-attachment-download';
            download.textContent = '下载';
            card.appendChild(link);
            card.appendChild(caption);
            card.appendChild(download);
            gallery.appendChild(card);
        });
        if (gallery.children.length) container.appendChild(gallery);
    }

    function addMessage(role, content, attachments) {
        const box = $('#workspace-messages');
        if (!box) return;
        const wrapper = document.createElement('div');
        wrapper.className = 'workspace-message ' + (role || 'assistant');
        const label = role === 'user' ? '你' : (role === 'system' ? '系统' : 'AI');
        wrapper.innerHTML = '<div class="workspace-message-meta">' + label + '</div>'
            + '<div class="workspace-message-bubble">' + escapeHtml(content || '') + '</div>';
        addAttachmentCards(wrapper, attachments);
        box.appendChild(wrapper);
        scrollMessages();
    }

    function addToolCall(name, args) {
        const box = $('#workspace-messages');
        if (!box) return;
        const node = document.createElement('div');
        node.className = 'workspace-tool-call';
        node.textContent = '已读取：' + name + (args && Object.keys(args).length
            ? '\n' + formatToolArgs(args) : '');
        box.appendChild(node);
        scrollMessages();
    }

    function toolDescription(name) {
        const descriptions = {
            run_analysis: '提交一个分析模块，需要你确认后执行。',
            run_pipeline: '提交一条完整分析流程，后台会按顺序执行。',
            propose_parameter_sweep: '生成候选参数方案，不会直接运行。',
            run_parameter_sweep: '启动参数搜索任务，可能需要较长时间。',
            start_goal_agent: '创建目标驱动的 Agent 会话。',
            continue_goal_agent: '继续已有目标，可能生成或采纳候选分支。',
        };
        return descriptions[name] || '该操作会改变项目状态，需要确认。';
    }

    function addProposedTool(tool, index) {
        const box = $('#workspace-messages');
        if (!box) return;
        const node = document.createElement('div');
        node.className = 'workspace-proposed-tool';
        node.dataset.toolIndex = String(index);
        node.innerHTML = '<div class="workspace-proposed-title">AI 建议执行：'
            + escapeHtml(tool.name || '未知工具') + '</div>'
            + '<div class="workspace-proposed-desc">' + escapeHtml(tool.description || toolDescription(tool.name)) + '</div>'
            + '<pre class="workspace-proposed-params">' + escapeHtml(formatToolArgs(tool.args)) + '</pre>'
            + '<div class="workspace-proposed-actions">'
            + '<button type="button" class="btn btn-primary btn-sm" data-approve-tool="' + index + '">确认执行</button>'
            + '<button type="button" class="btn btn-outline-secondary btn-sm" data-reject-tool="' + index + '">取消</button>'
            + '</div>';
        box.appendChild(node);
        scrollMessages();
    }

    function showTyping() {
        hideTyping();
        const box = $('#workspace-messages');
        if (!box) return;
        const node = document.createElement('div');
        node.id = 'workspace-typing';
        node.className = 'workspace-typing';
        node.textContent = 'AI 正在读取项目上下文并组织下一步…';
        box.appendChild(node);
        scrollMessages();
    }

    function hideTyping() {
        const node = $('#workspace-typing');
        if (node) node.remove();
    }

    function setBusy(value) {
        busy = value;
        const input = $('#workspace-input');
        const button = $('#workspace-send');
        if (input) input.disabled = value;
        if (button) button.disabled = value;
    }

    function setModelStatus(text, kind) {
        const node = $('#workspace-model-status');
        if (node) node.textContent = text;
        const dot = $('#workspace-status-dot');
        if (dot) {
            // Keep status presentation within the workspace's single blue
            // visual scale; the adjacent text carries the actual state.
            dot.style.background = kind === 'error' ? '#6f8fb9' : '#74aaff';
        }
    }

    async function loadConfig() {
        try {
            const response = await fetch('/api/chat/config');
            const data = await response.json();
            if (!response.ok || data.error) {
                setModelStatus(data.error || 'AI 配置不可用', 'error');
                return;
            }
            if (!data.configured) {
                setModelStatus('未配置 API Key · 请打开 AI API 设置', 'error');
                const unconfiguredBadge = $('#workspace-model-badge');
                if (unconfiguredBadge) unconfiguredBadge.textContent = '待配置 API';
                return;
            }
            setModelStatus((data.model || '未命名模型') + ' · ' + (data.provider || 'provider'));
            const badge = $('#workspace-model-badge');
            if (badge) badge.textContent = data.model || '当前模型';
        } catch (error) {
            setModelStatus('AI 状态读取失败', 'error');
        }
    }

    async function loadHistory() {
        try {
            const response = await fetch('/api/chat/history/' + encodeURIComponent(projectId));
            const data = await response.json();
            if (!response.ok || !Array.isArray(data.messages) || !data.messages.length) return;
            const box = $('#workspace-messages');
            if (!box) return;
            box.innerHTML = '';
            data.messages.forEach((message) => addMessage(message.role, message.content || '', message.attachments || []));
        } catch (error) {
            addMessage('system', '历史会话加载失败：' + error.message);
        }
    }

    async function loadContext() {
        try {
            const response = await fetch('/api/projects/' + encodeURIComponent(projectId) + '/current-context');
            const data = await response.json();
            const source = $('#workspace-baseline');
            if (!source) return;
            if (!data || data.source === 'none') {
                source.textContent = '尚未建立分析基线';
                return;
            }
            const labels = {
                accepted_branch: '已采纳候选分支',
                main_task: '主线任务输出',
                intermediate: '中间 h5ad',
                upload: '上传文件',
            };
            source.textContent = labels[data.source] || data.source;
        } catch (error) {
            const source = $('#workspace-baseline');
            if (source) source.textContent = '上下文读取失败';
        }
    }

    async function loadAgentState() {
        try {
            const response = await fetch('/api/projects/' + encodeURIComponent(projectId) + '/agent/jobs');
            const data = await response.json();
            const node = $('#workspace-agent-jobs');
            if (!node) return;
            const jobs = (data.jobs || []).filter((job) => ['pending', 'running'].includes(job.status));
            node.textContent = jobs.length ? ('Agent 运行中：' + jobs.length + ' 个任务') : '暂无 Agent 后台任务';
        } catch (error) {
            const node = $('#workspace-agent-jobs');
            if (node) node.textContent = 'Agent 状态暂不可用';
        }
    }

    function workspaceContext() {
        return {
            workspace: 'ai_main',
            active_assay: activeAssay,
            available_assays: ['scrna', 'bulk_rna'],
            planned_assays: ['wes', 'bulk_atac'],
        };
    }

    async function sendMessage(prefilled) {
        if (busy) return;
        const input = $('#workspace-input');
        const message = (prefilled || (input && input.value) || '').trim();
        if (!message) return;
        if (input) {
            input.value = '';
            input.style.height = 'auto';
        }
        addMessage('user', message);
        showTyping();
        setBusy(true);
        try {
            const response = await fetch('/api/chat', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    message: message,
                    project_id: projectId,
                    context: workspaceContext(),
                }),
            });
            const data = await response.json();
            hideTyping();
            setBusy(false);
            if (!response.ok || data.error) {
                setModelStatus('AI 调用失败', 'error');
                addMessage('assistant', '执行失败：' + (data.error || '未知错误'));
                return;
            }
            if (data.config && data.config.model) {
                setModelStatus(data.config.model + ' · ' + (data.config.message_count || 0) + ' 条上下文');
            }
            (data.tool_calls || []).forEach((call) => addToolCall(call.name, call.args || {}));
            $$('.workspace-proposed-tool').forEach((node) => node.remove());
            pendingTools = data.proposed_tools || [];
            pendingTools.forEach((tool, index) => addProposedTool(tool, index));
            if (data.reply || (data.attachments && data.attachments.length)) {
                addMessage('assistant', data.reply || '分析结果图片已生成：', data.attachments || []);
            }
            loadContext();
            loadAgentState();
        } catch (error) {
            hideTyping();
            setBusy(false);
            setModelStatus('网络错误', 'error');
            addMessage('assistant', '网络错误：' + error.message);
        }
    }

    async function approveTool(index) {
        const tool = pendingTools[index];
        if (!tool) return;
        const card = document.querySelector('[data-tool-index="' + index + '"]');
        const actions = card && card.querySelector('.workspace-proposed-actions');
        if (actions) actions.innerHTML = '<span class="text-muted small">执行中…</span>';
        try {
            const response = await fetch('/api/chat/approve', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    tool_name: tool.name,
                    args: tool.args || {},
                    project_id: projectId,
                }),
            });
            const data = await response.json();
            if (!response.ok || data.error) {
                addMessage('assistant', '工具执行失败：' + (data.error || '未知错误'));
                if (actions) actions.innerHTML = '<button type="button" class="btn btn-primary btn-sm" data-approve-tool="' + index + '">重试</button>';
                return;
            }
            addMessage('assistant', data.message || ('已提交：' + tool.name));
            if (card) card.remove();
            loadContext();
            loadAgentState();
        } catch (error) {
            addMessage('assistant', '工具执行网络错误：' + error.message);
        }
    }

    async function clearChat() {
        try {
            await fetch('/api/chat/clear/' + encodeURIComponent(projectId), {method: 'POST'});
            const box = $('#workspace-messages');
            if (box) box.innerHTML = '';
            addMessage('assistant', '对话已清空。可以直接说分析目标，我会先检查项目状态，再给出可确认的执行计划。');
            pendingTools = [];
        } catch (error) {
            addMessage('system', '清空失败：' + error.message);
        }
    }

    function newConversation() {
        clearChat().finally(() => {
            const input = $('#workspace-input');
            if (input) input.focus();
        });
    }

    function autosize() {
        const input = $('#workspace-input');
        if (!input) return;
        input.style.height = 'auto';
        input.style.height = Math.min(input.scrollHeight, 130) + 'px';
    }

    function selectAssay(key) {
        activeAssay = key;
        $$('.workspace-card').forEach((card) => card.classList.toggle('selected', card.dataset.assay === key));
        const card = document.querySelector('[data-assay="' + key + '"]');
        const prompt = card && card.dataset.prompt;
        const input = $('#workspace-input');
        if (input && prompt) {
            input.value = prompt;
            autosize();
            input.focus();
        }
        const label = $('#workspace-active-assay');
        if (label && card) label.textContent = card.dataset.title || key;
    }

    function switchContextTab(tab) {
        $$('.workspace-context-tab').forEach((button) => button.classList.toggle('active', button.dataset.tab === tab));
        $$('.workspace-context-view').forEach((view) => view.classList.toggle('active', view.dataset.view === tab));
    }

    function bindEvents() {
        const input = $('#workspace-input');
        if (input) {
            input.addEventListener('input', autosize);
            input.addEventListener('keydown', (event) => {
                if (event.key === 'Enter' && !event.shiftKey) {
                    event.preventDefault();
                    sendMessage();
                }
            });
        }
        document.addEventListener('click', (event) => {
            const quick = event.target.closest('[data-workspace-prompt]');
            if (quick) {
                event.preventDefault();
                sendMessage(quick.dataset.workspacePrompt);
                return;
            }
            const assay = event.target.closest('[data-select-assay]');
            if (assay) {
                event.preventDefault();
                selectAssay(assay.dataset.selectAssay);
                return;
            }
            const approve = event.target.closest('[data-approve-tool]');
            if (approve) {
                event.preventDefault();
                approveTool(Number(approve.dataset.approveTool));
                return;
            }
            const reject = event.target.closest('[data-reject-tool]');
            if (reject) {
                event.preventDefault();
                const card = reject.closest('.workspace-proposed-tool');
                if (card) card.remove();
                return;
            }
            const tab = event.target.closest('[data-context-tab]');
            if (tab) {
                event.preventDefault();
                switchContextTab(tab.dataset.contextTab);
            }
        });
        const clear = $('#workspace-clear');
        if (clear) clear.addEventListener('click', clearChat);
        const newChat = $('#workspace-new-chat');
        if (newChat) newChat.addEventListener('click', newConversation);
    }

    window.workspaceSend = sendMessage;
    window.workspaceSelectAssay = selectAssay;

    document.addEventListener('DOMContentLoaded', () => {
        bindEvents();
        loadConfig();
        loadHistory();
        loadContext();
        loadAgentState();
        setInterval(loadContext, 15000);
        setInterval(loadAgentState, 15000);
    });
}());
