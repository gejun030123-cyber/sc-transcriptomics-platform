# routes/chat.py
"""AI 对话 API"""
from flask import Blueprint, request, jsonify
from config import Config
from routes.auth import require_ai_token

chat_bp = Blueprint('chat', __name__)


MAX_HISTORY_PER_PROJECT = 50


class ChatHistoryStore:
    """聊天历史存储，每个项目最多保留 MAX_HISTORY_PER_PROJECT 条消息。"""
    def __init__(self, max_per_project=MAX_HISTORY_PER_PROJECT):
        self._store = {}
        self._max = max_per_project

    def get(self, project_id):
        return self._store.get(project_id, [])

    def set(self, project_id, messages):
        if len(messages) > self._max:
            messages = messages[-self._max:]
        self._store[project_id] = messages

    def clear(self, project_id):
        self._store.pop(project_id, None)


_chat_histories = ChatHistoryStore()


@chat_bp.route('/api/chat', methods=['POST'])
@require_ai_token
def chat_endpoint():
    """发送消息给 AI 助手"""
    if not Config.AI_API_KEY:
        return jsonify({"error": "AI API Key 未配置，请在 config.py 或环境变量 AI_API_KEY 中设置"}), 400

    data = request.get_json(silent=True) or {}
    message = data.get('message', '').strip()
    project_id = data.get('project_id', '')

    if not message:
        return jsonify({"error": "消息不能为空"}), 400

    if len(message) > 4000:
        return jsonify({"error": "消息过长（最多 4000 字符），请缩短后重试"}), 400

    if not project_id:
        return jsonify({"error": "需要指定项目 ID"}), 400

    # 获取聊天历史
    history = _chat_histories.get(project_id)
    history.append({"role": "user", "content": message})

    # 调用 AI
    try:
        from modules.ai_adapter import chat as ai_chat
        result = ai_chat(history, project_id=project_id)

        # 更新历史
        _chat_histories.set(project_id, result["messages"])

        return jsonify({
            "reply": result["reply"],
            "tool_calls": result["tool_calls"],
            "proposed_tools": result.get("proposed_tools", []),
        })
    except Exception as e:
        return jsonify({"error": f"AI 调用失败: {str(e)}"}), 500


@chat_bp.route('/api/chat/history/<pid>')
@require_ai_token
def chat_history(pid):
    """获取聊天历史"""
    history = _chat_histories.get(pid)
    messages = []
    for msg in history:
        if msg.get("role") in ("user", "assistant"):
            messages.append({
                "role": msg["role"],
                "content": msg.get("content", ""),
            })
    return jsonify({"messages": messages})


@chat_bp.route('/api/chat/clear/<pid>', methods=['POST'])
@require_ai_token
def clear_history(pid):
    """清空聊天历史"""
    _chat_histories.clear(pid)
    return jsonify({"status": "ok"})


@chat_bp.route('/api/chat/approve', methods=['POST'])
@require_ai_token
def approve_tool():
    """用户确认执行 AI 提议的工具调用"""
    data = request.get_json(silent=True) or {}
    tool_name = data.get('tool_name', '')
    args = data.get('args', {})
    project_id = data.get('project_id', '')

    if not tool_name:
        return jsonify({"error": "缺少 tool_name"}), 400
    if not project_id:
        return jsonify({"error": "缺少 project_id"}), 400

    try:
        from modules.ai_tools import execute_tool
        result = execute_tool(tool_name, args, project_id)
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": f"执行失败: {str(e)}"}), 500
