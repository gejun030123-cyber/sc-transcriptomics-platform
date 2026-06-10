# routes/chat.py
"""AI 对话 API"""
from flask import Blueprint, request, jsonify
from config import Config

chat_bp = Blueprint('chat', __name__)

# 内存中的聊天历史（按项目 ID 存储，重启清空）
_chat_histories = {}


@chat_bp.route('/api/chat', methods=['POST'])
def chat_endpoint():
    """发送消息给 AI 助手"""
    if not Config.AI_API_KEY:
        return jsonify({"error": "AI API Key 未配置，请在 config.py 或环境变量 AI_API_KEY 中设置"}), 400

    data = request.get_json(silent=True) or {}
    message = data.get('message', '').strip()
    project_id = data.get('project_id', '')

    if not message:
        return jsonify({"error": "消息不能为空"}), 400

    if not project_id:
        return jsonify({"error": "需要指定项目 ID"}), 400

    # 获取或创建聊天历史
    if project_id not in _chat_histories:
        _chat_histories[project_id] = []

    history = _chat_histories[project_id]
    history.append({"role": "user", "content": message})

    # 调用 AI
    try:
        from modules.ai_adapter import chat as ai_chat
        result = ai_chat(history, project_id=project_id)

        # 更新历史
        _chat_histories[project_id] = result["messages"]

        return jsonify({
            "reply": result["reply"],
            "tool_calls": result["tool_calls"],
        })
    except Exception as e:
        return jsonify({"error": f"AI 调用失败: {str(e)}"}), 500


@chat_bp.route('/api/chat/history/<pid>')
def chat_history(pid):
    """获取聊天历史"""
    history = _chat_histories.get(pid, [])
    messages = []
    for msg in history:
        if msg.get("role") in ("user", "assistant"):
            messages.append({
                "role": msg["role"],
                "content": msg.get("content", ""),
            })
    return jsonify({"messages": messages})


@chat_bp.route('/api/chat/clear/<pid>', methods=['POST'])
def clear_history(pid):
    """清空聊天历史"""
    _chat_histories.pop(pid, None)
    return jsonify({"status": "ok"})
