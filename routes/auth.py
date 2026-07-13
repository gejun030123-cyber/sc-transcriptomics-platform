# routes/auth.py
"""共享鉴权装饰器 — 供 chat、branches 等路由复用"""
from functools import wraps
from flask import request, jsonify
from config import Config


def require_ai_token(f):
    """API Token 认证装饰器。AI_API_TOKEN 为空时跳过认证。"""
    @wraps(f)
    def decorated(*args, **kwargs):
        if Config.AI_API_TOKEN:
            auth = request.headers.get('Authorization', '')
            if not auth.startswith('Bearer ') or auth[7:] != Config.AI_API_TOKEN:
                return jsonify({"error": "认证失败，无效的 API Token"}), 401
        return f(*args, **kwargs)
    return decorated
