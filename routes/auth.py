# routes/auth.py
"""共享鉴权功能。

``require_ai_token`` protects the AI-specific API routes.  The optional
platform access gate below protects the whole web application when it is
published through an external tunnel such as NAT123.
"""
import hmac
from functools import wraps
from datetime import timedelta
from urllib.parse import urlsplit

from flask import (
    jsonify,
    redirect,
    render_template,
    request,
    session,
    url_for,
)
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


def _safe_next_path(value):
    """Return a local redirect target and reject external redirect URLs."""
    if not value:
        return '/'
    parsed = urlsplit(value)
    if parsed.scheme or parsed.netloc or not value.startswith('/') or value.startswith('//'):
        return '/'
    return value


def register_platform_access_gate(app):
    """Register an optional shared-password gate for all application routes.

    The gate is disabled when ``PLATFORM_ACCESS_PASSWORD`` is empty so local
    development and the existing test suite retain their current behavior.
    When enabled, browser routes redirect to ``/login`` while API routes get a
    JSON 401 response, which keeps frontend fetch calls predictable.
    """
    app.permanent_session_lifetime = timedelta(
        hours=Config.PLATFORM_ACCESS_SESSION_HOURS
    )

    @app.route('/login', methods=['GET', 'POST'])
    def platform_login():
        next_path = _safe_next_path(
            request.args.get('next') or request.form.get('next')
        )
        error = None
        if request.method == 'POST':
            submitted = request.form.get('password', '')
            expected = Config.PLATFORM_ACCESS_PASSWORD
            if expected and hmac.compare_digest(submitted, expected):
                # Rotate the signed session contents on successful login and
                # retain only the access marker needed by this gate.
                session.clear()
                session.permanent = True
                session['platform_access_granted'] = True
                return redirect(next_path)
            error = '访问密码不正确。'
        return render_template('login.html', error=error, next_path=next_path)

    @app.route('/logout', methods=['GET', 'POST'])
    def platform_logout():
        session.clear()
        return redirect(url_for('platform_login'))

    @app.before_request
    def enforce_platform_access():
        # Empty password deliberately means the gate is disabled.  The
        # deployment instructions require setting it before NAT123 exposure.
        if not Config.PLATFORM_ACCESS_PASSWORD:
            return None
        if request.endpoint in {'platform_login', 'platform_logout', 'static'}:
            return None
        if session.get('platform_access_granted'):
            return None
        if request.path.startswith('/api/'):
            return jsonify({'error': '需要先登录访问平台'}), 401
        return redirect(url_for('platform_login', next=request.full_path))
