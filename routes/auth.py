"""Account sessions and project ownership at the HTTP boundary."""
import hmac
import threading
import time
from datetime import timedelta
from functools import wraps
from urllib.parse import urlsplit

from flask import abort, g, jsonify, redirect, render_template, request, session, url_for

from config import Config
from models import User, Project

_login_failures = {}
_login_lock = threading.Lock()


def require_ai_token(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if Config.AI_API_TOKEN:
            auth = request.headers.get('Authorization', '')
            if not auth.startswith('Bearer ') or not hmac.compare_digest(auth[7:], Config.AI_API_TOKEN):
                return jsonify({'error': '认证失败，无效的 API Token'}), 401
        return f(*args, **kwargs)
    return decorated


def _safe_next_path(value):
    if not value:
        return '/'
    parsed = urlsplit(value)
    if parsed.scheme or parsed.netloc or not value.startswith('/') or value.startswith('//'):
        return '/'
    return value


def current_user():
    return getattr(g, 'current_user', None)


def project_is_owned(pid):
    user = current_user()
    project = Project.get_by_id(pid) if pid else None
    return bool(user and project and project.owner_user_id == user.id)


def _same_origin_write():
    origin = request.headers.get('Origin') or request.headers.get('Referer')
    if origin:
        parsed = urlsplit(origin)
        return parsed.scheme == request.scheme and parsed.netloc == request.host
    return request.headers.get('Sec-Fetch-Site') == 'same-origin'


def _project_for_path():
    """Recognize project paths, including the chat API and query based APIs."""
    values = request.view_args or {}
    pid = values.get('pid') or values.get('project_id')
    if pid:
        return pid
    if request.endpoint in {'chat.chat_endpoint', 'chat.approve_tool'}:
        body = request.get_json(silent=True) or {}
        return body.get('project_id')
    if request.endpoint in {'api.list_presets', 'api.get_preset', 'api.delete_preset'}:
        return request.args.get('project_id')
    if request.endpoint == 'api.save_preset':
        body = request.get_json(silent=True) or {}
        return body.get('project_id')
    return None


def _object_project():
    values = request.view_args or {}
    from models import AnalysisTask, ResultFile, PipelineRun
    if values.get('task_id'):
        item = AnalysisTask.get_by_id(values['task_id'])
        return item.project_id if item else None
    if values.get('file_id'):
        item = ResultFile.get_by_id(values['file_id'])
        return item.project_id if item else None
    if values.get('run_id'):
        item = PipelineRun.get_by_id(values['run_id'])
        return item.project_id if item else None
    return None


def register_platform_access_gate(app):
    app.permanent_session_lifetime = timedelta(hours=Config.PLATFORM_ACCESS_SESSION_HOURS)

    @app.context_processor
    def account_context():
        return {'signed_in_user': current_user()}

    @app.route('/login', methods=['GET', 'POST'])
    def platform_login():
        next_path = _safe_next_path(request.args.get('next') or request.form.get('next'))
        error = None
        if request.method == 'POST':
            if not _same_origin_write():
                abort(403)
            identifier = request.form.get('username', '').strip()
            password = request.form.get('password', '')
            key = (request.remote_addr or '', identifier.casefold())
            now = time.monotonic()
            with _login_lock:
                recent = [stamp for stamp in _login_failures.get(key, ()) if now - stamp < 900]
                _login_failures[key] = recent
                limited = len(recent) >= 5
            if limited:
                return render_template('login.html', error='尝试次数过多，请 15 分钟后重试。', next_path=next_path), 429
            user = User.get_by_login(identifier)
            if user and user.is_active and user.check_password(password):
                with _login_lock:
                    _login_failures.pop(key, None)
                session.clear()
                session.permanent = True
                session['user_id'] = user.id
                session['session_version'] = user.session_version
                return redirect(next_path)
            with _login_lock:
                _login_failures.setdefault(key, []).append(now)
            error = '用户名或密码不正确。'
        return render_template('login.html', error=error, next_path=next_path)

    @app.route('/logout', methods=['POST'])
    def platform_logout():
        if not _same_origin_write():
            abort(403)
        session.clear()
        return redirect(url_for('platform_login'))

    @app.before_request
    def enforce_platform_access():
        g.current_user = None
        # Existing unit tests explicitly opt out; production has no anonymous mode.
        if not app.config.get('AUTH_REQUIRED', True):
            return None
        if request.endpoint == 'static':
            return None
        if request.endpoint in {'platform_login', 'accounts.register',
                                'accounts.forgot_password', 'accounts.reset_password'}:
            return None
        uid = session.get('user_id')
        user = User.get_by_id(uid) if uid else None
        if not user or not user.is_active or user.session_version != session.get('session_version'):
            session.clear()
            if request.path.startswith('/api/'):
                return jsonify({'error': '请先登录'}), 401
            return redirect(url_for('platform_login', next=request.full_path))
        g.current_user = user
        if request.method not in {'GET', 'HEAD', 'OPTIONS'} and not _same_origin_write():
            abort(403)
        pid = _project_for_path()
        if pid and not project_is_owned(pid):
            abort(404)
        object_pid = _object_project()
        if object_pid and not project_is_owned(object_pid):
            abort(404)
        if pid and object_pid and pid != object_pid:
            abort(404)
        return None
