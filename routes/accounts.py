"""Self registration and small-lab account management pages."""
import hashlib
import secrets
import sqlite3
import threading
import time
from datetime import datetime, timedelta, timezone

from flask import Blueprint, abort, flash, make_response, redirect, render_template, request, session, url_for
from werkzeug.security import generate_password_hash

from database import get_conn
from models import User
from routes.auth import _same_origin_write, current_user

accounts_bp = Blueprint('accounts', __name__)
_request_attempts = {}
_request_lock = threading.Lock()
_reset_attempts = {}
_reset_submit_attempts = {}
_reset_lock = threading.Lock()


def _utc_now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _no_store(response):
    response.headers['Cache-Control'] = 'no-store, private'
    response.headers['Referrer-Policy'] = 'no-referrer'
    return response


def _create_form_user():
    username = request.form.get('username', '').strip()
    email = request.form.get('email', '').strip()
    password = request.form.get('password', '')
    confirmation = request.form.get('password_confirm', '')
    if not email:
        raise ValueError('请填写邮箱地址。')
    if password != confirmation:
        raise ValueError('两次密码不一致。')
    try:
        return User.create(username, password, email=email)
    except sqlite3.IntegrityError:
        raise ValueError('邮箱或用户名已被使用，请换一个。') from None


@accounts_bp.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        if not _same_origin_write():
            abort(403)
        if request.content_length and request.content_length > 16_384:
            abort(413)
        key = request.remote_addr or ''
        now = time.monotonic()
        with _request_lock:
            recent = [stamp for stamp in _request_attempts.get(key, ()) if now - stamp < 900]
            if len(recent) >= 5:
                return render_template('register.html', error='注册次数过多，请 15 分钟后重试。'), 429
            _request_attempts[key] = recent + [now]
        try:
            _create_form_user()
        except ValueError as exc:
            return render_template('register.html', error=str(exc)), 400
        return render_template('register.html', registered=True)
    return render_template('register.html')


@accounts_bp.route('/admin/users', methods=['GET'])
def users():
    actor = current_user()
    if not actor or not actor.is_admin:
        abort(403)
    conn = get_conn()
    try:
        rows = conn.execute(
            'SELECT id, username, email, is_admin, is_active, created_at FROM users '
            'ORDER BY created_at DESC, username'
        ).fetchall()
        pending_reset_count = conn.execute(
            'SELECT COUNT(*) FROM password_reset_requests '
            'WHERE token_hash IS NULL AND revoked_at IS NULL'
        ).fetchone()[0]
    finally:
        conn.close()
    return render_template('admin_users.html', users=rows,
                           pending_reset_count=pending_reset_count)


@accounts_bp.route('/admin/users/create', methods=['POST'])
def create_user():
    actor = current_user()
    if not actor or not actor.is_admin:
        abort(403)
    try:
        user = _create_form_user()
    except ValueError as exc:
        flash(str(exc), 'danger')
    else:
        flash(f'已创建账号 {user.username}。请让用户登录后修改初始密码。', 'success')
    return redirect(url_for('accounts.users'))


@accounts_bp.route('/admin/users/<user_id>/<action>', methods=['POST'])
def set_user_active(user_id, action):
    actor = current_user()
    if not actor or not actor.is_admin:
        abort(403)
    if action not in {'disable', 'enable'}:
        abort(404)
    target = User.get_by_id(user_id)
    if not target:
        abort(404)
    if target.is_admin or target.id == actor.id:
        abort(403)
    active = int(action == 'enable')
    if target.is_active == bool(active):
        abort(409)
    conn = get_conn()
    try:
        updated = conn.execute(
            'UPDATE users SET is_active=?, session_version=session_version+1 WHERE id=? AND is_admin=0',
            (active, user_id),
        ).rowcount
        conn.commit()
    finally:
        conn.close()
    if not updated:
        abort(409)
    flash(f'已{"启用" if active else "停用"}账号：{target.username}', 'success')
    return redirect(url_for('accounts.users'))


@accounts_bp.route('/admin/users/<user_id>/reset-password', methods=['GET', 'POST'])
def reset_user_password(user_id):
    actor = current_user()
    if not actor or not actor.is_admin:
        abort(403)
    target = User.get_by_id(user_id)
    if not target:
        abort(404)
    if target.is_admin:
        abort(403)
    if request.method == 'POST':
        password = request.form.get('password', '')
        confirm = request.form.get('password_confirm', '')
        if len(password) < 12 or len(password) > 1024:
            return render_template('reset_user_password.html', target=target,
                                   error='新密码需为 12 至 1024 字符。'), 400
        if password != confirm:
            return render_template('reset_user_password.html', target=target,
                                   error='两次密码不一致。'), 400
        conn = get_conn()
        try:
            conn.execute('UPDATE users SET password_hash=?, session_version=session_version+1 WHERE id=?',
                         (generate_password_hash(password), target.id))
            conn.commit()
        finally:
            conn.close()
        flash(f'已重置 {target.username} 的密码，原会话已失效。请通过实验室内部渠道告知新密码。', 'success')
        return redirect(url_for('accounts.users'))
    return render_template('reset_user_password.html', target=target)


@accounts_bp.route('/account/profile', methods=['GET', 'POST'])
def profile():
    user = current_user()
    if not user:
        abort(401)
    if request.method == 'POST':
        if not user.check_password(request.form.get('current_password', '')):
            return render_template('account_profile.html', user=user,
                                   error='当前密码不正确。'), 400
        try:
            email = User.normalize_email(request.form.get('email', ''))
        except ValueError as exc:
            return render_template('account_profile.html', user=user, error=str(exc)), 400
        conn = get_conn()
        try:
            try:
                conn.execute('UPDATE users SET email=? WHERE id=?', (email, user.id))
                conn.commit()
            except sqlite3.IntegrityError:
                return render_template('account_profile.html', user=user,
                                       error='邮箱已被其他账号使用。'), 400
        finally:
            conn.close()
        flash('邮箱已更新。现在可以用用户名或邮箱登录。', 'success')
        return redirect(url_for('accounts.profile'))
    return render_template('account_profile.html', user=user)


@accounts_bp.route('/account/password', methods=['GET', 'POST'])
def change_password():
    user = current_user()
    if not user:
        abort(401)
    if request.method == 'POST':
        current = request.form.get('current_password', '')
        new = request.form.get('new_password', '')
        confirm = request.form.get('new_password_confirm', '')
        if not user.check_password(current):
            return render_template('change_password.html', error='当前密码不正确。'), 400
        if len(new) < 12 or len(new) > 1024:
            return render_template('change_password.html', error='新密码需为 12 至 1024 字符。'), 400
        if new != confirm:
            return render_template('change_password.html', error='两次新密码不一致。'), 400
        if new == current:
            return render_template('change_password.html', error='新密码不能与当前密码相同。'), 400
        conn = get_conn()
        try:
            conn.execute('UPDATE users SET password_hash=?, session_version=session_version+1 WHERE id=?',
                         (generate_password_hash(new), user.id))
            conn.commit()
        finally:
            conn.close()
        session.clear()
        return redirect(url_for('platform_login', changed='1'))
    return render_template('change_password.html')


@accounts_bp.route('/forgot-password', methods=['GET', 'POST'])
def forgot_password():
    """Accept a request without revealing whether the identifier exists."""
    if request.method == 'POST':
        if not _same_origin_write():
            abort(403)
        if request.content_length and request.content_length > 4096:
            abort(413)
        identifier = request.form.get('identifier', '').strip()
        key = (request.remote_addr or '', identifier.casefold())
        now = time.monotonic()
        with _reset_lock:
            recent = [stamp for stamp in _reset_attempts.get(key, ()) if now - stamp < 900]
            if len(recent) >= 5:
                return _no_store(make_response(render_template(
                    'forgot_password.html', error='提交次数过多，请 15 分钟后重试。'), 429))
            _reset_attempts[key] = recent + [now]
        user = User.get_by_login(identifier) if 0 < len(identifier) <= 254 else None
        if user and user.is_active and not user.is_admin:
            conn = get_conn()
            try:
                conn.execute('BEGIN IMMEDIATE')
                pending = conn.execute(
                    'SELECT 1 FROM password_reset_requests WHERE user_id=? '
                    'AND token_hash IS NULL AND revoked_at IS NULL LIMIT 1',
                    (user.id,),
                ).fetchone()
                if not pending:
                    conn.execute(
                        'INSERT INTO password_reset_requests (id, user_id, requested_at) VALUES (?, ?, ?)',
                        (secrets.token_hex(12), user.id, _utc_now()),
                    )
                conn.commit()
            finally:
                conn.close()
        return _no_store(make_response(render_template('forgot_password.html', submitted=True)))
    return _no_store(make_response(render_template('forgot_password.html')))


@accounts_bp.route('/admin/password-reset-requests')
def reset_requests():
    actor = current_user()
    if not actor or not actor.is_admin:
        abort(403)
    conn = get_conn()
    try:
        rows = conn.execute(
            'SELECT r.id, r.requested_at, r.generated_at, r.expires_at, r.used_at, '
            'r.revoked_at, r.token_hash, u.username, u.email, u.is_active, u.is_admin, '
            'reviewer.username AS generated_by_name '
            'FROM password_reset_requests r JOIN users u ON u.id=r.user_id '
            'LEFT JOIN users reviewer ON reviewer.id=r.generated_by '
            'ORDER BY r.requested_at DESC LIMIT 100'
        ).fetchall()
    finally:
        conn.close()
    now = _utc_now()
    items = []
    for row in rows:
        item = dict(row)
        if item['used_at']:
            item['status'] = 'used'
        elif item['revoked_at']:
            item['status'] = 'revoked'
        elif not item['token_hash']:
            item['status'] = 'pending'
        elif item['expires_at'] <= now:
            item['status'] = 'expired'
        else:
            item['status'] = 'generated'
        item.pop('token_hash')
        items.append(item)
    return _no_store(make_response(render_template('admin_reset_requests.html', requests=items)))


@accounts_bp.route('/admin/password-reset-requests/<request_id>/generate', methods=['POST'])
def generate_reset_link(request_id):
    actor = current_user()
    if not actor or not actor.is_admin:
        abort(403)
    token = secrets.token_urlsafe(32)
    token_hash = hashlib.sha256(token.encode('ascii')).hexdigest()
    now = datetime.now(timezone.utc)
    generated_at = now.isoformat(timespec='seconds')
    expires_at = (now + timedelta(minutes=30)).isoformat(timespec='seconds')
    conn = get_conn()
    try:
        conn.execute('BEGIN IMMEDIATE')
        row = conn.execute(
            'SELECT r.user_id, r.token_hash, r.revoked_at, u.username, u.is_active, u.is_admin '
            'FROM password_reset_requests r JOIN users u ON u.id=r.user_id WHERE r.id=?',
            (request_id,),
        ).fetchone()
        if (not row or row['token_hash'] or row['revoked_at']
                or not row['is_active'] or row['is_admin']):
            abort(409)
        conn.execute(
            'UPDATE password_reset_requests SET revoked_at=? '
            'WHERE user_id=? AND id<>? AND used_at IS NULL AND revoked_at IS NULL',
            (generated_at, row['user_id'], request_id),
        )
        conn.execute(
            'UPDATE password_reset_requests SET token_hash=?, generated_at=?, expires_at=?, generated_by=? '
            'WHERE id=?',
            (token_hash, generated_at, expires_at, actor.id, request_id),
        )
        conn.commit()
    finally:
        conn.close()
    return _no_store(make_response(render_template(
        'reset_link_once.html', token=token, username=row['username'], expires_at=expires_at)))


@accounts_bp.route('/admin/password-reset-requests/<request_id>/dismiss', methods=['POST'])
def dismiss_reset_request(request_id):
    actor = current_user()
    if not actor or not actor.is_admin:
        abort(403)
    conn = get_conn()
    try:
        changed = conn.execute(
            'UPDATE password_reset_requests SET revoked_at=? '
            'WHERE id=? AND used_at IS NULL AND revoked_at IS NULL',
            (_utc_now(), request_id),
        ).rowcount
        conn.commit()
    finally:
        conn.close()
    if not changed:
        abort(409)
    flash('申请或链接已作废。', 'success')
    return redirect(url_for('accounts.reset_requests'))


@accounts_bp.route('/reset-password', methods=['GET', 'POST'])
def reset_password():
    if request.method == 'GET':
        return _no_store(make_response(render_template(
            'reset_password.html', done=request.args.get('done') == '1')))
    if not _same_origin_write():
        abort(403)
    if request.content_length and request.content_length > 4096:
        abort(413)
    key = request.remote_addr or ''
    attempt_time = time.monotonic()
    with _reset_lock:
        recent = [stamp for stamp in _reset_submit_attempts.get(key, ())
                  if attempt_time - stamp < 900]
        if len(recent) >= 20:
            return _no_store(make_response(render_template('reset_password.html', invalid=True), 429))
        _reset_submit_attempts[key] = recent + [attempt_time]
    token = request.form.get('token', '')
    password = request.form.get('password', '')
    confirmation = request.form.get('password_confirm', '')
    if len(password) < 12 or len(password) > 1024:
        return _no_store(make_response(render_template(
            'reset_password.html', token=token, error='密码需为 12 至 1024 字符。'), 400))
    if password != confirmation:
        return _no_store(make_response(render_template(
            'reset_password.html', token=token, error='两次密码不一致。'), 400))
    if not token or len(token) > 128:
        return _no_store(make_response(render_template('reset_password.html', invalid=True), 400))
    token_hash = hashlib.sha256(token.encode('utf-8')).hexdigest()
    now = _utc_now()
    conn = get_conn()
    try:
        conn.execute('BEGIN IMMEDIATE')
        row = conn.execute(
            'SELECT r.id, r.expires_at, r.used_at, r.revoked_at, r.user_id, '
            'u.is_active, u.is_admin FROM password_reset_requests r '
            'JOIN users u ON u.id=r.user_id WHERE r.token_hash=?',
            (token_hash,),
        ).fetchone()
        invalid = (not row or row['used_at'] or row['revoked_at'] or not row['expires_at']
                   or row['expires_at'] <= now or not row['is_active'] or row['is_admin'])
        if invalid:
            conn.rollback()
        else:
            conn.execute('UPDATE users SET password_hash=?, session_version=session_version+1 WHERE id=?',
                         (generate_password_hash(password), row['user_id']))
            conn.execute('UPDATE password_reset_requests SET used_at=? WHERE id=?',
                         (now, row['id']))
            conn.commit()
    finally:
        conn.close()
    if invalid:
        return _no_store(make_response(render_template('reset_password.html', invalid=True), 400))
    session.clear()
    return redirect(url_for('accounts.reset_password', done='1'))
