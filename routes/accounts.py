"""Self registration and small-lab account management pages."""
import sqlite3
import threading
import time

from flask import Blueprint, abort, flash, redirect, render_template, request, session, url_for
from werkzeug.security import generate_password_hash

from database import get_conn
from models import User
from routes.auth import _same_origin_write, current_user

accounts_bp = Blueprint('accounts', __name__)
_request_attempts = {}
_request_lock = threading.Lock()


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
    finally:
        conn.close()
    return render_template('admin_users.html', users=rows)


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
