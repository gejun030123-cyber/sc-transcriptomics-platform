"""Administrator-mediated, one-use password recovery."""
import hashlib
import re

import pytest

from app import create_app
from config import Config
from database import get_conn
from models import User
from routes.accounts import _reset_attempts, _reset_submit_attempts

ORIGIN = {'Origin': 'http://localhost'}


@pytest.fixture
def recovery_app(tmp_path, monkeypatch):
    monkeypatch.setattr(Config, 'AUTH_REQUIRED', True)
    monkeypatch.setattr(Config, 'AI_API_TOKEN', '')
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'recovery.db'))
    monkeypatch.setattr(Config, 'DATA_DIR', str(tmp_path / 'data'))
    _reset_attempts.clear()
    _reset_submit_attempts.clear()
    app = create_app()
    app.config['TESTING'] = True
    admin = User.create('recoveradmin', 'admin-password-123', is_admin=True)
    user = User.create('recoveruser', 'old-password-123', email='User@Lab.Example')
    return app, admin, user


def login(client, identifier, password):
    return client.post('/login', data={'username': identifier, 'password': password},
                       headers=ORIGIN, follow_redirects=False)


def request_reset(client, identifier):
    return client.post('/forgot-password', data={'identifier': identifier}, headers=ORIGIN)


def pending_id():
    conn = get_conn()
    try:
        row = conn.execute('SELECT id FROM password_reset_requests WHERE token_hash IS NULL '
                           'AND revoked_at IS NULL ORDER BY requested_at DESC LIMIT 1').fetchone()
        return row['id'] if row else None
    finally:
        conn.close()


def generate(client, request_id):
    response = client.post(f'/admin/password-reset-requests/{request_id}/generate', headers=ORIGIN)
    assert response.status_code == 200
    html = response.get_data(as_text=True)
    match = re.search(r'/reset-password#([A-Za-z0-9_-]+)', html)
    assert match
    return match.group(1), html, response


def test_recovery_request_is_generic_and_admin_only(recovery_app):
    app, admin, user = recovery_app
    applicant = app.test_client()
    reviewer = app.test_client()
    assert applicant.get('/forgot-password').status_code == 200
    assert applicant.get('/reset-password').status_code == 200
    assert applicant.post('/forgot-password', data={'identifier': user.username},
                          headers={'Origin': 'https://other.example'}).status_code == 403
    unknown = request_reset(applicant, 'missing@lab.example')
    assert unknown.status_code == 200
    assert pending_id() is None
    wrong_case = request_reset(applicant, 'user@lab.example')
    assert wrong_case.status_code == 200
    assert pending_id() is None
    known = request_reset(applicant, 'User@Lab.Example')
    assert known.status_code == 200
    assert '申请已提交' in known.get_data(as_text=True)
    request_id = pending_id()
    assert request_id
    assert request_reset(applicant, user.username).status_code == 200
    conn = get_conn()
    try:
        assert conn.execute('SELECT COUNT(*) FROM password_reset_requests').fetchone()[0] == 1
    finally:
        conn.close()
    assert applicant.get('/admin/password-reset-requests').status_code == 302
    assert login(applicant, user.username, 'old-password-123').status_code == 302
    assert applicant.get('/admin/password-reset-requests').status_code == 403
    assert login(reviewer, admin.username, 'admin-password-123').status_code == 302
    assert user.username in reviewer.get('/admin/password-reset-requests').get_data(as_text=True)
    assert reviewer.post(f'/admin/password-reset-requests/{request_id}/generate',
                         headers={'Origin': 'https://other.example'}).status_code == 403


def test_one_use_link_resets_password_and_invalidates_sessions(recovery_app):
    app, admin, user = recovery_app
    applicant = app.test_client()
    other_session = app.test_client()
    reviewer = app.test_client()
    assert login(applicant, user.username, 'old-password-123').status_code == 302
    assert login(other_session, user.username, 'old-password-123').status_code == 302
    assert request_reset(applicant, user.username).status_code == 200
    request_id = pending_id()
    login(reviewer, admin.username, 'admin-password-123')
    token, html, response = generate(reviewer, request_id)
    assert '离开本页后无法再次查看' in html
    assert response.headers['Cache-Control'] == 'no-store, private'
    assert response.headers['Referrer-Policy'] == 'no-referrer'
    assert reviewer.post(f'/admin/password-reset-requests/{request_id}/generate',
                         headers=ORIGIN).status_code == 409
    assert token not in reviewer.get('/admin/password-reset-requests').get_data(as_text=True)
    conn = get_conn()
    try:
        row = conn.execute('SELECT token_hash FROM password_reset_requests WHERE id=?',
                           (request_id,)).fetchone()
        assert row['token_hash'] == hashlib.sha256(token.encode('ascii')).hexdigest()
    finally:
        conn.close()
    assert applicant.post('/reset-password', data={
        'token': token, 'password': 'new-password-123', 'password_confirm': 'new-password-123',
    }, headers={'Origin': 'https://other.example'}).status_code == 403
    response = applicant.post('/reset-password', data={
        'token': token, 'password': 'new-password-123', 'password_confirm': 'new-password-123',
    }, headers=ORIGIN)
    assert response.status_code == 302
    assert other_session.get('/').status_code == 302
    assert login(applicant, user.username, 'old-password-123').status_code == 200
    assert login(applicant, user.username, 'new-password-123').status_code == 302
    assert applicant.post('/reset-password', data={
        'token': token, 'password': 'another-password-123',
        'password_confirm': 'another-password-123',
    }, headers=ORIGIN).status_code == 400
    page = reviewer.get('/admin/password-reset-requests').get_data(as_text=True)
    assert '已使用' in page
    assert token not in page


def test_new_link_revokes_old_link_and_expiry_is_enforced(recovery_app):
    app, admin, user = recovery_app
    applicant = app.test_client()
    reviewer = app.test_client()
    login(reviewer, admin.username, 'admin-password-123')
    request_reset(applicant, user.username)
    first_id = pending_id()
    first_token, _, _ = generate(reviewer, first_id)
    request_reset(applicant, user.username)
    second_id = pending_id()
    assert second_id and second_id != first_id
    second_token, _, _ = generate(reviewer, second_id)
    conn = get_conn()
    try:
        assert conn.execute('SELECT revoked_at FROM password_reset_requests WHERE id=?',
                            (first_id,)).fetchone()['revoked_at']
        conn.execute("UPDATE password_reset_requests SET expires_at='2000-01-01T00:00:00+00:00' WHERE id=?",
                     (second_id,))
        conn.commit()
    finally:
        conn.close()
    for token in (first_token, second_token):
        assert applicant.post('/reset-password', data={
            'token': token, 'password': 'new-password-123',
            'password_confirm': 'new-password-123',
        }, headers=ORIGIN).status_code == 400
    assert login(applicant, user.username, 'old-password-123').status_code == 302
    assert '已过期' in reviewer.get('/admin/password-reset-requests').get_data(as_text=True)


def test_admin_can_dismiss_and_admin_account_uses_existing_cli(recovery_app):
    app, admin, user = recovery_app
    applicant = app.test_client()
    reviewer = app.test_client()
    assert request_reset(applicant, admin.username).status_code == 200
    assert pending_id() is None
    request_reset(applicant, user.username)
    request_id = pending_id()
    login(reviewer, admin.username, 'admin-password-123')
    assert reviewer.post(f'/admin/password-reset-requests/{request_id}/dismiss',
                         headers=ORIGIN).status_code == 302
    assert reviewer.post(f'/admin/password-reset-requests/{request_id}/generate',
                         headers=ORIGIN).status_code == 409
    assert '已作废' in reviewer.get('/admin/password-reset-requests').get_data(as_text=True)


def test_anonymous_recovery_rate_limits(recovery_app):
    app, _, _ = recovery_app
    client = app.test_client()
    for _ in range(5):
        assert request_reset(client, 'missing@lab.example').status_code == 200
    assert request_reset(client, 'missing@lab.example').status_code == 429
    # A shared lab proxy should not block another account's request.
    assert request_reset(client, 'other@lab.example').status_code == 200
    payload = {'token': 'invalid-token', 'password': 'new-password-123',
               'password_confirm': 'new-password-123'}
    for _ in range(20):
        assert client.post('/reset-password', data=payload, headers=ORIGIN).status_code == 400
    assert client.post('/reset-password', data=payload, headers=ORIGIN).status_code == 429
