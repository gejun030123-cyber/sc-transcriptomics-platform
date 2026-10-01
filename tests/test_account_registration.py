"""Web registration, email login, and account management."""
import sqlite3

import pytest
from app import create_app
from config import Config
from database import get_conn, init_db
from models import User
from routes.accounts import _request_attempts

ORIGIN = {'Origin': 'http://localhost'}


@pytest.fixture
def account_app(tmp_path, monkeypatch):
    monkeypatch.setattr(Config, 'AUTH_REQUIRED', True)
    monkeypatch.setattr(Config, 'AI_API_TOKEN', '')
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'accounts.db'))
    monkeypatch.setattr(Config, 'DATA_DIR', str(tmp_path / 'data'))
    _request_attempts.clear()
    app = create_app()
    app.config['TESTING'] = True
    admin = User.create('reviewadmin', 'admin-password-123', is_admin=True)
    return app, admin


def register(client, username='newresearcher', **overrides):
    data = {'username': username, 'email': f'{username}@lab.example',
            'password': 'research-password-123',
            'password_confirm': 'research-password-123'}
    data.update(overrides)
    return client.post('/register', data=data, headers=ORIGIN)


def login(client, identifier, password):
    return client.post('/login', data={'username': identifier, 'password': password},
                       headers=ORIGIN, follow_redirects=False)


def test_registration_logs_in_by_username_or_email_without_approval(account_app):
    app, _ = account_app
    client = app.test_client()
    assert client.get('/register').status_code == 200
    assert register(client, email='NewResearcher@Lab.Example').status_code == 200
    user = User.get_by_username('newresearcher')
    assert user.is_active and not user.is_admin
    assert user.email == 'NewResearcher@Lab.Example'
    assert login(client, 'newresearcher', 'research-password-123').status_code == 302
    assert client.get('/admin/users').status_code == 403
    assert client.post('/logout', headers=ORIGIN).status_code == 302
    assert login(client, 'newresearcher@lab.example', 'research-password-123').status_code == 200
    assert login(client, 'NewResearcher@Lab.Example', 'research-password-123').status_code == 302
    assert client.get('/').status_code == 200


def test_admin_can_create_disable_enable_and_reset(account_app):
    app, admin = account_app
    applicant = app.test_client()
    reviewer = app.test_client()
    assert register(applicant).status_code == 200
    user = User.get_by_username('newresearcher')
    assert applicant.get('/admin/users').status_code == 302
    assert login(reviewer, admin.username, 'admin-password-123').status_code == 302
    page = reviewer.get('/admin/users').get_data(as_text=True)
    assert user.email in page
    assert reviewer.post(f'/admin/users/{user.id}/disable', headers=ORIGIN).status_code == 302
    assert not User.get_by_id(user.id).is_active
    assert login(applicant, user.email, 'research-password-123').status_code == 200
    assert reviewer.post(f'/admin/users/{user.id}/enable', headers=ORIGIN).status_code == 302
    assert login(applicant, user.email, 'research-password-123').status_code == 302
    assert reviewer.get(f'/admin/users/{user.id}/reset-password').status_code == 200
    assert reviewer.post(f'/admin/users/{user.id}/reset-password', data={
        'password': 'replacement-password-123', 'password_confirm': 'replacement-password-123',
    }, headers=ORIGIN).status_code == 302
    assert applicant.get('/').status_code == 302
    assert login(applicant, user.email, 'replacement-password-123').status_code == 302
    assert reviewer.post('/admin/users/create', data={
        'username': 'createduser', 'email': 'created@lab.example',
        'password': 'temporary-password-123', 'password_confirm': 'temporary-password-123',
    }, headers=ORIGIN).status_code == 302
    assert login(app.test_client(), 'created@lab.example', 'temporary-password-123').status_code == 302
    assert reviewer.post(f'/admin/users/{admin.id}/disable', headers=ORIGIN).status_code == 403


def test_registration_validates_email_origin_and_rate_limit(account_app):
    app, _ = account_app
    client = app.test_client()
    assert client.post('/register', data={'username': 'bad'},
                       headers={'Origin': 'https://other.example'}).status_code == 403
    assert register(client, 'missing', email='').status_code == 400
    assert register(client, 'invalid', email='not-an-email').status_code == 400
    assert register(client, 'sameuser', email='one@lab.example').status_code == 200
    assert register(client, 'second', email='ONE@LAB.EXAMPLE').status_code == 200
    assert User.get_by_username('second').email == 'ONE@LAB.EXAMPLE'
    assert register(client, 'duplicate', email='one@lab.example').status_code == 400
    assert User.get_by_username('duplicate') is None
    assert register(client, 'extra').status_code == 429
    assert User.get_by_username('extra') is None


def test_profile_email_update_and_password_change(account_app):
    app, admin = account_app
    client = app.test_client()
    other = app.test_client()
    login(client, admin.username, 'admin-password-123')
    login(other, admin.username, 'admin-password-123')
    assert client.post('/account/profile', data={
        'email': 'admin@lab.example', 'current_password': 'wrong',
    }, headers=ORIGIN).status_code == 400
    assert client.post('/account/profile', data={
        'email': 'ADMIN@LAB.EXAMPLE', 'current_password': 'admin-password-123',
    }, headers=ORIGIN).status_code == 302
    assert User.get_by_id(admin.id).email == 'ADMIN@LAB.EXAMPLE'
    second = User.create('anotheruser', 'another-password-123', email='another@lab.example')
    assert client.post('/account/profile', data={
        'email': 'another@lab.example', 'current_password': 'admin-password-123',
    }, headers=ORIGIN).status_code == 400
    assert client.post('/account/profile', data={
        'email': 'ANOTHER@LAB.EXAMPLE', 'current_password': 'admin-password-123',
    }, headers=ORIGIN).status_code == 302
    assert User.get_by_id(admin.id).email == 'ANOTHER@LAB.EXAMPLE'
    assert User.get_by_id(second.id).email == 'another@lab.example'
    assert client.post('/account/password', data={
        'current_password': 'admin-password-123', 'new_password': 'new-admin-password-123',
        'new_password_confirm': 'new-admin-password-123',
    }, headers=ORIGIN).status_code == 302
    assert client.get('/admin/users').status_code == 302
    assert other.get('/admin/users').status_code == 302
    assert login(client, 'another@lab.example', 'new-admin-password-123').status_code == 200
    assert login(client, 'ANOTHER@LAB.EXAMPLE', 'new-admin-password-123').status_code == 302


def test_existing_users_migrate_with_optional_email(tmp_path, monkeypatch):
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'legacy.db'))
    conn = sqlite3.connect(Config.DB_PATH)
    conn.execute('CREATE TABLE users (id TEXT PRIMARY KEY, username TEXT NOT NULL UNIQUE, '
                 'password_hash TEXT NOT NULL, is_admin INTEGER NOT NULL DEFAULT 0, '
                 'is_active INTEGER NOT NULL DEFAULT 1, session_version INTEGER NOT NULL DEFAULT 0, '
                 "sc_source_roots_json TEXT NOT NULL DEFAULT '[]', "
                 "wes_source_roots_json TEXT NOT NULL DEFAULT '[]', "
                 'created_at DATETIME DEFAULT CURRENT_TIMESTAMP)')
    conn.execute("INSERT INTO users(id, username, password_hash) VALUES ('old', 'olduser', 'hash')")
    conn.commit()
    conn.close()
    init_db()
    migrated = User.get_by_id('old')
    assert migrated.email is None
    assert migrated.is_active
    db = get_conn()
    try:
        assert db.execute("SELECT email FROM users WHERE id='old'").fetchone()[0] is None
    finally:
        db.close()
