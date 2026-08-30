"""Tests for the optional shared access gate used by tunnel deployments."""
from app import create_app
from config import Config


def test_platform_access_gate_protects_browser_and_api(monkeypatch):
    monkeypatch.setattr(Config, 'PLATFORM_ACCESS_PASSWORD', 'lab-test-password')
    app = create_app()
    client = app.test_client()

    browser_response = client.get('/', follow_redirects=False)
    assert browser_response.status_code == 302
    assert '/login?next=/' in browser_response.headers['Location']

    api_response = client.get('/api/system/status')
    assert api_response.status_code == 401
    assert api_response.get_json()['error'] == '需要先登录访问平台'

    bad_login = client.post('/login', data={'password': 'wrong'})
    assert bad_login.status_code == 200
    assert '访问密码不正确' in bad_login.get_data(as_text=True)

    good_login = client.post(
        '/login',
        data={'password': 'lab-test-password', 'next': '/projects/new'},
        follow_redirects=False,
    )
    assert good_login.status_code == 302
    assert good_login.headers['Location'].endswith('/projects/new')
    assert client.get('/').status_code == 200

    logout = client.post('/logout', follow_redirects=False)
    assert logout.status_code == 302
    assert client.get('/', follow_redirects=False).status_code == 302


def test_platform_access_gate_rejects_external_next(monkeypatch):
    monkeypatch.setattr(Config, 'PLATFORM_ACCESS_PASSWORD', 'lab-test-password')
    app = create_app()
    client = app.test_client()

    response = client.post(
        '/login?next=https://evil.example',
        data={'password': 'lab-test-password'},
        follow_redirects=False,
    )
    assert response.status_code == 302
    assert response.headers['Location'].endswith('/')
