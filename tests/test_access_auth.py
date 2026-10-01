"""Individual login and cross-user project isolation."""
import os

import pytest
from app import create_app
from config import Config
from models import AnalysisTask, Project, ResultFile, User


@pytest.fixture
def accounts(tmp_path, monkeypatch):
    monkeypatch.setattr(Config, 'AUTH_REQUIRED', True)
    monkeypatch.setattr(Config, 'AI_API_TOKEN', '')
    monkeypatch.setattr(Config, 'DATA_DIR', str(tmp_path / 'data'))
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'users.db'))
    app = create_app()
    app.config['TESTING'] = True
    alice = User.create('alice', 'alice-password-123')
    bob = User.create('bob', 'bob-password-123')
    pa = Project(name='Alice private project', owner_user_id=alice.id)
    pb = Project(name='Bob private project', owner_user_id=bob.id)
    pa.save()
    pb.save()
    for project in (pa, pb):
        os.makedirs(Config.results_dir(project.id), exist_ok=True)
    return app, alice, bob, pa, pb


def login(client, username, password):
    return client.post('/login', data={'username': username, 'password': password},
                       headers={'Origin': 'http://localhost'}, follow_redirects=False)


def test_personal_login_and_project_lists(accounts):
    app, _, _, pa, pb = accounts
    alice_client = app.test_client()
    assert alice_client.get('/').status_code == 302
    assert alice_client.get('/api/system/status').status_code == 401
    assert login(alice_client, 'alice', 'wrong').status_code == 200
    assert login(alice_client, 'alice', 'alice-password-123').status_code == 302
    page = alice_client.get('/projects').get_data(as_text=True)
    assert pa.name in page
    assert pb.name not in page
    assert alice_client.get(f'/projects/{pb.id}').status_code == 404
    assert alice_client.post(f'/projects/{pb.id}/delete', headers={'Origin': 'http://localhost'}).status_code == 404
    assert alice_client.post('/logout', headers={'Origin': 'http://localhost'}).status_code == 302
    assert alice_client.get('/api/system/status').status_code == 401


def test_task_file_and_chat_ids_cannot_cross_users(accounts):
    app, _, _, pa, pb = accounts
    task = AnalysisTask(project_id=pb.id, module_name='qc', status='completed')
    task.save()
    path = os.path.join(Config.results_dir(pb.id), 'private.csv')
    with open(path, 'w') as handle:
        handle.write('gene,value\nA,1\n')
    result = ResultFile.create(task.id, pb.id, 'csv', 'table', 'private', path)
    alice_client = app.test_client()
    login(alice_client, 'alice', 'alice-password-123')
    assert alice_client.get(f'/api/tasks/{task.id}/status').status_code == 404
    assert alice_client.get(f'/api/result-file/{result.id}', query_string={'project_id': pa.id}).status_code == 404
    assert alice_client.get('/api/data-info', query_string={'file_path': path}).status_code == 403
    assert alice_client.get(f'/api/chat/history/{pb.id}').status_code == 404
    assert alice_client.post('/api/chat/approve', json={'project_id': pb.id, 'tool_name': 'get_project_status'},
                             headers={'Origin': 'http://localhost'}).status_code == 404


def test_ai_settings_and_presets_are_user_scoped(accounts):
    app, _, _, pa, pb = accounts
    alice_client = app.test_client()
    bob_client = app.test_client()
    login(alice_client, 'alice', 'alice-password-123')
    login(bob_client, 'bob', 'bob-password-123')
    payload = {'api_url': 'https://alice.example/v1', 'model': 'alice-model',
               'provider': 'openai', 'api_key': 'alice-private-key'}
    assert alice_client.post('/api/settings/ai', json=payload,
                             headers={'Origin': 'http://localhost'}).status_code == 200
    assert bob_client.get('/api/settings/ai').get_json()['api_url'] != payload['api_url']
    assert alice_client.post('/api/presets', json={'name': 'private', 'project_id': pb.id},
                             headers={'Origin': 'http://localhost'}).status_code == 404
    assert alice_client.get('/api/presets', query_string={'project_id': pb.id}).status_code == 404
    assert alice_client.post('/api/presets', json={'name': 'global', 'scope': 'global'},
                             headers={'Origin': 'http://localhost'}).status_code == 403
    assert alice_client.post('/api/presets', json={'name': 'mine', 'project_id': pa.id},
                             headers={'Origin': 'http://localhost'}).status_code == 200


def test_cross_site_write_and_external_next_are_rejected(accounts):
    app, _, _, pa, _ = accounts
    client = app.test_client()
    assert login(client, 'alice', 'alice-password-123').status_code == 302
    assert client.post(f'/projects/{pa.id}/delete', headers={'Origin': 'https://evil.example'}).status_code == 403
    response = client.post('/login?next=https://evil.example',
                           data={'username': 'alice', 'password': 'alice-password-123'},
                           headers={'Origin': 'http://localhost'})
    assert response.status_code == 302
    assert response.headers['Location'].endswith('/')


def test_uploaded_capture_bed_is_private(accounts):
    import io
    app, _, _, pa, pb = accounts
    alice_client = app.test_client()
    bob_client = app.test_client()
    login(alice_client, 'alice', 'alice-password-123')
    login(bob_client, 'bob', 'bob-password-123')
    response = alice_client.post(
        f'/api/projects/{pa.id}/wes/capture-kits/upload',
        data={'capture_kit_id': 'alice_kit', 'name': 'Alice kit', 'version': 'v1',
              'calling_bed': (io.BytesIO(b'chr1\t0\t10\n'), 'alice.bed')},
        content_type='multipart/form-data',
        headers={'Origin': 'http://localhost'},
    )
    assert response.status_code == 201, response.get_json()
    own = alice_client.get(f'/api/projects/{pa.id}/wes/capture-kits').get_json()
    other = bob_client.get(f'/api/projects/{pb.id}/wes/capture-kits').get_json()
    assert 'alice_kit' in str(own)
    assert 'alice_kit' not in str(other)


def test_claim_existing_projects_backs_up_database_and_moves_ai_settings(tmp_path, monkeypatch):
    import sys
    from database import get_conn, init_db
    from manage_users import main
    db_path = tmp_path / 'legacy.db'
    monkeypatch.setattr(Config, 'DB_PATH', str(db_path))
    init_db()
    project = Project(name='Legacy project')
    project.save()
    conn = get_conn()
    try:
        conn.execute("INSERT INTO platform_settings (key, value) VALUES ('ai_model', 'legacy-model')")
        conn.commit()
    finally:
        conn.close()
    passwords = iter(['admin-password-123', 'admin-password-123'])
    monkeypatch.setattr('getpass.getpass', lambda _: next(passwords))
    monkeypatch.setattr(sys, 'argv', ['manage_users.py', 'create', 'labadmin', '--admin', '--claim-existing'])
    main()
    admin = User.get_by_username('labadmin')
    assert Project.get_by_id(project.id).owner_user_id == admin.id
    assert (tmp_path / 'legacy.db.before-user-ownership.bak').exists()
    conn = get_conn()
    try:
        assert conn.execute("SELECT value FROM user_ai_settings WHERE user_id=? AND key='ai_model'", (admin.id,)).fetchone()['value'] == 'legacy-model'
        assert conn.execute("SELECT COUNT(*) AS n FROM platform_settings").fetchone()['n'] == 0
    finally:
        conn.close()


def test_server_source_grants_limit_each_user(accounts, tmp_path, monkeypatch):
    import json
    from flask import g
    from database import get_conn
    from modules.workflows.assets import register_data_asset
    app, alice, _, pa, _ = accounts
    sc_root = tmp_path / 'sc-sources'
    wes_root = tmp_path / 'wes-sources'
    for root in (sc_root, wes_root):
        (root / 'alice').mkdir(parents=True)
        (root / 'bob').mkdir()
    own_wes = wes_root / 'alice' / 'sample.fastq.gz'
    other_wes = wes_root / 'bob' / 'sample.fastq.gz'
    own_wes.write_bytes(b'own')
    other_wes.write_bytes(b'other')
    monkeypatch.setenv('SC_BATCH_SOURCE_ROOTS', str(sc_root))
    monkeypatch.setenv('WES_SOURCE_ROOTS', str(wes_root))
    conn = get_conn()
    try:
        conn.execute('UPDATE users SET sc_source_roots_json=?, wes_source_roots_json=? WHERE id=?',
                     (json.dumps([str(sc_root / 'alice')]), json.dumps([str(wes_root / 'alice')]), alice.id))
        conn.commit()
    finally:
        conn.close()
    with app.test_request_context(f'/api/projects/{pa.id}/wes/preflight'):
        g.current_user = User.get_by_id(alice.id)
        assert Config.validate_sc_batch_source_path(str(sc_root / 'alice'))
        assert Config.validate_wes_source_path(str(own_wes))
        with pytest.raises(ValueError):
            Config.validate_sc_batch_source_path(str(sc_root / 'bob'))
        with pytest.raises(ValueError):
            Config.validate_wes_source_path(str(other_wes))
        register_data_asset(pa.id, artifact_kind='fastq', file_path=str(own_wes))
        with pytest.raises(ValueError):
            register_data_asset(pa.id, artifact_kind='fastq', file_path=str(other_wes))


def test_disabled_account_invalidates_existing_session(accounts):
    from database import get_conn
    app, alice, _, pa, _ = accounts
    client = app.test_client()
    login(client, 'alice', 'alice-password-123')
    assert client.get(f'/projects/{pa.id}').status_code == 200
    conn = get_conn()
    try:
        conn.execute('UPDATE users SET is_active=0, session_version=session_version+1 WHERE id=?', (alice.id,))
        conn.commit()
    finally:
        conn.close()
    assert client.get(f'/projects/{pa.id}', follow_redirects=False).status_code == 302
    assert client.get('/api/system/status').status_code == 401


def test_existing_admin_can_claim_unowned_legacy_projects(tmp_path, monkeypatch):
    import sys
    from database import init_db
    from manage_users import main
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'later.db'))
    init_db()
    admin = User.create('lateradmin', 'admin-password-123', is_admin=True)
    legacy = Project(name='Unowned')
    legacy.save()
    monkeypatch.setattr(sys, 'argv', ['manage_users.py', 'claim-existing', admin.username])
    main()
    assert Project.get_by_id(legacy.id).owner_user_id == admin.id


def test_source_grant_cli_updates_account_and_invalidates_session(accounts, tmp_path, monkeypatch):
    import json
    import sys
    from manage_users import main
    app, alice, _, pa, _ = accounts
    root = tmp_path / 'approved'
    own = root / 'alice'
    own.mkdir(parents=True)
    monkeypatch.setenv('SC_BATCH_SOURCE_ROOTS', str(root))
    client = app.test_client()
    login(client, 'alice', 'alice-password-123')
    monkeypatch.setattr(sys, 'argv', ['manage_users.py', 'grant-source-root', 'alice', 'sc', str(own)])
    main()
    assert str(own) in json.loads(User.get_by_id(alice.id).sc_source_roots_json)
    assert client.get(f'/projects/{pa.id}', follow_redirects=False).status_code == 302
    monkeypatch.setattr(sys, 'argv', ['manage_users.py', 'revoke-source-root', 'alice', 'sc', str(own)])
    main()
    assert json.loads(User.get_by_id(alice.id).sc_source_roots_json) == []
