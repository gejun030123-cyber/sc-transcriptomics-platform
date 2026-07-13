# tests/conftest.py
"""pytest 共享 fixtures — 自动设置项目根目录到 sys.path 和环境变量"""
import os
import sys

# 环境变量：确保裸环境测试可复现
os.environ.setdefault('NUMBA_DISABLE_JIT', '1')
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/numba_cache_test')

# 将项目根目录加入 sys.path
_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

import pytest
from database import init_db, get_conn


@pytest.fixture
def test_project(tmp_path, monkeypatch):
    """创建测试用 project（数据和目录隔离到 tmp_path），返回 project ID."""
    import uuid
    from config import Config

    pid = 'test_' + uuid.uuid4().hex[:8]

    # 隔离 DATA_DIR 到临时目录
    data_dir = tmp_path / 'data'
    data_dir.mkdir(exist_ok=True)

    # 隔离 DB_PATH 到临时目录
    db_path = tmp_path / 'test.db'
    monkeypatch.setattr(Config, 'DATA_DIR', str(data_dir))
    monkeypatch.setattr(Config, 'DB_PATH', str(db_path))

    # 重新初始化数据库（使用临时路径）
    init_db()

    # 创建测试项目目录
    proj_dir = Config.project_dir(pid)
    os.makedirs(proj_dir, exist_ok=True)
    os.makedirs(os.path.join(proj_dir, 'intermediate'), exist_ok=True)
    os.makedirs(os.path.join(proj_dir, 'uploads'), exist_ok=True)
    os.makedirs(os.path.join(proj_dir, 'branches'), exist_ok=True)

    conn = get_conn()
    try:
        conn.execute(
            "INSERT OR REPLACE INTO projects (id, name, description, status) "
            "VALUES (?, ?, ?, ?)",
            (pid, 'Test Project', 'Test project for unit tests', 'active')
        )
        conn.commit()
    finally:
        conn.close()

    yield pid
    # tmp_path 会在测试结束后自动清理，无需手动删除文件
