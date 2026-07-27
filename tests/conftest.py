# tests/conftest.py
"""pytest 共享 fixtures — 自动设置项目根目录到 sys.path 和环境变量"""
import os
import sys
import tempfile

# Test temporary files and compilation caches must not consume the root
# filesystem.  Tests intentionally override the application setting so a
# developer's production RUNTIME_TMP_DIR is never mixed with test artifacts.
_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_test_tmp_root = os.path.abspath(os.environ.get(
    'PYTEST_RUNTIME_TMP_DIR',
    os.path.join(_project_root, 'data', 'runtime_tmp', 'pytest'),
))
_test_numba_cache = os.path.abspath(os.environ.get(
    'PYTEST_NUMBA_CACHE_DIR', os.path.join(_test_tmp_root, 'numba_cache'),
))
_test_mpl_cache = os.path.abspath(os.environ.get(
    'PYTEST_MPLCONFIGDIR', os.path.join(_test_tmp_root, 'mplconfig'),
))
for _directory in (_test_tmp_root, _test_numba_cache, _test_mpl_cache):
    os.makedirs(_directory, exist_ok=True)

os.environ['RUNTIME_TMP_DIR'] = _test_tmp_root
for _variable in ('TMPDIR', 'TMP', 'TEMP'):
    os.environ[_variable] = _test_tmp_root
tempfile.tempdir = _test_tmp_root

# 环境变量：确保裸环境测试可复现
os.environ.setdefault('NUMBA_DISABLE_JIT', '1')
os.environ['NUMBA_CACHE_DIR'] = _test_numba_cache
os.environ['MPLCONFIGDIR'] = _test_mpl_cache

# 将项目根目录加入 sys.path
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
