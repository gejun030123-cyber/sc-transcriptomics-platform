# tests/conftest.py
"""pytest 共享 fixtures — 自动设置项目根目录到 sys.path 和环境变量"""
import os
import sys
import tempfile

# Test temporary files and compilation caches must not consume the root
# filesystem.  Tests intentionally override the application setting so a
# developer's production RUNTIME_TMP_DIR is never mixed with test artifacts.
_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 默认测试临时目录必须位于仓库内可写、且不经过符号链接的位置：
# data 是指向实验室数据卷的符号链接（只读挂载时 pytest 解析真实路径后会
# 因只读而失败），因此测试产物默认放在仓库根目录的 .test_tmp/ 下。
_test_tmp_root = os.path.abspath(os.environ.get(
    'PYTEST_RUNTIME_TMP_DIR',
    os.path.join(_project_root, '.test_tmp', 'pytest'),
))
_test_numba_cache = os.path.abspath(os.environ.get(
    'PYTEST_NUMBA_CACHE_DIR', os.path.join(_test_tmp_root, 'numba_cache'),
))
_test_mpl_cache = os.path.abspath(os.environ.get(
    'PYTEST_MPLCONFIGDIR', os.path.join(_test_tmp_root, 'mplconfig'),
))
# Tests that do not use the project fixture still open the platform database.
# Keep that default database on the test volume; the production .env continues
# to point at the lab data disk.  This also makes pytest work in a restricted
# CI/sandbox where the lab mount is readable but not writable.
_test_data_root = os.path.abspath(os.environ.get(
    'PYTEST_PLATFORM_DATA_DIR', os.path.join(_test_tmp_root, 'platform-data'),
))
_test_db_path = os.path.abspath(os.environ.get(
    'PYTEST_PLATFORM_DB_PATH', os.path.join(_test_tmp_root, 'platform-test.db'),
))
os.environ.setdefault('DATA_DIR', _test_data_root)
os.environ.setdefault('DB_PATH', _test_db_path)
for _directory in (_test_tmp_root, _test_numba_cache, _test_mpl_cache, _test_data_root):
    os.makedirs(_directory, exist_ok=True)

os.environ['RUNTIME_TMP_DIR'] = _test_tmp_root
for _variable in ('TMPDIR', 'TMP', 'TEMP'):
    os.environ[_variable] = _test_tmp_root
tempfile.tempdir = _test_tmp_root

# 环境变量：确保裸环境测试可复现
os.environ.setdefault('NUMBA_DISABLE_JIT', '1')
# Test clients must not inherit a developer's deployment access gate or an
# enabled WES executor from the shell/.env.  Individual tests explicitly
# monkeypatch these values when they exercise either gate.
os.environ['PLATFORM_ACCESS_PASSWORD'] = ''
os.environ['WES_EXECUTOR_ENABLED'] = 'false'
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
