import os
import re
import tempfile


def _load_dotenv(path=None):
    """Load simple KEY=VALUE pairs from .env without overriding real env vars."""
    env_path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), '.env')
    if not os.path.exists(env_path):
        return
    try:
        with open(env_path, encoding='utf-8') as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line or line.startswith('#') or '=' not in line:
                    continue
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                if not key or key in os.environ:
                    continue
                if (value.startswith('"') and value.endswith('"')) or (
                    value.startswith("'") and value.endswith("'")
                ):
                    value = value[1:-1]
                os.environ[key] = value
    except OSError:
        return


_load_dotenv()


class Config:
    _BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    @staticmethod
    def _load_secret():
        key = os.environ.get('SECRET_KEY')
        if key:
            return key
        secret_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'instance', '.secret_key')
        try:
            with open(secret_file) as f:
                return f.read().strip()
        except FileNotFoundError:
            import secrets
            key = secrets.token_hex(24)
            os.makedirs(os.path.dirname(secret_file), exist_ok=True)
            with open(secret_file, 'w') as f:
                f.write(key)
            os.chmod(secret_file, 0o600)
            return key

    SECRET_KEY = _load_secret.__func__()
    DATA_DIR = os.environ.get('DATA_DIR', os.path.join(_BASE_DIR, 'data'))
    # Keep transient analysis files on the project data volume by default.
    # ``/tmp`` is commonly mounted on the root filesystem and can fill up when
    # 10x compatibility staging or plotting jobs handle large matrices.
    RUNTIME_TMP_DIR = os.path.abspath(os.environ.get(
        'RUNTIME_TMP_DIR', os.path.join(DATA_DIR, 'runtime_tmp')
    ))
    NUMBA_CACHE_DIR = os.path.abspath(os.environ.get(
        'NUMBA_CACHE_DIR', os.path.join(RUNTIME_TMP_DIR, 'numba_cache')
    ))
    MPLCONFIG_DIR = os.path.abspath(os.environ.get(
        'MPLCONFIGDIR', os.path.join(RUNTIME_TMP_DIR, 'mplconfig')
    ))
    CACHE_DIR = os.path.abspath(os.environ.get(
        'CACHE_DIR', os.path.join(DATA_DIR, 'cache')
    ))
    DB_PATH = os.environ.get('DB_PATH', os.path.join(_BASE_DIR, 'instance', 'bioinfo.db'))
    CELLMARKER_PATH = os.environ.get('CELLMARKER_PATH', os.path.join(os.path.dirname(_BASE_DIR), 'CellMarker_Augmented_2021.txt'))

    MAX_WORKERS = int(os.environ.get('MAX_WORKERS', '2'))
    CHUNK_SIZE_MB = 50
    PLOTLY_MAX_CELLS = 50000
    MIN_FREE_RAM_GB = float(os.environ.get('MIN_FREE_RAM_GB', '4'))
    CUDA_DEVICES = os.environ.get('CUDA_DEVICES', '0,1')

    # Server-resident 10x data are intentionally opt-in.  A web request may
    # only refer to a directory beneath one of these administrator configured
    # roots; arbitrary absolute paths and symlinks remain disallowed.  Example:
    # SC_BATCH_SOURCE_ROOTS=/home/oelab/data/GJ:/mnt/sc_data
    SC_BATCH_SOURCE_ROOTS = os.environ.get('SC_BATCH_SOURCE_ROOTS', '')

    # AI 对话配置（支持 OpenAI compatible 和 Anthropic messages compatible API）
    AI_API_KEY = os.environ.get('AI_API_KEY', '')
    AI_API_URL = os.environ.get('AI_API_URL', 'https://token-plan-cn.xiaomimimo.com/anthropic')
    AI_MODEL = os.environ.get('AI_MODEL', 'mimo-v2.5-pro')
    AI_API_TOKEN = os.environ.get('AI_API_TOKEN', '')

    _PID_RE = re.compile(r'^[a-zA-Z0-9_-]{1,64}$')

    @classmethod
    def _validate_pid(cls, pid):
        if not cls._PID_RE.match(pid):
            raise ValueError(f"Invalid project ID: {pid!r}")
        return pid

    @classmethod
    def project_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid)

    @classmethod
    def uploads_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'uploads')

    @classmethod
    def results_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'results')

    @classmethod
    def intermediate_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'intermediate')

    @classmethod
    def plots_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'plots')

    @classmethod
    def runtime_tmp_dir(cls):
        """Return a writable, project-owned directory for disposable files."""
        path = os.path.abspath(cls.RUNTIME_TMP_DIR)
        os.makedirs(path, exist_ok=True)
        return path

    @classmethod
    def configure_runtime_tmpdir(cls):
        """Route Python and common Unix temporary-file lookups away from ``/tmp``.

        This runs while ``config`` is imported, before analysis libraries are
        loaded.  ``RUNTIME_TMP_DIR`` remains an explicit deployment override
        for installations that keep data on a separate mounted volume.
        """
        path = cls.runtime_tmp_dir()
        for variable in ('TMPDIR', 'TMP', 'TEMP'):
            os.environ[variable] = path
        for variable, cache_dir in (
            ('NUMBA_CACHE_DIR', cls.NUMBA_CACHE_DIR),
            ('MPLCONFIGDIR', cls.MPLCONFIG_DIR),
        ):
            os.makedirs(cache_dir, exist_ok=True)
            os.environ[variable] = cache_dir
        # ``tempfile`` caches its selected directory, so update it explicitly
        # as well as the environment variables.
        tempfile.tempdir = path
        return path

    @classmethod
    def branches_dir(cls, pid):
        """候选分支输出目录，隔离于主线 intermediate/."""
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'branches')

    @classmethod
    def branch_dir(cls, pid, branch_id):
        """单个候选分支目录."""
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'branches', branch_id)

    @classmethod
    def _validate_path(cls, path, pid):
        """验证路径在项目目录内，拒绝符号链接。返回 (is_valid: bool, error: str|None)."""
        cls._validate_pid(pid)
        proj_dir = os.path.realpath(cls.project_dir(pid))
        try:
            real_path = os.path.realpath(path)
        except (OSError, ValueError) as e:
            return False, f"路径解析失败: {e}"
        if os.path.islink(path):
            return False, "不支持符号链接文件"
        if not real_path.startswith(proj_dir + os.sep) and real_path != proj_dir:
            return False, f"路径不在项目目录内: {path}"
        return True, None

    @classmethod
    def sc_batch_source_roots(cls):
        """Return canonical administrator-approved roots for batch 10x input.

        Read the environment on each call so a long-running development server
        can be configured without importing arbitrary filesystem paths into the
        application at module import time.  Empty by default is the safe mode.
        """
        raw = os.environ.get('SC_BATCH_SOURCE_ROOTS', cls.SC_BATCH_SOURCE_ROOTS)
        roots = []
        for item in str(raw or '').split(os.pathsep):
            item = item.strip()
            if not item:
                continue
            try:
                root = os.path.realpath(item)
            except (OSError, ValueError):
                continue
            if os.path.isdir(root) and root not in roots:
                roots.append(root)
        return tuple(roots)

    @classmethod
    def validate_sc_batch_source_path(cls, path, *, require_directory=True):
        """Validate a server-side batch source against ``SC_BATCH_SOURCE_ROOTS``.

        This is deliberately separate from ``_validate_path``: batch input is
        read-only and can live on a mounted data volume, while all platform
        outputs still remain beneath the project directory.
        """
        if not path:
            raise ValueError('缺少服务器数据目录')
        roots = cls.sc_batch_source_roots()
        if not roots:
            raise ValueError(
                '服务器目录批量导入尚未启用；管理员需在 .env 中设置 '
                'SC_BATCH_SOURCE_ROOTS（例如 /home/oelab/data/GJ）'
            )
        try:
            supplied = os.path.abspath(str(path))
            resolved = os.path.realpath(supplied)
        except (OSError, ValueError) as exc:
            raise ValueError(f'服务器目录解析失败: {exc}') from exc
        if os.path.islink(supplied):
            raise ValueError('服务器批量导入不接受符号链接路径')
        allowed = any(resolved == root or resolved.startswith(root + os.sep)
                      for root in roots)
        if not allowed:
            raise ValueError('服务器目录不在允许的 SC_BATCH_SOURCE_ROOTS 下')
        if require_directory and not os.path.isdir(resolved):
            raise ValueError(f'服务器目录不存在或不是目录: {path}')
        return resolved


# Apply the non-root temporary directory before Scanpy, Matplotlib, or worker
# modules can initialize their own caches.
Config.configure_runtime_tmpdir()
