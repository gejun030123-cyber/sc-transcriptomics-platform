import os
import re


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
    DB_PATH = os.environ.get('DB_PATH', os.path.join(_BASE_DIR, 'instance', 'bioinfo.db'))
    CELLMARKER_PATH = os.environ.get('CELLMARKER_PATH', os.path.join(os.path.dirname(_BASE_DIR), 'CellMarker_Augmented_2021.txt'))

    MAX_WORKERS = int(os.environ.get('MAX_WORKERS', '2'))
    CHUNK_SIZE_MB = 50
    PLOTLY_MAX_CELLS = 50000
    MIN_FREE_RAM_GB = float(os.environ.get('MIN_FREE_RAM_GB', '4'))
    CUDA_DEVICES = os.environ.get('CUDA_DEVICES', '0,1')

    # AI 对话配置（支持任意 OpenAI 兼容 API）
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
