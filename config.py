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
