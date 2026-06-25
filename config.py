import os

class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or os.urandom(24).hex()

    _BASE_DIR = os.path.dirname(os.path.abspath(__file__))
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

    @classmethod
    def project_dir(cls, pid):
        return os.path.join(cls.DATA_DIR, 'projects', pid)

    @classmethod
    def uploads_dir(cls, pid):
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'uploads')

    @classmethod
    def results_dir(cls, pid):
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'results')

    @classmethod
    def intermediate_dir(cls, pid):
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'intermediate')

    @classmethod
    def plots_dir(cls, pid):
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'plots')
