import os

class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY', 'bioinfo-lab-dev-key')
    DATA_DIR = '/data/GJ/platform/data'
    DB_PATH = '/data/GJ/platform/instance/bioinfo.db'
    CELLMARKER_PATH = '/data/GJ/CellMarker_Augmented_2021.txt'
    MAX_WORKERS = 2
    CHUNK_SIZE_MB = 50
    PLOTLY_MAX_CELLS = 50000
    MIN_FREE_RAM_GB = 50
    CUDA_DEVICES = '0,1'

    # AI 对话配置（支持任意 OpenAI 兼容 API）
    AI_API_KEY = os.environ.get('AI_API_KEY', '')
    AI_API_URL = os.environ.get('AI_API_URL', 'https://api.anthropic.com/v1')
    AI_MODEL = os.environ.get('AI_MODEL', 'claude-sonnet-4-20250514')
