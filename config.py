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
