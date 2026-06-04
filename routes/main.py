import psutil
import shutil
from flask import Blueprint, render_template
from models import Project
from config import Config

main_bp = Blueprint('main', __name__)

@main_bp.route('/')
def index():
    projects = Project.get_all()
    disk = shutil.disk_usage(Config.DATA_DIR)
    mem = psutil.virtual_memory()
    return render_template('index.html', projects=projects,
        disk_free_gb=round(disk.free / (1024**3), 1),
        ram_avail_gb=round(mem.available / (1024**3), 1),
        ram_percent=mem.percent)
