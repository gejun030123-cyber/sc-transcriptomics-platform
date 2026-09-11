import psutil
import shutil
from flask import Blueprint, redirect, render_template, url_for
from models import Project
from config import Config

main_bp = Blueprint('main', __name__)

@main_bp.route('/')
def index():
    """Open the most recent research project directly in the AI workspace.

    A chat without an active project cannot safely carry assay, sample, input
    and result context.  For that reason the platform enters the latest
    project workspace by default; a no-project state keeps the same agent
    shell and guides the researcher to create that context first.
    """
    projects = Project.get_all()
    if projects:
        return redirect(url_for('workspace.index', pid=projects[0].id))
    return render_template('workspace_onboarding.html')


@main_bp.route('/projects')
def project_dashboard():
    """Project management remains available without competing with the AI home."""
    projects = Project.get_all()
    disk = shutil.disk_usage(Config.DATA_DIR)
    mem = psutil.virtual_memory()
    return render_template('index.html', projects=projects,
        disk_free_gb=round(disk.free / (1024**3), 1),
        ram_avail_gb=round(mem.available / (1024**3), 1),
        ram_percent=mem.percent)
