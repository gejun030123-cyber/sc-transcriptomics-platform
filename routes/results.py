import json
import os
from flask import Blueprint, render_template, send_file, flash, redirect, url_for
from models import Project, AnalysisTask, ResultFile
from routes.analysis import MODULE_DISPLAY_MAP, STATUS_MAP
from modules import PIPELINE_ORDER
from config import Config


def _validate_path(file_path):
    if not file_path:
        return False
    abs_path = os.path.abspath(file_path)
    data_dir = os.path.abspath(Config.DATA_DIR)
    return abs_path.startswith(data_dir + os.sep)


results_bp = Blueprint('results', __name__)

@results_bp.route('/<pid>/task/<task_id>')
def task_detail(pid, task_id):
    p = Project.get_by_id(pid)
    t = AnalysisTask.get_by_id(task_id)
    if not p or not t:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    files = ResultFile.get_by_task(task_id)
    for f in files:
        try:
            f._file_size = os.path.getsize(f.file_path) if f.file_path and os.path.exists(f.file_path) else 0
        except Exception:
            f._file_size = 0
    plotly_files = [f for f in files if f.file_type == 'plotly_json']
    image_files = [f for f in files if f.file_type in ('png', 'svg', 'jpg', 'jpeg')]
    csv_files = [f for f in files if f.file_type == 'csv']
    result_data = {}
    try:
        result_data = json.loads(t.result_json) if t.result_json else {}
    except Exception:
        pass

    # 计算下一步模块
    next_module = None
    current_module = t.module_name
    try:
        idx = PIPELINE_ORDER.index(current_module)
        if idx < len(PIPELINE_ORDER) - 1:
            next_module = PIPELINE_ORDER[idx + 1]
    except ValueError:
        pass

    return render_template('analysis_result.html', project=p, task=t,
                          plotly_files=plotly_files, image_files=image_files, csv_files=csv_files,
                          result_data=result_data,
                          module_display=MODULE_DISPLAY_MAP.get(t.module_name, t.module_name),
                          status_cn=STATUS_MAP.get(t.status, t.status),
                          next_module=next_module,
                          module_display_map=MODULE_DISPLAY_MAP)

@results_bp.route('/<pid>/results')
def results_gallery(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    return render_template('results_gallery.html', project=p, tasks=tasks,
                          module_display_map=MODULE_DISPLAY_MAP, status_map=STATUS_MAP)

@results_bp.route('/<pid>/results/file/<file_id>')
def view_file(pid, file_id):
    f = ResultFile.get_by_id(file_id)
    if not f:
        flash('文件未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    if not _validate_path(f.file_path):
        flash('文件路径不合法', 'danger')
        return redirect(url_for('main.index'))
    if f.file_type == 'csv':
        return send_file(f.file_path, as_attachment=True)
    return send_file(f.file_path)

@results_bp.route('/<pid>/download/<task_id>')
def download_adata(pid, task_id):
    t = AnalysisTask.get_by_id(task_id)
    if not t or not t.output_adata_path:
        flash('文件未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    if not _validate_path(t.output_adata_path):
        flash('文件路径不合法', 'danger')
        return redirect(url_for('main.index'))
    return send_file(t.output_adata_path, as_attachment=True)
