import json
from flask import Blueprint, render_template, send_file, flash, redirect, url_for
from models import Project, AnalysisTask, ResultFile

results_bp = Blueprint('results', __name__)

@results_bp.route('/<pid>/task/<task_id>')
def task_detail(pid, task_id):
    p = Project.get_by_id(pid)
    t = AnalysisTask.get_by_id(task_id)
    if not p or not t:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    files = ResultFile.get_by_task(task_id)
    plotly_files = [f for f in files if f.file_type == 'plotly_json']
    csv_files = [f for f in files if f.file_type == 'csv']
    result_data = {}
    try:
        result_data = json.loads(t.result_json) if t.result_json else {}
    except Exception:
        pass
    return render_template('analysis_result.html', project=p, task=t,
                          plotly_files=plotly_files, csv_files=csv_files,
                          result_data=result_data)

@results_bp.route('/<pid>/results')
def results_gallery(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    return render_template('results_gallery.html', project=p, tasks=tasks)

@results_bp.route('/<pid>/results/file/<file_id>')
def view_file(pid, file_id):
    f = ResultFile.get_by_id(file_id)
    if not f:
        flash('文件未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    if f.file_type == 'csv':
        return send_file(f.file_path, as_attachment=True)
    return send_file(f.file_path)

@results_bp.route('/<pid>/download/<task_id>')
def download_adata(pid, task_id):
    t = AnalysisTask.get_by_id(task_id)
    if not t or not t.output_adata_path:
        flash('文件未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    return send_file(t.output_adata_path, as_attachment=True)
