import os
import shutil
from flask import Blueprint, render_template, request, redirect, url_for, flash
from models import Project, AnalysisTask, ResultFile
from config import Config
from routes.analysis import MODULE_DISPLAY_MAP, STATUS_MAP

projects_bp = Blueprint('projects', __name__)

@projects_bp.route('/new', methods=['GET', 'POST'])
def new():
    if request.method == 'POST':
        name = request.form.get('name', '').strip()
        if not name:
            flash('请输入项目名称', 'danger')
            return redirect(url_for('projects.new'))
        p = Project(name=name, description=request.form.get('description', ''))
        proj_dir = os.path.join(Config.DATA_DIR, 'projects', p.id)
        for sub in ['uploads', 'intermediate', 'results', 'plots']:
            os.makedirs(os.path.join(proj_dir, sub), exist_ok=True)
        p.save()
        flash(f'项目 "{name}" 已创建', 'success')
        return redirect(url_for('projects.detail', pid=p.id))
    return render_template('project_new.html')

@projects_bp.route('/<pid>')
def detail(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = p.get_tasks()
    files = ResultFile.get_by_project(pid)
    uploads_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'uploads')
    uploaded = []
    if os.path.isdir(uploads_dir):
        for f in os.listdir(uploads_dir):
            fpath = os.path.join(uploads_dir, f)
            if os.path.isfile(fpath):
                uploaded.append({'name': f, 'size_mb': round(os.path.getsize(fpath) / (1024**2), 1)})
    return render_template('project_detail.html', project=p, tasks=tasks, result_files=files, uploaded=uploaded,
                          module_display_map=MODULE_DISPLAY_MAP, status_map=STATUS_MAP)

@projects_bp.route('/<pid>/delete', methods=['POST'])
def delete(pid):
    p = Project.get_by_id(pid)
    if p:
        proj_dir = os.path.join(Config.DATA_DIR, 'projects', pid)
        if os.path.isdir(proj_dir):
            shutil.rmtree(proj_dir)
        p.delete()
        flash('项目已删除', 'success')
    return redirect(url_for('main.index'))
