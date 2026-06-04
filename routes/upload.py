import os
from flask import Blueprint, render_template, request, redirect, url_for, flash, jsonify
from werkzeug.utils import secure_filename
from models import Project
from config import Config

upload_bp = Blueprint('upload', __name__)

ALLOWED_EXT = {'.h5ad', '.h5', '.csv', '.txt', '.mtx', '.gz'}

@upload_bp.route('/<pid>/upload', methods=['GET', 'POST'])
def upload(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('请选择文件', 'danger')
            return redirect(url_for('upload.upload', pid=pid))
        file = request.files['file']
        if file.filename == '':
            flash('请选择文件', 'danger')
            return redirect(url_for('upload.upload', pid=pid))
        fname = secure_filename(file.filename)
        ext = os.path.splitext(fname)[1].lower()
        if ext not in ALLOWED_EXT and not fname.endswith('.mtx.gz') and not fname.endswith('.tsv.gz'):
            flash(f'不支持的文件格式: {ext}', 'danger')
            return redirect(url_for('upload.upload', pid=pid))
        uploads_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'uploads')
        os.makedirs(uploads_dir, exist_ok=True)
        fpath = os.path.join(uploads_dir, fname)
        file.save(fpath)
        p.status = 'data_ready'
        p.save()
        flash(f'文件 "{fname}" 上传成功', 'success')
        return redirect(url_for('projects.detail', pid=pid))
    return render_template('upload.html', project=p)
