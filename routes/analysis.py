import os
import json
from flask import Blueprint, render_template, request, redirect, url_for, flash
from models import Project, AnalysisTask
from config import Config

analysis_bp = Blueprint('analysis', __name__)

MODULE_LIST = [
    {'name': 'qc', 'display': 'Quality Control', 'desc': 'Filter cells, remove doublets'},
    {'name': 'preprocess', 'display': 'Preprocessing', 'desc': 'Normalize, select HVGs'},
    {'name': 'dimred', 'display': 'Dimensionality Reduction', 'desc': 'PCA, UMAP'},
    {'name': 'batch_correct', 'display': 'Batch Correction', 'desc': 'Harmony, ComBat, SysVI'},
    {'name': 'clustering', 'display': 'Clustering', 'desc': 'Leiden clustering'},
    {'name': 'annotation', 'display': 'Cell Annotation', 'desc': 'Marker-based annotation'},
    {'name': 'deg', 'display': 'Differential Expression', 'desc': 'DEG analysis'},
    {'name': 'trajectory', 'display': 'Trajectory', 'desc': 'Pseudotime analysis'},
    {'name': 'proportion', 'display': 'Proportions', 'desc': 'Cell proportion analysis'},
]

PARAM_SCHEMAS = {
    'qc': [
        {'key': 'mito_perc', 'label': 'Max MT fraction', 'type': 'number', 'default': 0.2, 'step': 0.01},
        {'key': 'nUMIs', 'label': 'Min UMIs', 'type': 'number', 'default': 500},
        {'key': 'detected_genes', 'label': 'Min detected genes', 'type': 'number', 'default': 250},
        {'key': 'batch_key', 'label': 'Batch column', 'type': 'text', 'default': 'batch'},
    ],
    'preprocess': [
        {'key': 'n_top_genes', 'label': 'Number of HVGs', 'type': 'number', 'default': 2000},
        {'key': 'target_sum', 'label': 'Target sum (normalization)', 'type': 'number', 'default': 10000},
    ],
    'dimred': [
        {'key': 'n_comps', 'label': 'PCA components', 'type': 'number', 'default': 50},
        {'key': 'use_mde', 'label': 'Use MDE (faster UMAP)', 'type': 'checkbox', 'default': False},
    ],
    'batch_correct': [
        {'key': 'method', 'label': 'Method', 'type': 'select', 'options': ['harmony', 'combat', 'sysvi'], 'default': 'harmony'},
        {'key': 'batch_key', 'label': 'Batch column', 'type': 'text', 'default': 'batch'},
        {'key': 'n_pcs', 'label': 'Number of PCs', 'type': 'number', 'default': 50},
        {'key': 'max_epochs', 'label': 'Max epochs (SysVI only)', 'type': 'number', 'default': 200},
    ],
    'clustering': [
        {'key': 'resolutions', 'label': 'Leiden resolutions (comma-separated)', 'type': 'text', 'default': '0.6,0.8,1.0'},
        {'key': 'n_neighbors', 'label': 'Number of neighbors', 'type': 'number', 'default': 15},
    ],
    'annotation': [
        {'key': 'method', 'label': 'Method', 'type': 'select', 'options': ['auto_marker', 'manual'], 'default': 'auto_marker'},
        {'key': 'cluster_key', 'label': 'Cluster column', 'type': 'text', 'default': 'leiden'},
        {'key': 'resolution', 'label': 'Leiden resolution', 'type': 'text', 'default': '0.8'},
    ],
    'deg': [
        {'key': 'groupby', 'label': 'Group by', 'type': 'text', 'default': 'celltype'},
        {'key': 'reference', 'label': 'Reference', 'type': 'text', 'default': 'rest'},
        {'key': 'method', 'label': 'Method', 'type': 'select', 'options': ['wilcoxon', 't-test', 'logreg'], 'default': 'wilcoxon'},
        {'key': 'n_genes', 'label': 'Top N genes to show', 'type': 'number', 'default': 20},
    ],
    'trajectory': [
        {'key': 'method', 'label': 'Method', 'type': 'select', 'options': ['diffusion_map', 'slingshot'], 'default': 'diffusion_map'},
        {'key': 'cluster_key', 'label': 'Cluster column', 'type': 'text', 'default': 'leiden'},
    ],
    'proportion': [
        {'key': 'groupby', 'label': 'Group by', 'type': 'text', 'default': 'celltype'},
        {'key': 'batch_key', 'label': 'Batch column', 'type': 'text', 'default': 'batch'},
    ],
}

@analysis_bp.route('/<pid>/analyze/<module_name>', methods=['GET', 'POST'])
def analyze(pid, module_name):
    p = Project.get_by_id(pid)
    if not p:
        flash('Project not found', 'danger')
        return redirect(url_for('main.index'))
    mod_info = next((m for m in MODULE_LIST if m['name'] == module_name), None)
    if not mod_info:
        flash(f'Module "{module_name}" not found', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    schema = PARAM_SCHEMAS.get(module_name, [])
    tasks = AnalysisTask.get_by_project(pid)
    completed_tasks = [t for t in tasks if t.status == 'completed' and t.output_adata_path]
    if request.method == 'POST':
        params = {}
        for param in schema:
            val = request.form.get(param['key'])
            if param['type'] == 'number':
                params[param['key']] = float(val) if val else param['default']
            elif param['type'] == 'checkbox':
                params[param['key']] = request.form.get(param['key']) == 'on'
            else:
                params[param['key']] = val or param['default']
        input_path = request.form.get('input_path', '')
        if not input_path:
            flash('Please select input data', 'danger')
            return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
        task = AnalysisTask(project_id=pid, module_name=module_name,
                           params_json=json.dumps(params))
        task.save()
        p.status = 'processing'
        p.save()
        from worker import submit_task
        submit_task(task.id, pid, module_name, params,
                   os.path.join(Config.DATA_DIR, 'projects', pid), input_path)
        return redirect(url_for('results.task_detail', pid=pid, task_id=task.id))
    uploads_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'uploads')
    uploaded_h5ad = []
    if os.path.isdir(uploads_dir):
        for f in os.listdir(uploads_dir):
            if f.endswith('.h5ad'):
                uploaded_h5ad.append({'name': f, 'path': os.path.join(uploads_dir, f)})
    return render_template('analysis_select.html', project=p, module=mod_info,
                          schema=schema, completed_tasks=completed_tasks,
                          uploaded_h5ad=uploaded_h5ad, all_modules=MODULE_LIST)
