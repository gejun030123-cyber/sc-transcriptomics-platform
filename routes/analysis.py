import os
import json
from flask import Blueprint, render_template, request, redirect, url_for, flash
from models import Project, AnalysisTask
from config import Config
from modules.schemas import (
    SC_MODULE_LIST, BULK_MODULE_LIST, MODULE_LIST,
    MODULE_DISPLAY_MAP, SC_MODULE_NAMES, BULK_MODULE_NAMES,
    STATUS_MAP, PARAM_SCHEMAS, parse_form_params, list_upload_files,
)

analysis_bp = Blueprint('analysis', __name__)


@analysis_bp.route('/<pid>/sc-analysis')
def sc_analysis_list(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    return render_template('sc_analysis.html', project=p, modules=SC_MODULE_LIST, tasks=tasks)


@analysis_bp.route('/<pid>/bulk-analysis')
def bulk_analysis_list(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    return render_template('bulk_analysis.html', project=p, modules=BULK_MODULE_LIST, tasks=tasks)


@analysis_bp.route('/<pid>/analyze/<module_name>', methods=['GET', 'POST'])
def analyze(pid, module_name):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    mod_info = next((m for m in MODULE_LIST if m['name'] == module_name), None)
    if not mod_info:
        flash(f'模块 "{module_name}" 不存在', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    schema = PARAM_SCHEMAS.get(module_name, [])
    tasks = AnalysisTask.get_by_project(pid)
    completed_tasks = [t for t in tasks if t.status == 'completed' and t.output_adata_path]

    is_bulk = module_name in BULK_MODULE_NAMES

    if request.method == 'POST':
        params = parse_form_params(schema, request.form)
        # 处理自动检测的分组映射
        auto_mapping = request.form.get('_auto_group_mapping')
        if auto_mapping:
            try:
                params['_auto_group_mapping'] = json.loads(auto_mapping)
            except Exception:
                pass
        input_path = request.form.get('input_path', '')
        if not input_path:
            flash('请选择输入数据', 'danger')
            return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
        is_valid, err = Config._validate_path(input_path, pid)
        if not is_valid or not os.path.isfile(input_path):
            flash(err or '输入文件不存在', 'danger')
            return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
        # 注入 _visualization 和 _filters 到 params
        viz_json = request.form.get('_visualization', '')
        filters_json = request.form.get('_filters', '')
        if viz_json:
            try:
                params['_visualization'] = json.loads(viz_json)
            except json.JSONDecodeError:
                pass
        if filters_json:
            try:
                params['_filters'] = json.loads(filters_json)
            except json.JSONDecodeError:
                pass
        task = AnalysisTask(project_id=pid, module_name=module_name,
                           params_json=json.dumps(params))
        task.save()
        p.status = 'processing'
        p.save()
        from worker import submit_task
        submit_task(task.id, pid, module_name, params,
                   Config.project_dir(pid), input_path)
        return redirect(url_for('results.task_detail', pid=pid, task_id=task.id))

    uploaded_files = list_upload_files(pid, is_bulk)

    # For enrichment module, auto-populate input_source with DEG results
    if module_name == 'bulk_enrichment':
        results_dir = Config.results_dir(pid)
        deg_csv = os.path.join(results_dir, 'bulk_deg_results.csv')
        if os.path.exists(deg_csv):
            for param in schema:
                if param['key'] == 'input_source':
                    param['default'] = deg_csv

    sidebar_modules = BULK_MODULE_LIST if is_bulk else SC_MODULE_LIST
    return render_template('analysis_select.html', project=p, module=mod_info,
                          schema=schema, completed_tasks=completed_tasks,
                          uploaded_h5ad=uploaded_files, all_modules=sidebar_modules,
                          module_display_map=MODULE_DISPLAY_MAP)
