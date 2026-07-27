import os
import json
import copy
import re
from flask import Blueprint, render_template, request, redirect, url_for, flash
from models import Project, AnalysisTask, ResultFile
from config import Config
from modules.schemas import (
    SC_MODULE_LIST, BULK_MODULE_LIST, MODULE_LIST,
    MODULE_DISPLAY_MAP, SC_MODULE_NAMES, BULK_MODULE_NAMES,
    STATUS_MAP, PARAM_SCHEMAS, parse_form_params, list_upload_files,
)

analysis_bp = Blueprint('analysis', __name__)


def _comparison_label_for_deg_result(result_file, task):
    """Recover a human-readable contrast label from a registered DEG table."""
    label = str(result_file.label or '').strip()
    match = re.search(r'\(([^()]*\bvs\b[^()]*)\)', label, flags=re.IGNORECASE)
    if match:
        return match.group(1).strip()
    suffix_match = re.search(r'bulk_deg_results_(\d+)\.csv$', os.path.basename(result_file.file_path), flags=re.I)
    if suffix_match:
        try:
            result = json.loads(task.result_json or '{}')
            comparisons = list(result.get('comparisons') or [])
            index = int(suffix_match.group(1))
            if 0 <= index < len(comparisons):
                return str(comparisons[index]).strip()
        except (TypeError, ValueError, json.JSONDecodeError):
            pass
    return label or os.path.basename(result_file.file_path)


def _bulk_deg_enrichment_sources(project_id):
    """List valid single-comparison DEG tables for ORA/GSEA selection.

    Merged multi-comparison tables are deliberately excluded because combining
    contrasts would make a single enrichment result biologically ambiguous.
    """
    options = []
    for result_file in ResultFile.get_by_project(project_id):
        if result_file.file_type != 'csv':
            continue
        filename = os.path.basename(result_file.file_path).lower()
        if not filename.startswith('bulk_deg_results') or not filename.endswith('.csv'):
            continue
        task = AnalysisTask.get_by_id(result_file.task_id)
        if not task or task.module_name != 'bulk_deg' or task.status != 'completed':
            continue
        valid, _ = Config._validate_path(result_file.file_path, project_id)
        if not valid or not os.path.isfile(result_file.file_path):
            continue
        options.append({
            'path': result_file.file_path,
            'label': _comparison_label_for_deg_result(result_file, task),
            'comparison': _comparison_label_for_deg_result(result_file, task),
            'task_id': task.id,
        })
    return sorted(options, key=lambda item: (item['label'], item['path']))


def _path_has_module_requirements(path, requirements):
    """Check the lightweight AnnData contract needed by a module input."""
    if not requirements:
        return True
    if not path or not str(path).lower().endswith('.h5ad') or not os.path.isfile(path):
        return False
    adata = None
    try:
        import anndata
        adata = anndata.read_h5ad(path, backed='r')
        for requirement in requirements:
            if str(requirement).startswith('X_'):
                if requirement not in adata.obsm:
                    return False
            elif requirement not in adata.obs.columns:
                return False
        return True
    except Exception:
        return False
    finally:
        try:
            if adata is not None and getattr(adata, 'isbacked', False):
                adata.file.close()
        except Exception:
            pass


def build_input_options(module_name, completed_tasks, uploaded_files):
    """Build an ordered, dependency-aware input list for the analysis form.

    The previous UI showed every path without provenance, making it easy to
    accidentally feed a downstream output back into an earlier module.  Keep
    valid upstream history available, but select the newest compatible output
    by default and explain why it is recommended.
    """
    from modules import MODULE_REGISTRY, PIPELINE_ORDER

    module_index = {name: index for index, name in enumerate(PIPELINE_ORDER)}
    current_index = module_index.get(module_name, len(PIPELINE_ORDER))
    module_cls = MODULE_REGISTRY.get(module_name)
    requirements = list(getattr(module_cls, 'INPUT_REQUIRES', []) or [])
    preferred_modules = {
        'qc': ('convert_10x',),
        'normalize': ('qc',),
        'hvg': ('normalize',),
        'dimred': ('hvg',),
        'batch_correct': ('dimred', 'hvg'),
        'clustering': ('dimred',),
        'subcluster': ('clustering',),
        'qc_reassess': ('clustering',),
        'annotation': ('qc_reassess', 'clustering'),
        'sc_timecourse': ('annotation', 'qc_reassess', 'clustering'),
        'deg': ('annotation', 'qc_reassess', 'clustering'),
        'trajectory': ('qc_reassess', 'clustering'),
        'proportion': ('annotation', 'qc_reassess', 'clustering'),
        'cell_communication': ('annotation',),
        'bulk_normalize': ('bulk_qc',),
        'bulk_pca': ('bulk_normalize',),
        'bulk_deg': ('bulk_normalize',),
        'bulk_heatmap': ('bulk_deg',),
        'bulk_enrichment': ('bulk_deg',),
        'bulk_timecourse': ('bulk_normalize',),
        'bulk_deg_integration': ('bulk_deg',),
    }.get(module_name, ())
    preference_rank = {name: index for index, name in enumerate(preferred_modules)}

    def task_time(task):
        return str(task.finished_at or task.started_at or '')

    upstream_tasks = []
    seen_paths = set()
    for task in completed_tasks:
        if task.module_name == module_name or not task.output_adata_path:
            continue
        path = os.path.abspath(task.output_adata_path)
        if path in seen_paths or not os.path.isfile(path):
            continue
        task_index = module_index.get(task.module_name, -1)
        if task.module_name not in preference_rank and (task_index < 0 or task_index >= current_index):
            continue
        if not _path_has_module_requirements(path, requirements):
            continue
        seen_paths.add(path)
        upstream_tasks.append((task, path, preference_rank.get(task.module_name, len(preferred_modules))))

    # The timestamp is an ISO string, so sort preferred modules first and then
    # newest output within the same module without relying on database order.
    ordered_tasks = []
    for rank in sorted({item[2] for item in upstream_tasks}):
        same_rank = [item for item in upstream_tasks if item[2] == rank]
        same_rank.sort(key=lambda item: task_time(item[0]), reverse=True)
        ordered_tasks.extend(same_rank)
    upstream_tasks = ordered_tasks

    options = []
    if upstream_tasks:
        recommended_task, recommended_path, _ = upstream_tasks[0]
        options.append({
            'path': recommended_path,
            'name': os.path.basename(recommended_path),
            'source': 'task',
            'task_id': recommended_task.id,
            'module_name': recommended_task.module_name,
            'module_label': MODULE_DISPLAY_MAP.get(recommended_task.module_name, recommended_task.module_name),
            'finished_at': task_time(recommended_task),
            'recommended': True,
            'reason': (
                f"来自“{MODULE_DISPLAY_MAP.get(recommended_task.module_name, recommended_task.module_name)}”的最新输出"
                + (f"，满足本模块所需字段：{', '.join(requirements)}" if requirements else '')
            ),
        })

    for task, path, _ in upstream_tasks:
        if options and path == options[0]['path']:
            continue
        options.append({
            'path': path,
            'name': os.path.basename(path),
            'source': 'task',
            'task_id': task.id,
            'module_name': task.module_name,
            'module_label': MODULE_DISPLAY_MAP.get(task.module_name, task.module_name),
            'finished_at': task_time(task),
            'recommended': False,
            'reason': '历史上游分析输出，可手动切换',
        })

    for uploaded in uploaded_files:
        path = os.path.abspath(uploaded.get('path', ''))
        if not path or path in seen_paths or not os.path.isfile(path):
            continue
        seen_paths.add(path)
        options.append({
            'path': path,
            'name': uploaded.get('name') or os.path.basename(path),
            'source': 'upload',
            'task_id': '',
            'module_name': '',
            'module_label': '上传文件',
            'finished_at': '',
            'recommended': not options,
            'reason': '没有找到满足依赖的上游输出，从上传文件开始',
        })

    # If a module has no prior task, keep the first upload as the default.
    if options and not any(item['recommended'] for item in options):
        options[0]['recommended'] = True
        options[0]['reason'] = options[0]['reason'] or '默认使用此项目输入文件'

    recommended = next((item for item in options if item['recommended']), None)
    return options, recommended


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
    # Defaults can be project-specific (for example a single available DEG
    # table), so never mutate the shared schema object.
    schema = copy.deepcopy(PARAM_SCHEMAS.get(module_name, []))
    tasks = AnalysisTask.get_by_project(pid)
    completed_tasks = [t for t in tasks if t.status == 'completed' and t.output_adata_path]

    is_bulk = module_name in BULK_MODULE_NAMES
    uploaded_files = list_upload_files(pid, is_bulk)
    input_options, recommended_input = build_input_options(
        module_name, completed_tasks, uploaded_files,
    )

    if request.method == 'POST':
        params = parse_form_params(schema, request.form)
        if module_name == 'bulk_enrichment':
            custom_genes = str(params.get('custom_genes', '')).strip()
            selected_deg = str(params.get('input_source', '')).strip()
            if not custom_genes and not selected_deg:
                flash('请选择一个 DEG 比较结果后再进行通路富集；多比较结果不能合并为一次富集。', 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
            if not custom_genes:
                selected_source = next(
                    (item for item in _bulk_deg_enrichment_sources(pid)
                     if os.path.abspath(item['path']) == os.path.abspath(selected_deg)),
                    None,
                )
                if selected_source is None:
                    flash('请选择页面列出的单个 DEG 比较结果。', 'danger')
                    return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
                params['input_comparison'] = selected_source['comparison']
        if module_name == 'bulk_heatmap' and params.get('gene_import_source') == 'deg':
            selected_deg = str(params.get('deg_comparison', '')).strip()
            if selected_deg:
                selected_source = next(
                    (item for item in _bulk_deg_enrichment_sources(pid)
                     if os.path.abspath(item['path']) == os.path.abspath(selected_deg)),
                    None,
                )
                if selected_source is None:
                    flash('请选择页面列出的单个 DEG 比较结果。', 'danger')
                    return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
                params['deg_comparison_label'] = selected_source['comparison']
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
        allowed_input_paths = {os.path.abspath(item['path']) for item in input_options}
        if allowed_input_paths and os.path.abspath(input_path) not in allowed_input_paths:
            flash('请选择页面列出的输入文件；系统已过滤失效或不属于当前项目的路径。', 'danger')
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

        # Keep the normal form on the same hard experimental-design boundary
        # shown by the preflight card.  Only prerequisites that would make the
        # module fail or invalidate its declared statistical method are
        # blocked here; descriptive/exploratory warnings remain runnable.
        if module_name in {'bulk_deg', 'sc_timecourse', 'batch_correct', 'bulk_normalize'}:
            try:
                from modules.design_preflight import preflight_blockers
                from modules.io_utils import read_expression_matrix

                blockers = preflight_blockers(
                    read_expression_matrix(input_path), module_name, params, input_path,
                )
            except Exception:
                # The worker will retain its existing detailed parser/input
                # validation.  Do not replace a module-specific error with a
                # generic preflight parsing failure.
                blockers = []
            if blockers:
                flash('分析前检查未通过：' + '；'.join(blockers[:2]), 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
        task = AnalysisTask(project_id=pid, module_name=module_name,
                           params_json=json.dumps(params))
        task.save()
        p.status = 'processing'
        p.save()
        from worker import submit_task
        submit_task(task.id, pid, module_name, params,
                   Config.project_dir(pid), input_path)
        return redirect(url_for('results.task_detail', pid=pid, task_id=task.id))

    enrichment_deg_sources = []
    heatmap_deg_sources = []
    # List every valid per-comparison DEG table.  When several contrasts are
    # present, the user must choose one rather than silently enriching an
    # arbitrary or merged result.
    if module_name == 'bulk_enrichment':
        enrichment_deg_sources = _bulk_deg_enrichment_sources(pid)
        if len(enrichment_deg_sources) == 1:
            for param in schema:
                if param['key'] == 'input_source':
                    param['default'] = enrichment_deg_sources[0]['path']
    if module_name == 'bulk_heatmap':
        heatmap_deg_sources = _bulk_deg_enrichment_sources(pid)

    sidebar_modules = BULK_MODULE_LIST if is_bulk else SC_MODULE_LIST
    return render_template('analysis_select.html', project=p, module=mod_info,
                          schema=schema, completed_tasks=completed_tasks,
                          uploaded_h5ad=uploaded_files, input_options=input_options,
                          recommended_input=recommended_input,
                          enrichment_deg_sources=enrichment_deg_sources,
                          heatmap_deg_sources=heatmap_deg_sources,
                          all_modules=sidebar_modules,
                          module_display_map=MODULE_DISPLAY_MAP)
