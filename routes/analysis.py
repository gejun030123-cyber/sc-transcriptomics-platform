import os
import json
import copy
import re
import pandas as pd
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
        try:
            task_params = json.loads(task.params_json or '{}')
        except (TypeError, json.JSONDecodeError):
            task_params = {}
        try:
            task_summary = json.loads(task.result_json or '{}')
        except (TypeError, json.JSONDecodeError):
            task_summary = {}
        source_groupby = str(
            task_summary.get('groupby') or task_params.get('groupby') or ''
        ).strip()
        if source_groupby.lower() in {'无', 'none', 'null'}:
            source_groupby = ''
        options.append({
            'path': result_file.file_path,
            'label': _comparison_label_for_deg_result(result_file, task),
            'comparison': _comparison_label_for_deg_result(result_file, task),
            'groupby': source_groupby,
            'task_id': task.id,
        })
    return sorted(options, key=lambda item: (item['label'], item['path']))


def _sc_deg_enrichment_sources(project_id):
    """List completed, task-bound single-cell DEG sources for enrichment.

    Unlike the retired package-directory scan, this follows the DEG task's
    recorded internal source files.  It prevents a newer unrelated CSV from
    silently changing the biological contrast used for pathway analysis.
    """
    options = []
    for task in AnalysisTask.get_by_project(project_id):
        if task.status != 'completed' or task.module_name not in {
            'deg', 'sc_cell_deg', 'sc_pseudobulk_deg',
        }:
            continue
        try:
            summary = json.loads(task.result_json or '{}')
        except (TypeError, ValueError, json.JSONDecodeError):
            continue
        source_files = summary.get('deg_source_files') or []
        if isinstance(source_files, str):
            source_files = [source_files]
        valid_files = []
        for source_path in source_files:
            valid, _ = Config._validate_path(str(source_path), project_id)
            if valid and os.path.isfile(source_path):
                valid_files.append(str(source_path))
        if not valid_files:
            continue
        clusters = set()
        for source_path in valid_files:
            try:
                header = pd.read_csv(source_path, nrows=0)
                if "cluster" not in header.columns:
                    clusters.add("All")
                    continue
                values = pd.read_csv(source_path, usecols=["cluster"])["cluster"]
                clusters.update(
                    str(value).strip() for value in values.dropna()
                    if str(value).strip()
                )
            except (OSError, ValueError, pd.errors.ParserError):
                # The worker revalidates the exact DEG files before execution.
                # Do not expose an unreadable historical file as a selectable
                # cluster list in the run form.
                continue
        if not clusters:
            clusters.add("All")
        source_level = str(summary.get('deg_source_level') or '').strip().lower()
        if source_level not in {'cell_level', 'pseudobulk'}:
            continue
        contract = str(summary.get('deg_source_task_contract') or '').strip()
        label_level = '样本级 pseudobulk' if source_level == 'pseudobulk' else '细胞级探索性'
        grouping_mode = str(summary.get('analysis_grouping') or '').strip()
        group_label = str(summary.get('analysis_group_label') or '').strip()
        if not group_label:
            group_label = '细胞类型' if grouping_mode == 'annotated_celltype' else 'cluster'
        options.append({
            'task_id': task.id,
            'files': valid_files,
            'source_level': source_level,
            'contract': contract,
            'clusters': sorted(clusters, key=lambda value: (value != 'All', value)),
            'group_label': group_label,
            'grouping_mode': grouping_mode,
            'label': f'{label_level} · {group_label} 分组 · {contract or task.module_name} [{task.id}]',
            'output_adata_path': os.path.abspath(task.output_adata_path or ''),
        })
    return sorted(options, key=lambda item: item['task_id'], reverse=True)


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


BASE_GRN_EXTENSIONS = (
    '.parquet', '.pq', '.csv', '.tsv', '.txt', '.gz',
    '.pickle', '.pkl', '.gpickle', '.oracle', '.celloracle',
)


def _list_base_grn_files(pid):
    """列出项目 uploads 中可作为 CellOracle base GRN 的文件。"""
    from config import Config
    uploads_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'uploads')
    files = []
    if os.path.isdir(uploads_dir):
        for name in sorted(os.listdir(uploads_dir)):
            if name.lower().endswith(BASE_GRN_EXTENSIONS):
                files.append({
                    'name': name,
                    'path': os.path.join(uploads_dir, name),
                })
    return files


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
        # Clustering with the web default ``use_corrected=true`` must consume
        # the batch-corrected output when one exists.  Keeping dimred as the
        # fallback preserves the explicit PCA-only path when no correction was
        # run (or when the user disables ``use_corrected``).
        'clustering': ('batch_correct', 'dimred'),
        'subcluster': ('clustering',),
        'qc_reassess': ('clustering',),
        'annotation': ('qc_reassess', 'clustering'),
        'functional_state': ('annotation',),
        'scenic': ('annotation',),
        'sc_timecourse': ('annotation', 'qc_reassess', 'clustering'),
        'deg': ('annotation', 'qc_reassess', 'clustering'),
        # The default DE contract is condition testing within reviewed cell
        # types.  Cluster outputs remain valid fallbacks for the explicit
        # compatibility mode, but do not hide an available annotation output.
        'sc_cell_deg': ('annotation', 'qc_reassess', 'clustering'),
        'sc_pseudobulk_deg': ('annotation', 'qc_reassess', 'clustering'),
        'sc_cell_go': ('sc_pseudobulk_deg', 'sc_cell_deg', 'deg', 'annotation', 'qc_reassess', 'clustering'),
        'trajectory': ('qc_reassess', 'clustering'),
        'proportion': ('annotation', 'qc_reassess', 'clustering'),
        'neighborhood_da': ('annotation', 'qc_reassess', 'clustering'),
        'cell_communication': ('annotation',),
        'virtual_ko': ('annotation', 'qc_reassess', 'clustering'),
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
        selected_sc_deg_source = None
        if module_name == 'bulk_enrichment':
            method = str(params.get('method', 'ORA') or 'ORA').upper()
            database = str(params.get('database', 'GO_BP') or 'GO_BP')
            batch_databases = bool(params.get('batch_databases', False))
            selected_databases = str(params.get('databases', '') or '').strip()
            custom_genes = str(params.get('custom_genes', '')).strip()
            custom_geneset = str(params.get('custom_geneset_text', '')).strip()
            selected_deg = str(params.get('input_source', '')).strip()
            if batch_databases and not selected_databases:
                flash('请至少选择一个批量基因集数据库。', 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
            if database == 'Custom_GMT' and not custom_geneset:
                flash('选择 Custom_GMT 时必须粘贴 Human GMT 内容。', 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
            if method == 'GSEA' and custom_genes:
                flash('GSEA 需要完整排序 DEG 表，不能使用无排序的自定义基因列表。', 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
            if not custom_genes and not selected_deg:
                flash('请选择一个 DEG 比较结果后再进行通路富集；多比较结果不能合并为一次富集。', 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
            if selected_deg:
                selected_source = next(
                    (item for item in _bulk_deg_enrichment_sources(pid)
                     if os.path.abspath(item['path']) == os.path.abspath(selected_deg)),
                    None,
                )
                if selected_source is None:
                    flash('请选择页面列出的单个 DEG 比较结果。', 'danger')
                    return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
                if custom_genes:
                    params['input_comparison'] = 'Custom genes'
                    params['background_comparison'] = selected_source['comparison']
                else:
                    params['input_comparison'] = selected_source['comparison']
        if module_name == 'bulk_heatmap' and params.get('gene_import_source') in {'top_var', 'deg'}:
            gene_source = params.get('gene_import_source')
            selected_deg = str(params.get('deg_comparison', '')).strip()
            if gene_source == 'deg' and not selected_deg:
                flash('使用 DEG 基因来源时，请明确选择一个 DEG 比较结果。', 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
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
                if not str(params.get('groupby', '') or '').strip() and selected_source.get('groupby'):
                    params['groupby'] = selected_source['groupby']
        if module_name == 'sc_cell_go':
            source_task_id = str(params.get('deg_source_task_id', '') or '').strip()
            selected_sc_deg_source = next(
                (item for item in _sc_deg_enrichment_sources(pid)
                 if item['task_id'] == source_task_id),
                None,
            )
            if selected_sc_deg_source is None:
                flash('请明确选择一个已完成的单细胞 DEG 任务后再运行富集。', 'danger')
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
            params['deg_source_files'] = selected_sc_deg_source['files']
            params['deg_source_level'] = selected_sc_deg_source['source_level']
            params['source_level'] = selected_sc_deg_source['source_level']
            params['source_task_id'] = selected_sc_deg_source['task_id']
            requested_clusters = [
                item.strip() for item in re.split(
                    r'[,;\n]+', str(params.get('target_clusters', '') or ''),
                ) if item.strip()
            ]
            unknown_clusters = sorted(
                set(requested_clusters) - set(selected_sc_deg_source.get('clusters') or []),
            )
            if unknown_clusters:
                flash(
                    '所选 DEG 任务中不存在以下 cluster：' + '、'.join(unknown_clusters[:5]),
                    'danger',
                )
                return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
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
        if selected_sc_deg_source and selected_sc_deg_source['output_adata_path']:
            if os.path.abspath(input_path) != selected_sc_deg_source['output_adata_path']:
                flash('富集输入 AnnData 必须与所选 DEG 任务的输入一致，避免将通路结果绑定到不同数据版本。', 'danger')
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
        task = AnalysisTask(project_id=pid, module_name=module_name)
        # Internal artifacts and downstream source tables are keyed by task
        # ID.  This keeps concurrent/rerun outputs distinct without exposing
        # a growing set of mutable canonical filenames.
        params['_analysis_id'] = task.id
        task.params_json = json.dumps(params, ensure_ascii=False)
        task.save()
        p.status = 'processing'
        p.save()
        from worker import submit_task
        submit_task(task.id, pid, module_name, params,
                   Config.project_dir(pid), input_path)
        return redirect(url_for('results.task_detail', pid=pid, task_id=task.id))

    enrichment_deg_sources = []
    heatmap_deg_sources = []
    sc_enrichment_deg_sources = []
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
    if module_name == 'sc_cell_go':
        sc_enrichment_deg_sources = _sc_deg_enrichment_sources(pid)
        if len(sc_enrichment_deg_sources) == 1:
            for param in schema:
                if param['key'] == 'deg_source_task_id':
                    param['default'] = sc_enrichment_deg_sources[0]['task_id']

    sidebar_modules = BULK_MODULE_LIST if is_bulk else SC_MODULE_LIST
    return render_template('analysis_select.html', project=p, module=mod_info,
                          schema=schema, completed_tasks=completed_tasks,
                          uploaded_h5ad=uploaded_files, input_options=input_options,
                          recommended_input=recommended_input,
                          enrichment_deg_sources=enrichment_deg_sources,
                          heatmap_deg_sources=heatmap_deg_sources,
                          sc_enrichment_deg_sources=sc_enrichment_deg_sources,
                          base_grn_files=_list_base_grn_files(pid),
                          all_modules=sidebar_modules,
                          module_display_map=MODULE_DISPLAY_MAP)
