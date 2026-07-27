# modules/ai_tools.py
"""AI 工具函数 — 平台分析功能的接口"""
import os
import json
from config import Config


BULK_UPLOAD_EXTENSIONS = ('.h5ad', '.csv', '.txt', '.tsv', '.xlsx', '.xls')


def _latest_bulk_upload(project_id):
    """Return the newest supported Bulk upload, or None."""
    uploads_dir = Config.uploads_dir(project_id)
    if not os.path.isdir(uploads_dir):
        return None
    files = [
        os.path.join(uploads_dir, name)
        for name in os.listdir(uploads_dir)
        if name.lower().endswith(BULK_UPLOAD_EXTENSIONS)
        and os.path.isfile(os.path.join(uploads_dir, name))
    ]
    return max(files, key=os.path.getmtime) if files else None


def resolve_current_analysis_input(project_id, module_name=''):
    """Resolve the current analysis input, including raw Bulk uploads."""
    ctx = resolve_current_adata_path(project_id, allow_upload=True)
    if ctx['path']:
        return ctx
    if module_name.startswith('bulk_'):
        bulk_path = _latest_bulk_upload(project_id)
        if bulk_path:
            return {
                'source': 'bulk_upload',
                'path': bulk_path,
                'accepted_branch_id': None,
                'task_id': None,
                'warnings': [],
            }
    return ctx


def execute_tool(name, args, project_id):
    """执行指定的工具函数"""
    if name == "run_analysis":
        return _run_analysis(args, project_id)
    elif name == "run_pipeline":
        return _run_pipeline(args, project_id)
    elif name == "get_project_status":
        return _get_project_status(project_id)
    elif name == "get_task_results":
        return _get_task_results(args.get("task_id", ""))
    elif name == "list_modules":
        return _list_modules(args.get("pipeline_type", "all"))
    elif name == "recommend_analysis_config":
        return _recommend_analysis_config(args, project_id)
    elif name == "inspect_analysis_state":
        return _inspect_analysis_state(project_id)
    elif name == "inspect_adata":
        return _inspect_adata(args, project_id)
    elif name == "get_cluster_summary":
        return _get_cluster_summary(args, project_id)
    elif name == "score_cell_type_signature":
        return _score_cell_type_signature(args, project_id)
    elif name == "list_builtin_markers":
        return _list_builtin_markers()
    elif name == "propose_parameter_sweep":
        return _propose_parameter_sweep(args, project_id)
    elif name == "run_parameter_sweep":
        return _run_parameter_sweep(args, project_id)
    elif name == "start_goal_agent":
        return _start_goal_agent(args, project_id)
    elif name == "continue_goal_agent":
        return _continue_goal_agent(args, project_id)
    return {"error": f"未知工具: {name}"}


def _run_pipeline(args, project_id):
    """Submit an AI-requested sequential pipeline as one background run."""
    from models import PipelineRun, Project
    from worker import submit_pipeline_run
    from modules import MODULE_REGISTRY, SC_MODULE_NAMES, BULK_MODULE_NAMES, validate_pipeline_order
    from modules.schemas import PARAM_SCHEMAS, filter_active_params

    project = Project.get_by_id(project_id)
    if not project:
        return {"error": "项目不存在"}

    analysis_type = str(args.get("analysis_type", "")).strip()
    modules = args.get("modules", [])
    raw_params = args.get("params", {}) or {}
    if analysis_type not in {"sc", "bulk"}:
        return {"error": "analysis_type 必须是 sc 或 bulk"}
    if not isinstance(modules, list) or not modules:
        return {"error": "模块列表不能为空"}
    if not isinstance(raw_params, dict):
        return {"error": "params 必须按模块名组成对象"}

    permitted_modules = SC_MODULE_NAMES if analysis_type == "sc" else BULK_MODULE_NAMES
    if len(set(modules)) != len(modules):
        return {"error": "模块列表不能包含重复模块"}
    for module_name in modules:
        if module_name not in MODULE_REGISTRY:
            return {"error": f"未知模块: {module_name}"}
        if module_name not in permitted_modules:
            return {"error": f"模块 {module_name} 不属于 {analysis_type} 流程"}
    is_valid, errors = validate_pipeline_order(modules)
    if not is_valid:
        return {"error": "模块顺序不满足依赖约束", "details": errors}

    input_path = str(args.get("input_path", "") or "")
    input_source = "manual"
    if not input_path:
        context = resolve_current_analysis_input(project_id, modules[0])
        input_path = context.get("path") or ""
        input_source = context.get("source", "auto")
    path_error = _validate_project_path(input_path, project_id)
    if path_error:
        return {"error": f"input_path: {path_error}"}

    params_by_module = {}
    for module_name in modules:
        supplied = raw_params.get(module_name, {})
        if not isinstance(supplied, dict):
            return {"error": f"模块 {module_name} 的参数必须是对象"}
        cleaned, error = _validate_analysis_params(module_name, supplied)
        if error:
            return {"error": f"模块 {module_name} 参数校验失败: {error}"}
        schema = PARAM_SCHEMAS.get(module_name, [])
        defaults = {}
        for field in schema:
            default = field.get('default')
            if field.get('type') == 'select' and field.get('options'):
                if default not in field['options']:
                    default = field['options'][0]
            defaults[field['key']] = default
        defaults.update(cleaned)
        params_by_module[module_name] = filter_active_params(schema, defaults)

    pipeline_name = str(args.get("name", "") or "").strip()
    if not pipeline_name:
        pipeline_name = f"AI {analysis_type.upper()} 全流程"
    project_dir = Config.project_dir(project_id)
    pipeline_run = PipelineRun(
        project_id=project_id,
        name=pipeline_name,
        analysis_type=analysis_type,
        input_path=input_path,
        modules_json=json.dumps(modules, ensure_ascii=False),
        params_json=json.dumps(params_by_module, ensure_ascii=False),
    )
    pipeline_run.save()
    submit_pipeline_run(
        run_id=pipeline_run.id,
        project_id=project_id,
        modules=modules,
        params_by_module=params_by_module,
        project_dir=project_dir,
        input_path=input_path,
    )
    return {
        "status": "submitted",
        "pipeline_run_id": pipeline_run.id,
        "modules": modules,
        "input_path": input_path,
        "input_source": input_source,
        "status_url": f"/api/pipeline-runs/{pipeline_run.id}/status",
        "pipeline_url": f"/projects/{project_id}/pipeline-runs/{pipeline_run.id}",
        "message": f"全流程已在后台提交（{len(modules)} 个模块），将按顺序自动执行。",
    }


def _run_analysis(args, project_id):
    """执行分析模块"""
    from models import AnalysisTask
    from worker import submit_task
    from modules import MODULE_REGISTRY

    module_name = args.get("module_name", "")
    if module_name not in MODULE_REGISTRY:
        return {"error": f"未知模块: {module_name}，可用: {list(MODULE_REGISTRY.keys())}"}

    project_dir = Config.project_dir(project_id)
    if not os.path.isdir(project_dir):
        return {"error": "项目不存在"}

    # 确定输入文件
    input_path = args.get("input_path", "")
    input_source = "manual"
    if not input_path:
        # 优先使用已有 h5ad；Bulk 首次分析则回退到原始表格上传文件。
        ctx = resolve_current_analysis_input(project_id, module_name)
        if ctx['path']:
            input_path = ctx['path']
            input_source = ctx['source']

    if not input_path or not os.path.exists(input_path):
        return {"error": "未找到可用的输入文件，请先上传数据"}

    # Validate input_path is within project directory
    if input_path:
        real_input = os.path.realpath(input_path)
        real_project = os.path.realpath(project_dir)
        if not real_input.startswith(real_project + os.sep) and real_input != real_project:
            return {"error": "输入文件路径必须在项目目录内"}

    # Reject symlinks and special files
    if os.path.islink(input_path):
        return {"error": "不支持符号链接文件"}
    if not os.path.isfile(input_path):
        return {"error": "输入路径不是普通文件"}

    params = args.get("params", {})
    params, err = _validate_analysis_params(module_name, params)
    if err:
        return {"error": f"参数校验失败: {err}"}

    # AI 选择样本名自动分组时，由受信任后端注入映射；
    # 不允许模型直接伪造下划线内部参数。
    auto_group_fields = ('groupby', 'group_column')
    if any(params.get(field) == '_auto_group_' for field in auto_group_fields):
        mapping = _infer_auto_group_mapping(input_path)
        if not mapping:
            return {"error": "无法从样本名生成可靠的自动分组映射"}
        params['_auto_group_mapping'] = mapping

    compatibility_error = _validate_method_compatibility(module_name, input_path, params)
    if compatibility_error:
        return {"error": f"方法与数据不兼容: {compatibility_error}"}

    # AI-triggered execution must respect the same hard experimental-design
    # boundary as the regular form.  Keep this best-effort around malformed
    # legacy test/placeholder files; module-specific validation remains the
    # fallback for inputs that cannot be profiled here.
    if module_name in {'bulk_deg', 'sc_timecourse', 'batch_correct', 'bulk_normalize'}:
        try:
            from modules.design_preflight import preflight_blockers
            from modules.io_utils import read_expression_matrix

            blockers = preflight_blockers(
                read_expression_matrix(input_path), module_name, params, input_path,
            )
        except Exception:
            blockers = []
        if blockers:
            return {"error": "分析前检查未通过: " + "；".join(blockers[:2])}

    # 创建任务
    task = AnalysisTask(
        project_id=project_id,
        module_name=module_name,
        params_json=json.dumps(params, ensure_ascii=False),
    )
    task.save()

    # 提交到线程池执行
    submit_task(
        task.id, project_id, module_name,
        params, project_dir, input_path,
    )

    return {
        "status": "submitted",
        "task_id": task.id,
        "module": module_name,
        "input_path": input_path,
        "input_source": input_source,
        "message": f"分析任务已提交：{module_name}，任务 ID: {task.id}",
    }


def _get_project_status(project_id):
    """获取项目状态"""
    from models import AnalysisTask

    project_dir = Config.project_dir(project_id)
    if not os.path.isdir(project_dir):
        return {"error": "项目不存在"}

    tasks = AnalysisTask.get_by_project(project_id, include_branches=True)

    # 已上传文件
    uploads_dir = os.path.join(project_dir, 'uploads')
    uploaded = []
    uploaded_details = []
    if os.path.isdir(uploads_dir):
        uploaded = sorted(
            f for f in os.listdir(uploads_dir)
            if f.lower().endswith(BULK_UPLOAD_EXTENSIONS)
        )
        for name in uploaded:
            path = os.path.join(uploads_dir, name)
            uploaded_details.append({
                "name": name,
                "path": path,
                "extension": os.path.splitext(name)[1].lower(),
                "size_mb": round(os.path.getsize(path) / (1024 ** 2), 2),
            })

    # 已完成任务（排除分支任务，避免候选任务污染主线列表）
    completed_tasks = []
    for t in tasks:
        if t.status == 'completed' and not t.branch_id:
            completed_tasks.append({
                "id": t.id,
                "module": t.module_name,
                "finished_at": str(t.finished_at) if t.finished_at else None,
            })

    # 分支任务单独列出
    branch_tasks = []
    for t in tasks:
        if t.branch_id:
            branch_tasks.append({
                "id": t.id,
                "module": t.module_name,
                "branch_id": t.branch_id,
                "status": t.status,
            })

    # 最新可用文件
    intermediate_dir = os.path.join(project_dir, 'intermediate')
    intermediates = []
    if os.path.isdir(intermediate_dir):
        intermediates = [f for f in os.listdir(intermediate_dir) if f.endswith('.h5ad')]

    return {
        "project_id": project_id,
        "uploaded_files": uploaded,
        "uploaded_file_details": uploaded_details,
        "completed_tasks": completed_tasks,
        "intermediate_files": intermediates,
        "n_tasks_total": len(tasks),
        "n_tasks_completed": len(completed_tasks),
        "n_branch_tasks": len(branch_tasks),
    }


def _get_task_results(task_id):
    """获取任务结果摘要"""
    from models import AnalysisTask, ResultFile

    t = AnalysisTask.get_by_id(task_id)
    if not t:
        return {"error": "任务不存在"}

    result_data = {}
    try:
        result_data = json.loads(t.result_json) if t.result_json else {}
    except Exception:
        pass

    files = ResultFile.get_by_task(task_id)
    result_files = []
    for f in files:
        result_files.append({
            "label": f.label,
            "type": f.file_type,
            "category": f.category,
        })

    return {
        "task_id": t.id,
        "module": t.module_name,
        "status": t.status,
        "result_data": result_data,
        "result_files": result_files,
    }


def _list_modules(pipeline_type="all"):
    """列出可用模块"""
    from modules import MODULE_REGISTRY

    modules = {}
    for name, cls in MODULE_REGISTRY.items():
        is_bulk = name.startswith('bulk_')
        is_sc = not is_bulk and name != 'convert_10x'

        if pipeline_type == 'bulk' and not is_bulk:
            continue
        if pipeline_type == 'sc' and not is_sc:
            continue

        modules[name] = {
            "display": getattr(cls, 'DISPLAY_NAME', name),
            "description": getattr(cls, 'DESCRIPTION', ''),
        }

    return {"modules": modules}


def _matrix_value_profile(matrix, n_obs, n_vars):
    """Summarize matrix values without densifying large sparse matrices."""
    import numpy as np
    from scipy import sparse

    total_values = max(1, int(n_obs) * int(n_vars))
    if sparse.issparse(matrix):
        observed = np.asarray(matrix.data, dtype=float)
        zero_fraction = 1.0 - float(matrix.nnz) / total_values
    else:
        flat = np.asarray(matrix, dtype=float).reshape(-1)
        zero_fraction = float(np.mean(flat == 0)) if flat.size else 0.0
        observed = flat
    if observed.size > 200000:
        step = max(1, observed.size // 200000)
        observed = observed[::step][:200000]
    finite = observed[np.isfinite(observed)]
    if finite.size:
        integer_fraction = float(np.mean(np.isclose(finite, np.round(finite))))
        nonnegative_fraction = float(np.mean(finite >= 0))
        value_min = float(np.min(finite))
        value_max = float(np.max(finite))
        positive = finite[finite > 0]
        median_positive = float(np.median(positive)) if positive.size else 0.0
    else:
        integer_fraction = 0.0
        nonnegative_fraction = 0.0
        value_min = value_max = median_positive = 0.0

    measurement_type = (
        'raw_counts'
        if nonnegative_fraction >= 0.999 and integer_fraction >= 0.995
        else 'continuous_expression'
    )
    return {
        'measurement_type': measurement_type,
        'integer_fraction': round(integer_fraction, 4),
        'nonnegative_fraction': round(nonnegative_fraction, 4),
        'zero_fraction': round(zero_fraction, 4),
        'value_min': round(value_min, 4),
        'value_max': round(value_max, 4),
        'median_positive': round(median_positive, 4),
    }


def _profile_analysis_input(input_path):
    """Build a bounded data profile for deterministic method/parameter selection."""
    import pandas as pd
    from modules.io_utils import (
        read_expression_matrix, infer_sample_group_candidates,
        rank_obs_grouping_candidates,
    )
    from modules.sc_timecourse import _time_value

    adata = read_expression_matrix(input_path)
    value_profile = _matrix_value_profile(adata.X, adata.n_obs, adata.n_vars)
    normalization = dict(adata.uns.get('normalization', {})) if hasattr(adata, 'uns') else {}
    measurement_type = (
        'log_transformed'
        if normalization.get('is_log_transformed')
        else value_profile['measurement_type']
    )
    source_measurement_type = value_profile['measurement_type']
    if 'raw' in adata.layers:
        source_measurement_type = _matrix_value_profile(
            adata.layers['raw'], adata.n_obs, adata.n_vars)['measurement_type']
    counts_layer_profile = (
        _matrix_value_profile(adata.layers['counts'], adata.n_obs, adata.n_vars)
        if 'counts' in adata.layers else None
    )

    obs_columns = [str(c) for c in adata.obs.columns[:50]]
    time_name_tokens = ('time', 'day', 'hour', 'week', 'stage', 'minute')
    technical_batch_names = {'batch', 'technical_batch', 'sequencing_batch', 'library_batch'}
    time_candidates = []
    for column in adata.obs.columns:
        if str(column).lower() in technical_batch_names:
            continue
        values = adata.obs[column].dropna().astype(str).unique().tolist()
        if not (3 <= len(values) <= 30):
            continue
        is_numeric = pd.api.types.is_numeric_dtype(adata.obs[column])
        looks_temporal = bool(values) and all(_time_value(value) is not None for value in values)
        name_is_temporal = any(token in str(column).lower() for token in time_name_tokens)
        if is_numeric or looks_temporal or name_is_temporal:
            time_candidates.append(str(column))
    # Do not auto-select donor or technical-batch labels as sample IDs.  A
    # donor may recur across timepoints and a batch may be technical; both
    # would otherwise create invalid independent-replicate assumptions.
    sample_priority = ('sample_id', 'sample', 'library_id', 'orig.ident')
    lower_to_original = {str(column).lower(): str(column) for column in adata.obs.columns}
    sample_candidates = [lower_to_original[name] for name in sample_priority
                         if name in lower_to_original]

    grouping_candidates = []
    preferred_group_columns = rank_obs_grouping_candidates(
        adata, purpose='groupby',
    )
    for column in preferred_group_columns:
        series = adata.obs[column].astype(str)
        counts = series.value_counts()
        if 2 <= len(counts) <= 50 and int(counts.min()) >= 2:
            grouping_candidates.append({
                'source': 'obs',
                'column': str(column),
                'label': str(column),
                'values': sorted(counts.index.tolist()),
                'group_sizes': {str(k): int(v) for k, v in counts.sort_index().items()},
            })

    batch_candidates = []
    for column in rank_obs_grouping_candidates(
        adata, purpose='batch_key', include_technical=True,
    ):
        series = adata.obs[column].astype(str)
        counts = series.value_counts()
        if 2 <= len(counts) <= 50 and int(counts.min()) >= 2:
            batch_candidates.append({
                'source': 'obs',
                'column': str(column),
                'label': str(column),
                'values': sorted(counts.index.tolist()),
                'group_sizes': {str(k): int(v) for k, v in counts.sort_index().items()},
            })
    if not grouping_candidates:
        for candidate in infer_sample_group_candidates(adata.obs_names.tolist()):
            grouping_candidates.append({
                'source': 'sample_names',
                'column': '_auto_group_',
                **candidate,
            })

    return {
        'n_samples_or_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'measurement_type': measurement_type,
        'source_measurement_type': source_measurement_type,
        **{key: value for key, value in value_profile.items() if key != 'measurement_type'},
        'normalization': normalization,
        'obs_columns': obs_columns,
        'sample_names': [str(name) for name in adata.obs_names[:100]],
        'grouping_candidates': grouping_candidates[:8],
        'batch_candidates': batch_candidates[:8],
        'time_candidates': time_candidates[:8],
        'sample_candidates': sample_candidates,
        'counts_layer_present': bool('counts' in adata.layers),
        'counts_layer_is_raw': bool(
            counts_layer_profile
            and counts_layer_profile['measurement_type'] == 'raw_counts'
        ),
    }


def _infer_auto_group_mapping(input_path):
    """Return the recommended sample-name mapping for a trusted project input."""
    from modules.io_utils import read_expression_matrix, infer_sample_group_candidates

    adata = read_expression_matrix(input_path)
    candidates = infer_sample_group_candidates(adata.obs_names.tolist())
    return candidates[0]['mapping'] if candidates else {}


def _validate_method_compatibility(module_name, input_path, params):
    """Reject high-risk method/data combinations before an AI-submitted task is created."""
    if module_name not in {'bulk_normalize', 'bulk_deg'}:
        return None
    try:
        profile = _profile_analysis_input(input_path)
    except Exception as exc:
        return f'无法确认输入数据类型: {exc}'

    measurement = profile['measurement_type']
    source_measurement = profile.get('source_measurement_type', measurement)
    normalization_method = str(profile.get('normalization', {}).get('method', '')).lower()
    compromised_normalization = (
        measurement == 'log_transformed'
        and source_measurement == 'continuous_expression'
        and normalization_method in {'deseq2', 'tmm', 'cpm', 'vst', 'rlog'}
    )
    method = str(params.get('method', '') or '').lower()
    count_only_normalizers = {'deseq2', 'tmm', 'cpm', 'vst', 'rlog'}
    count_only_deg = {'deseq2', 'edger', 'limma'}

    if module_name == 'bulk_normalize':
        if measurement == 'continuous_expression' and method in count_only_normalizers:
            return f'连续 FPKM/TPM 类表达值不应使用 {method}；请选择 log2'
        if measurement == 'log_transformed':
            return '输入已是 log/标准化尺度，不应重复运行 bulk_normalize'

    if module_name == 'bulk_deg':
        if compromised_normalization:
            return f'当前中间文件用 {normalization_method} 处理了连续表达值；请从原始上传文件重新运行 log2'
        if measurement in {'continuous_expression', 'log_transformed'} and method in count_only_deg:
            return f'{method} 需要原始整数 counts，当前为 {measurement}'
        if measurement == 'continuous_expression' and method == 't-test':
            return '连续表达值需先运行 bulk_normalize(method=log2)，再做 Welch t-test'
    return None


def _infer_reference_comparisons(group_values):
    """Infer treatment-vs-control comparisons while preserving suffix strata."""
    values = [str(value) for value in group_values]
    controls = [
        value for value in values
        if value.lower().split('_', 1)[0] in {'ctr', 'ctrl', 'control', 'vehicle', 'untreated'}
    ]
    if not controls:
        return []
    comparisons = []
    for value in values:
        if value in controls:
            continue
        suffix = value.split('_', 1)[1] if '_' in value else ''
        matched = next(
            (control for control in controls
             if (control.split('_', 1)[1] if '_' in control else '') == suffix),
            controls[0],
        )
        comparisons.append(f'{value}-vs-{matched}')
    return comparisons[:20]


def _recommendation_input(args, project_id, module_name):
    explicit = str(args.get('input_path', '') or '').strip()
    if explicit:
        return explicit, 'manual'
    if module_name in {'bulk_qc', 'bulk_normalize'}:
        raw_upload = _latest_bulk_upload(project_id)
        if raw_upload:
            return raw_upload, 'bulk_upload'
    ctx = resolve_current_analysis_input(project_id, module_name)
    return ctx.get('path'), ctx.get('source', 'none')


def _recommend_analysis_config(args, project_id):
    """[Read-only] Recommend method and parameters from the actual input profile."""
    from modules import MODULE_REGISTRY
    from modules.schemas import PARAM_SCHEMAS

    module_name = str(args.get('module_name', '') or '').strip()
    objective = str(args.get('objective', '') or '').strip()
    if module_name not in MODULE_REGISTRY:
        return {'error': f'未知模块: {module_name}'}

    input_path, input_source = _recommendation_input(args, project_id, module_name)
    if not input_path:
        return {'error': '未找到可用输入文件'}
    path_error = _validate_project_path(input_path, project_id)
    if path_error:
        return {'error': path_error}

    try:
        profile = _profile_analysis_input(input_path)
    except Exception as exc:
        return {'error': f'数据画像读取失败: {exc}'}

    schema = PARAM_SCHEMAS.get(module_name, [])
    params = {entry['key']: entry.get('default') for entry in schema}
    rationale = []
    warnings = []
    prerequisites = []
    alternatives = []
    should_run = True

    def recommend(key, value, reason):
        if any(entry['key'] == key for entry in schema):
            params[key] = value
            rationale.append({'parameter': key, 'value': value, 'reason': reason})

    groups = profile.get('grouping_candidates') or []
    preferred_group = groups[0] if groups else None
    group_column = preferred_group.get('column') if preferred_group else ''
    group_sizes = preferred_group.get('group_sizes', {}) if preferred_group else {}
    min_group_size = min(group_sizes.values()) if group_sizes else 0
    measurement = profile['measurement_type']
    source_measurement = profile.get('source_measurement_type', measurement)
    normalization_method = str(profile.get('normalization', {}).get('method', '')).lower()
    compromised_normalization = (
        measurement == 'log_transformed'
        and source_measurement == 'continuous_expression'
        and normalization_method in {'deseq2', 'tmm', 'cpm', 'vst', 'rlog'}
    )
    n_obs = profile['n_samples_or_cells']
    n_genes = profile['n_genes']

    if module_name == 'bulk_qc':
        if measurement == 'continuous_expression':
            recommend('filter_strategy', 'custom', '连续表达值不适用 raw read-count 硬阈值')
            recommend('min_counts', 0, '避免用 count 文库大小阈值错误过滤 FPKM/TPM 样本')
            recommend('min_genes', max(100, min(5000, int(n_genes * 0.2))), '仅保留基础检测基因数质控')
        if group_column:
            recommend('group_column', group_column, '使用检测到的样本分组生成分组 QC')
        recommend('detect_outliers', True, '保留 PCA 离群告警，不自动删除样本')

    elif module_name == 'bulk_normalize':
        if measurement == 'raw_counts':
            recommend('method', 'deseq2', '近似整数 raw counts 适合中位比率归一化')
            alternatives.append({'method': 'tmm', 'when': '样本间组成偏差明显时'})
        elif measurement == 'continuous_expression':
            recommend('method', 'log2', '连续 FPKM/TPM 类表达值仅做 log2(x+1)，保留真实分布差异')
            alternatives.append({'method': 'log2_quantile', 'when': '仅在确认所有样本应具有相同分布时'})
        else:
            should_run = False
            warnings.append('输入已是 log/标准化尺度，不建议重复标准化')
        recommend('min_expr_samples', min_group_size or min(3, n_obs), '使用最小生物学分组样本数过滤低表达基因')

    elif module_name == 'bulk_deg':
        if group_column:
            recommend('groupby', group_column, '使用数据中最可解释的分组候选')
            comparisons = _infer_reference_comparisons(preferred_group.get('values', []))
            if comparisons:
                recommend('comparisons', ';'.join(comparisons), '在相同层次内将各处理与对照比较')
        else:
            warnings.append('未检测到可靠分组，执行前需指定 groupby/group1/group2')
            should_run = False
        if measurement == 'raw_counts':
            recommend('method', 'deseq2', '原始整数 count 适用负二项分布模型')
            alternatives.append({'method': 'edger', 'when': '小样本或组成偏差明显，并使用 raw counts'})
        elif measurement in {'continuous_expression', 'log_transformed'}:
            recommend('method', 't-test', '连续表达值在 log2 尺度上使用 Welch t-test')
            if measurement == 'continuous_expression':
                prerequisites.append({'module': 'bulk_normalize', 'params': {'method': 'log2'}, 'reason': 'DEG 前先进行 log2(x+1)'})
                should_run = False
        if compromised_normalization:
            should_run = False
            prerequisites.append({'module': 'bulk_normalize', 'params': {'method': 'log2'}, 'reason': '原始连续表达值曾被 count-only 方法处理，需从原始上传重新 log2'})
            warnings.append(f'当前中间文件用 {normalization_method} 处理了连续表达值，不建议直接用于 DEG')
        objective_lower = objective.lower()
        strict = any(word in objective_lower for word in ('严格', '严谨', '验证', 'strict'))
        exploratory = any(word in objective_lower for word in ('探索', '敏感', '召回', 'explor'))
        fc_threshold = 2.0 if strict else (1.5 if exploratory or (min_group_size and min_group_size <= 3) else 2.0)
        fdr_threshold = 0.01 if strict else (0.1 if exploratory else 0.05)
        recommend('fc_threshold', fc_threshold, '根据用户目标和最小分组样本数设置效应量阈值')
        recommend('pval_threshold', fdr_threshold, '验证目标使用严格 FDR，探索目标允许 FDR 0.1 但必须标注为候选')
        recommend('padj_method', 'fdr_bh', '使用 Benjamini-Hochberg 控制 FDR')
        recommend('top_n', 50, '保留足够候选基因用于图表和复核')

    elif module_name == 'bulk_pca':
        recommend('n_comps', max(2, min(10, n_obs - 1)), '主成分数不超过样本数-1')
        if group_column:
            recommend('color_by', group_column if group_column != '_auto_group_' else '_auto_group', '按推荐分组着色检查样本结构')
        if measurement == 'continuous_expression':
            prerequisites.append({'module': 'bulk_normalize', 'params': {'method': 'log2'}, 'reason': 'PCA 前先稳定连续表达值方差'})
            should_run = False
        elif compromised_normalization:
            prerequisites.append({'module': 'bulk_normalize', 'params': {'method': 'log2'}, 'reason': '从原始连续表达文件重新 log2 后再做 PCA'})
            warnings.append(f'当前中间文件用 {normalization_method} 处理了连续表达值')
            should_run = False

    elif module_name == 'bulk_heatmap':
        has_deg = bool([name for name in os.listdir(Config.results_dir(project_id)) if name.startswith('bulk_deg_results')]) if os.path.isdir(Config.results_dir(project_id)) else False
        recommend('gene_import_source', 'deg' if has_deg else 'top_var', '有 DEG 结果时优先显示差异基因，否则显示高变基因')
        recommend('top_n', 50, '在可读性和信息量之间取平衡')
        if group_column:
            recommend('groupby', group_column if group_column != '_auto_group_' else '_auto_group', '显示样本分组注释')
        if measurement == 'continuous_expression':
            prerequisites.append({'module': 'bulk_normalize', 'params': {'method': 'log2'}, 'reason': '热图前先对连续表达值做 log2(x+1)'})
            should_run = False
        elif compromised_normalization:
            prerequisites.append({'module': 'bulk_normalize', 'params': {'method': 'log2'}, 'reason': '从原始连续表达文件重新 log2 后再生成热图'})
            warnings.append(f'当前中间文件用 {normalization_method} 处理了连续表达值')
            should_run = False

    elif module_name == 'bulk_enrichment':
        recommend('method', 'GSEA', '当单基因 DEG 较少时，全基因排序 GSEA 比 ORA 稳健')
        recommend('database', 'GO_BP', '先检查生物学过程层面的整体变化')
        recommend('pvalue_cutoff', 0.05, '使用标准富集 FDR 阈值')
        recommend('top_n', 20, '展示主要通路并控制图表长度')
        results_dir = Config.results_dir(project_id)
        deg_files = sorted(
            [os.path.join(results_dir, name) for name in os.listdir(results_dir)
             if name.startswith('bulk_deg_results') and name.endswith('.csv')],
            key=os.path.getmtime,
            reverse=True,
        ) if os.path.isdir(results_dir) else []
        if deg_files:
            recommend('input_source', deg_files[0], '使用最新的差异分析全基因排序表')
        else:
            prerequisites.append({'module': 'bulk_deg', 'reason': 'GSEA 需要含 gene/log2FC 的差异分析结果'})
            should_run = False

    elif module_name == 'bulk_timecourse':
        time_columns = []
        for candidate in groups:
            if candidate.get('source') == 'obs' and candidate.get('column', '').lower() in {'time', 'minute', 'hour', 'day'}:
                time_columns.append(candidate['column'])
        if time_columns:
            recommend('time_column', time_columns[0], '使用检测到的时间列')
        else:
            should_run = False
            warnings.append('未检测到时间列，不应执行时序分析')

    elif module_name == 'sc_timecourse':
        time_candidates = profile.get('time_candidates', [])
        sample_candidates = profile.get('sample_candidates', [])
        obs_columns = set(profile.get('obs_columns', []))
        celltype_key = next((key for key in ('celltype', 'cell_type', 'annotation', 'leiden')
                             if key in obs_columns), '')
        condition_key = next((key for key in ('condition', 'treatment', 'group', 'treatment_group')
                              if key in obs_columns), '')
        if not time_candidates:
            should_run = False
            warnings.append('未检测到至少 3 个真实时间点列；不能用伪时间替代采样时间。')
        else:
            recommend('timepoint_key', time_candidates[0], '使用检测到的真实采样时间列')
        if not sample_candidates:
            should_run = False
            warnings.append('未检测到 sample_id/生物学重复列；请先补充元数据，避免把细胞当作独立重复。')
        else:
            recommend('sample_key', sample_candidates[0], '以独立生物学样本作为组成与伪 bulk 的统计单位')
            recommend('min_replicates_per_timepoint', 2, '每时间点至少 2 个样本才启用样本级 Kruskal 筛选')
        if celltype_key:
            recommend('celltype_key', celltype_key, '优先使用已注释细胞类型；无注释时以 leiden 作为探索性分组')
        else:
            warnings.append('未检测到 celltype/annotation/leiden 列，需先完成聚类或注释。')
            should_run = False
        if condition_key:
            recommend('condition_key', condition_key, '按条件分层展示时间趋势，避免把处理效应混入单一曲线')
        if not profile.get('counts_layer_is_raw'):
            warnings.append("未检测到非负近似整数的原始 counts 层；可进行样本级组成检验，但基因动态仅为描述性均值，建议从 QC 输出重新开始。")
        if time_candidates and time_candidates[0] == 'batch':
            warnings.append('时间点当前位于 batch 列；不要再将此列用于 Harmony/BBKNN/scVI 校正，否则可能移除真实时间信号。')

    elif module_name == 'annotation':
        # When the user's objective names an organoid context, expose the
        # corresponding tissue panel instead of silently leaving the broad
        # Universal panel selected. Maturity is evaluated separately from the
        # observed expression modules and any time metadata.
        objective_text = objective.lower()
        organoid_aliases = (
            ('intestinal', ('肠道', '肠', 'intestinal', 'gut')),
            ('cerebral', ('脑', '大脑', 'cerebral', 'brain', 'neural')),
            ('kidney', ('肾', 'kidney', 'renal')),
            ('liver', ('肝', 'liver', 'hepatic')),
            ('lung', ('肺', 'lung', 'pulmonary')),
            ('pancreatic', ('胰', 'pancrea', 'pancreatic')),
            ('cardiac', ('心脏', '心肌', 'cardiac', 'heart')),
        )
        detected_organoid = next(
            (organ_type for organ_type, aliases in organoid_aliases
             if any(alias in objective_text for alias in aliases)),
            None,
        )
        if detected_organoid:
            recommend('marker_set', 'Organoid', '目标描述包含类器官组织，使用对应的组织/发育 marker panel')
            recommend('organoid_type', detected_organoid, '从目标描述识别类器官类型，避免误用 Universal panel')
            recommend('method', 'multi_evidence', '类器官细胞状态和阶段差异较大，保留 Marker、一致性和复核证据')
            warnings.append('类器官 marker 仅作为第一轮候选注释；成熟度改由表达模块和样本时间元数据独立评估。')
        if any(token in objective_text for token in ('celltypist', 'model zoo', '参考模型', '参考注释')):
            recommend('use_celltypist_reference', True, '用户明确要求 CellTypist；作为本地参考证据运行，不覆盖 Marker 最终标签')
            reference_models = {
                'intestinal': 'Cells_Intestinal_Tract.pkl',
                'cerebral': 'Developing_Human_Brain.pkl',
                'lung': 'Cells_Fetal_Lung.pkl',
                'liver': 'Healthy_Human_Liver.pkl',
                'pancreatic': 'Fetal_Human_Pancreas.pkl',
                'cardiac': 'Healthy_Adult_Heart.pkl',
                'kidney': 'Developing_Human_Organs.pkl',
            }
            if detected_organoid:
                recommend(
                    'celltypist_model',
                    reference_models.get(detected_organoid, 'Pan_Fetal_Human.pkl'),
                    '选择与目标类器官组织最接近的人类参考模型；肾脏使用广义发育器官模型作为近似参考',
                )
            else:
                recommend('celltypist_model', 'Immune_All_Low.pkl', '未识别组织时仅使用通用免疫模型作为可选交叉证据')
            recommend('celltypist_mode', 'prob match', '保留 Unassigned/低置信度结果，避免将类器官细胞强行映射到参考标签')

    elif module_name == 'normalize':
        recommend('method', 'log1p', '标准 UMI 单细胞数据默认使用 log1p CPM')
        recommend('target_sum', 10000, '使用 Scanpy 常用的每细胞总数')

    elif module_name == 'hvg':
        recommend('n_top_genes', min(3000, max(1000, int(n_genes * 0.1))), '根据数据基因总数调整 HVG 规模')
        recommend('hvg_flavor', 'seurat_v3', '使用方差稳定的常用 HVG 方法')

    elif module_name == 'dimred':
        recommend('n_comps', 30 if n_obs < 5000 else 50, '根据细胞规模设置 PCA 维度')
        recommend('umap_n_neighbors', 15 if n_obs < 20000 else 30, '大数据增大邻居数以保留全局结构')

    elif module_name == 'clustering':
        rare_target = any(word in objective.lower() for word in ('稀有', '少数', 'rare'))
        resolutions = ('0.8,1.2,1.5' if rare_target else
                       ('0.4,0.6,0.8' if n_obs < 5000 else ('0.6,0.8,1.0' if n_obs < 50000 else '0.8,1.0,1.2')))
        recommend('resolutions', resolutions, '根据细胞数生成可对比的多分辨率候选')
        recommend('n_neighbors', 10 if rare_target else (15 if n_obs < 20000 else 30), '稀有细胞目标强调局部结构，否则按数据规模平衡局部与全局结构')
        recommend('auto_select_resolution', True, '运行多分辨率评分并保留人工复核')

    elif module_name == 'batch_correct':
        batch_candidate = (profile.get('batch_candidates') or [None])[0]
        if not batch_candidate:
            should_run = False
            warnings.append('未检测到可靠批次列，不应盲目进行批次校正')
        else:
            recommend('batch_key', batch_candidate['column'], '使用检测到的批次列')
            recommend('method', 'harmony', '作为稳健、快速的 PCA 空间批次校正基线')
            alternatives.extend([
                {'method': 'bbknn', 'when': '目标是构建批次平衡邻居图'},
                {'method': 'scvi', 'when': '数据量大且需要生成模型潜在空间'},
            ])

    elif module_name == 'deg':
        recommend('method', 'wilcoxon', '单细胞表达稀疏且非正态，优先使用 Wilcoxon')
        if group_column:
            recommend('groupby', group_column, '使用已有聚类或注释列')

    cleaned_params, validation_error = _validate_analysis_params(module_name, params)
    if validation_error:
        return {'error': validation_error}

    return {
        'module_name': module_name,
        'objective': objective,
        'input_path': input_path,
        'input_source': input_source,
        'data_profile': profile,
        'should_run': should_run,
        'recommended_params': cleaned_params,
        'rationale': rationale,
        'alternatives': alternatives,
        'prerequisites': prerequisites,
        'warnings': warnings,
        'needs_confirmation': False,
        'next_action': '向用户展示方法、参数、理由和风险；用户确认后再调用 run_analysis。',
    }


def _validate_analysis_params(module_name, params):
    """校验 AI 传入的分析参数，移除未知键，返回 (cleaned_params, error_msg)。"""
    from modules.schemas import PARAM_SCHEMAS, filter_active_params

    schema_list = PARAM_SCHEMAS.get(module_name, [])
    valid_keys = {s['key'] for s in schema_list}

    if not params:
        return {}, None

    cleaned = {}
    for key, value in params.items():
        if key.startswith('_'):
            continue  # 内置参数不允许 AI 设置
        if key not in valid_keys:
            continue  # 静默移除未知参数
        schema_entry = next((s for s in schema_list if s['key'] == key), None)
        if schema_entry:
            expected_type = schema_entry.get('type', 'text')
            if expected_type == 'number':
                try:
                    fval = float(value)
                    cleaned[key] = int(fval) if fval == int(fval) else fval
                except (ValueError, TypeError):
                    return None, f"参数 '{key}' 应为数字，收到: {value}"
            elif expected_type == 'checkbox':
                if isinstance(value, str):
                    cleaned[key] = value.lower() in ('true', '1', 'yes', 'on')
                else:
                    cleaned[key] = bool(value)
            else:
                cleaned[key] = str(value)
        else:
            cleaned[key] = value

    return filter_active_params(schema_list, cleaned), None


# ============================================================
# 阶段 2：只读分析状态工具
# ============================================================

def _validate_project_path(path, project_id):
    """验证路径安全：在项目目录内、不是符号链接、是普通文件。返回 error 或 None."""
    project_dir = Config.project_dir(project_id)
    if not os.path.isdir(project_dir):
        return "项目不存在"
    if not path:
        return "缺少路径参数"
    if os.path.islink(path):
        return "不支持符号链接文件"
    try:
        real_path = os.path.realpath(path)
        real_project = os.path.realpath(project_dir)
    except (OSError, ValueError) as e:
        return f"路径解析失败: {e}"
    if not real_path.startswith(real_project + os.sep) and real_path != real_project:
        return f"路径不在项目目录内: {path}"
    if not os.path.isfile(path):
        return "路径不是普通文件"
    return None


def resolve_current_adata_path(project_id, allow_upload=True):
    """
    统一解析项目当前 h5ad 路径。所有 AI 工具和普通分析入口共用此函数。

    优先级：
    1. accepted branch 的 output_adata_path（需 branch 存在、accepted=1、deleted=0）
    2. 最新主线 completed task（branch_id IS NULL）
    3. intermediate/ 最新 .h5ad
    4. uploads/ 最新 .h5ad（仅当 allow_upload=True）

    返回：dict with keys: source, path, accepted_branch_id, task_id, warnings
    """
    from models import Project, AnalysisTask

    result = {
        'source': 'none',
        'path': None,
        'accepted_branch_id': None,
        'task_id': None,
        'warnings': [],
    }

    # 1. 检查 accepted branch
    project = Project.get_by_id(project_id)
    if project:
        try:
            import json as _json
            meta = _json.loads(project.metadata_json or '{}')
            accepted_path = meta.get('current_adata_path', '')
            accepted_bid = meta.get('accepted_branch_id', '')
            if accepted_path and accepted_bid:
                from database import get_conn
                conn = get_conn()
                try:
                    branch_row = conn.execute(
                        "SELECT id FROM analysis_branches WHERE id=? AND accepted=1 AND deleted=0",
                        (accepted_bid,)
                    ).fetchone()
                finally:
                    conn.close()
                if branch_row and os.path.isfile(accepted_path):
                    result['source'] = 'accepted_branch'
                    result['path'] = accepted_path
                    result['accepted_branch_id'] = accepted_bid
                    return result
                elif not branch_row:
                    result['warnings'].append(f'已采纳 branch {accepted_bid} 已不存在或被删除')
                elif not os.path.isfile(accepted_path):
                    result['warnings'].append(f'已采纳 branch 输出文件不存在: {accepted_path}')
        except Exception as e:
            result['warnings'].append(f'读取项目 metadata 失败: {e}')

    # 2. 最新主线任务
    tasks = AnalysisTask.get_by_project(project_id)
    for t in tasks:
        if t.status == 'completed' and t.output_adata_path and not t.branch_id:
            if os.path.isfile(t.output_adata_path):
                result['source'] = 'main_task'
                result['path'] = t.output_adata_path
                result['task_id'] = t.id
                return result

    # 3. intermediate/
    project_dir = Config.project_dir(project_id)
    intermediate_dir = os.path.join(project_dir, 'intermediate')
    if os.path.isdir(intermediate_dir):
        h5ads = sorted(
            [f for f in os.listdir(intermediate_dir) if f.endswith('.h5ad')],
            key=lambda f: os.path.getmtime(os.path.join(intermediate_dir, f)),
            reverse=True,
        )
        if h5ads:
            result['source'] = 'intermediate'
            result['path'] = os.path.join(intermediate_dir, h5ads[0])
            return result

    # 4. uploads/
    if allow_upload:
        uploads_dir = os.path.join(project_dir, 'uploads')
        if os.path.isdir(uploads_dir):
            h5ads = sorted(
                [f for f in os.listdir(uploads_dir) if f.endswith('.h5ad')],
                key=lambda f: os.path.getmtime(os.path.join(uploads_dir, f)),
                reverse=True,
            )
            if h5ads:
                result['source'] = 'upload'
                result['path'] = os.path.join(uploads_dir, h5ads[0])
                return result

    return result


def _find_latest_adata(project_id):
    """查找项目中最新的 h5ad 文件。优先级：accepted branch > 主线任务 > intermediate > uploads."""
    from models import AnalysisTask

    # 1. 优先检查已采纳 branch （从项目 metadata 中读取）
    from models import Project
    project = Project.get_by_id(project_id)
    if project:
        latest = project.get_latest_adata_path()
        if latest and os.path.isfile(latest):
            return latest

    # 2. 回退到 intermediate 目录
    project_dir = Config.project_dir(project_id)
    intermediate_dir = os.path.join(project_dir, 'intermediate')
    if os.path.isdir(intermediate_dir):
        h5ads = sorted(
            [f for f in os.listdir(intermediate_dir) if f.endswith('.h5ad')],
            key=lambda f: os.path.getmtime(os.path.join(intermediate_dir, f)),
            reverse=True,
        )
        if h5ads:
            return os.path.join(intermediate_dir, h5ads[0])

    # 3. 回退到 uploads 目录
    uploads_dir = os.path.join(project_dir, 'uploads')
    if os.path.isdir(uploads_dir):
        h5ads = sorted(
            [f for f in os.listdir(uploads_dir) if f.endswith('.h5ad')],
            key=lambda f: os.path.getmtime(os.path.join(uploads_dir, f)),
            reverse=True,
        )
        if h5ads:
            return os.path.join(uploads_dir, h5ads[0])

    return None


def _inspect_analysis_state(project_id):
    """
    [只读] 检查项目当前分析状态。
    返回最新 h5ad 路径、已完成任务、可用聚类键/嵌入/注释等信息。
    """
    from models import AnalysisTask

    project_dir = Config.project_dir(project_id)
    if not os.path.isdir(project_dir):
        return {"error": "项目不存在"}

    tasks = AnalysisTask.get_by_project(project_id, include_branches=True)
    completed_tasks = []
    branch_tasks = []
    for t in tasks:
        if t.status == 'completed' and not t.branch_id:
            completed_tasks.append({
                "id": t.id,
                "module": t.module_name,
                "finished_at": str(t.finished_at) if t.finished_at else None,
            })
        elif t.branch_id:
            branch_tasks.append({
                "id": t.id,
                "module": t.module_name,
                "branch_id": t.branch_id,
                "status": t.status,
            })

    latest_adata_path = _find_latest_adata(project_id)

    # 读取 AnnData 获取元数据
    cluster_keys = []
    embedding_keys = []
    annotation_keys = []
    obs_columns = []
    n_obs = 0
    n_vars = 0
    warnings = []

    if latest_adata_path and os.path.isfile(latest_adata_path):
        try:
            import scanpy as sc
            adata = sc.read(latest_adata_path, backed='r')
            n_obs = adata.n_obs
            n_vars = adata.n_vars
            obs_columns = list(adata.obs.columns)

            # 识别聚类键
            for col in obs_columns:
                if col.startswith('leiden') or col.startswith('louvain'):
                    cluster_keys.append(col)
            # 识别注释列
            for col in ['celltype', 'cell_type', 'annotation', 'CellType']:
                if col in obs_columns:
                    annotation_keys.append(col)
            # 识别嵌入
            if hasattr(adata, 'obsm') and adata.obsm is not None:
                for key in adata.obsm.keys():
                    embedding_keys.append(key)
        except Exception as e:
            warnings.append(f"无法读取 h5ad 元数据: {e}")

    result = {
        "project_id": project_id,
        "latest_adata": latest_adata_path,
        "completed_tasks": completed_tasks,
        "branch_tasks": branch_tasks,
        "n_branch_tasks": len(branch_tasks),
        "available_cluster_keys": cluster_keys,
        "available_embeddings": embedding_keys,
        "available_annotations": annotation_keys,
        "obs_columns": obs_columns[:50],  # 最多 50 列
        "n_obs": n_obs,
        "n_vars": n_vars,
        "warnings": warnings,
    }
    return result


def _inspect_adata(args, project_id):
    """
    [只读] 检查指定 AnnData 或 Bulk 表达矩阵的详细结构。
    输入：{"adata_path": "..."}（兼容旧参数名）
    """
    adata_path = args.get("adata_path", "") or args.get("input_path", "")
    if not adata_path:
        ctx = resolve_current_analysis_input(project_id, 'bulk_qc')
        if ctx['path']:
            adata_path = ctx['path']
        else:
            return {"error": "未找到当前数据文件，请先上传数据或指定路径"}
    err = _validate_project_path(adata_path, project_id)
    if err:
        return {"error": err}

    try:
        if adata_path.lower().endswith('.h5ad'):
            import scanpy as sc
            adata = sc.read(adata_path, backed='r')
            data_type = 'anndata'
        else:
            from modules.io_utils import read_expression_matrix
            adata = read_expression_matrix(adata_path)
            data_type = 'bulk_expression_matrix'

        # 检测可能的聚类键
        cluster_keys = []
        for col in adata.obs.columns:
            if col.startswith('leiden') or col.startswith('louvain'):
                cluster_keys.append(col)

        # 嵌入键
        obsm_keys = list(adata.obsm.keys()) if hasattr(adata, 'obsm') and adata.obsm is not None else []

        # 层
        layers = list(adata.layers.keys()) if hasattr(adata, 'layers') and adata.layers is not None else []

        # 前5行预览
        preview = []
        try:
            obs_head = adata.obs.head(5)
            for idx, row in obs_head.iterrows():
                preview.append({col: str(row[col]) for col in obs_head.columns[:10]})
        except Exception:
            pass

        return {
            "data_type": data_type,
            "input_file": os.path.basename(adata_path),
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "sample_names": [str(name) for name in adata.obs_names[:100]],
            "gene_names_preview": [str(name) for name in adata.var_names[:20]],
            "obs_columns": list(adata.obs.columns)[:50],
            "var_columns": list(adata.var.columns)[:20],
            "layers": layers,
            "obsm": obsm_keys,
            "cluster_keys": cluster_keys,
            "embedding_keys": obsm_keys,
            "sample_preview": preview,
        }
    except Exception as e:
        return {"error": f"无法读取表达数据: {e}"}


def _get_cluster_summary(args, project_id):
    """
    [只读] 获取指定聚类键的每个 cluster 概况。
    输入：{"adata_path": "...", "cluster_key": "leiden"}
    """
    adata_path = args.get("adata_path", "")
    if not adata_path:
        ctx = resolve_current_adata_path(project_id, allow_upload=True)
        if ctx['path']:
            adata_path = ctx['path']
        else:
            return {"error": "未找到当前 h5ad 文件"}
    cluster_key = args.get("cluster_key", "leiden")

    err = _validate_project_path(adata_path, project_id)
    if err:
        return {"error": err}

    try:
        import scanpy as sc
        import numpy as np

        adata = sc.read(adata_path, backed='r')

        if cluster_key not in adata.obs.columns:
            return {
                "error": f"cluster_key '{cluster_key}' 不在 obs 列中",
                "available_columns": list(adata.obs.columns),
            }
        from modules.io_utils import obs_grouping_info
        grouping = obs_grouping_info(
            adata, cluster_key, max_categories=50,
            max_numeric_categories=20, require_multiple=False,
        )
        if not grouping['valid']:
            return {
                "error": f"cluster_key '{cluster_key}' 不是有效的分类聚类列：{grouping['reason']}",
                "cluster_key": cluster_key,
            }

        clusters = sorted(adata.obs[cluster_key].unique().astype(str))
        total_cells = adata.n_obs

        # 检测 batch/sample 列
        batch_col = None
        for col in ['batch', 'sample', 'sample_id', 'orig.ident']:
            if col in adata.obs.columns:
                batch_col = col
                break

        # 检测 UMAP 坐标
        umap_key = None
        for key in ['X_umap', 'UMAP']:
            if hasattr(adata, 'obsm') and adata.obsm is not None and key in adata.obsm:
                umap_key = key
                break

        result_clusters = []
        for cl in clusters:
            mask = adata.obs[cluster_key].astype(str) == cl
            n_cells = int(mask.sum())

            info = {"cluster": cl, "n_cells": n_cells, "pct_cells": round(n_cells / total_cells, 4)}

            # batch/sample 分布
            if batch_col and batch_col in adata.obs.columns:
                try:
                    counts = adata[mask].obs[batch_col].value_counts().to_dict()
                    info["batch_distribution"] = {str(k): int(v) for k, v in counts.items()}
                except Exception:
                    info["batch_distribution"] = {}

            # UMAP 中心
            if umap_key:
                try:
                    coords = adata[mask].obsm[umap_key]
                    info["umap_center"] = [round(float(np.mean(coords[:, i])), 4) for i in range(min(2, coords.shape[1]))]
                except Exception:
                    info["umap_center"] = []

            # QC 均值
            qc_means = {}
            for qc_col in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'percent_mito']:
                if qc_col in adata.obs.columns:
                    try:
                        qc_means[qc_col] = round(float(adata[mask].obs[qc_col].mean()), 4)
                    except Exception:
                        pass
            info["qc_means"] = qc_means

            result_clusters.append(info)

        return {
            "cluster_key": cluster_key,
            "total_cells": total_cells,
            "n_clusters": len(result_clusters),
            "clusters": result_clusters,
        }
    except Exception as e:
        return {"error": f"获取 cluster 摘要失败: {e}"}


# ============================================================
# 阶段 3：细胞类型 marker 打分工具
# ============================================================

def _score_cell_type_signature(args, project_id):
    """
    [只读] 对指定 AnnData 的聚类进行细胞类型签名评分。
    输入：
        {
            "adata_path": "...",
            "cluster_key": "leiden",
            "target_cell_type": "microglia",
            "positive_markers": [...],   // 可选
            "negative_markers": [...]    // 可选
        }
    """
    adata_path = args.get("adata_path", "")
    if not adata_path:
        adata_path = _find_latest_adata(project_id)
        if not adata_path:
            return {"error": "未找到可用 h5ad 文件，请先运行分析或指定路径"}

    err = _validate_project_path(adata_path, project_id)
    if err:
        return {"error": err}

    cluster_key = args.get("cluster_key", "leiden")
    target_cell_type = args.get("target_cell_type", "")
    positive_markers = args.get("positive_markers", [])
    negative_markers = args.get("negative_markers", [])

    if not target_cell_type:
        return {"error": "缺少 target_cell_type 参数"}

    from modules.evaluators.sc_cluster import score_cell_type_signature as do_score

    return do_score(
        adata_path=adata_path,
        cluster_key=cluster_key,
        target_cell_type=target_cell_type,
        positive_markers=positive_markers if positive_markers else None,
        negative_markers=negative_markers if negative_markers else None,
    )


def _list_builtin_markers():
    """[只读] 列出内置细胞类型和类器官 marker panel."""
    from modules.cell_markers import list_builtin_cell_types
    from modules.annotation import ORGANOID_MARKER_SETS, ORGANOID_MARKER_SET_LABELS

    organoid_panels = {
        organoid_type: {
            'display': ORGANOID_MARKER_SET_LABELS.get(organoid_type, organoid_type),
            'cell_types': marker_dict,
        }
        for organoid_type, marker_dict in ORGANOID_MARKER_SETS.items()
    }
    return {
        "cell_types": list_builtin_cell_types(),
        "organoid_panels": organoid_panels,
    }


# ============================================================
# 阶段 5：参数搜索工具
# ============================================================

def _propose_parameter_sweep(args, project_id):
    """
    [需确认] 生成参数搜索候选列表。
    输入：
        {
            "goal_type": "target_cluster_refinement",
            "target_cell_type": "microglia",
            "current_state": {...},  // 可选，来自 inspect_analysis_state
            "max_candidates": 6
        }
    """
    goal_type = args.get('goal_type', 'target_cluster_refinement')
    target_cell_type = args.get('target_cell_type', '')
    max_candidates = min(int(args.get('max_candidates', 6)), 12)  # 硬上限 12

    # 获取当前状态
    current_state = args.get('current_state', {})
    if not current_state:
        current_state = _inspect_analysis_state(project_id)

    latest_adata = current_state.get('latest_adata') or _find_latest_adata(project_id)

    if not latest_adata:
        return {'error': '未找到可用 h5ad 文件'}

    candidates = _build_sweep_candidates(goal_type, target_cell_type, current_state, max_candidates)

    return {
        'base_checkpoint': latest_adata,
        'goal_type': goal_type,
        'target_cell_type': target_cell_type,
        'n_candidates': len(candidates),
        'candidates': candidates,
        'needs_confirmation': True,
        'message': f'提议 {len(candidates)} 个候选参数组合用于优化 {target_cell_type} 分群。运行需用户确认。',
    }


def _build_sweep_candidates(goal_type, target_cell_type, state, max_candidates):
    """构建参数搜索候选列表（确定性算法）."""
    candidates = []
    idx = 0

    if goal_type == 'target_cluster_refinement':
        resolutions = ['0.6', '0.8', '1.0', '1.2', '1.5']
        n_neighbors_options = ['10', '15', '30']

        # Resolution sweep
        for res in resolutions:
            if idx >= max_candidates:
                break
            candidates.append({
                'name': f"clustering_res_{res}",
                'modules': ['clustering'],
                'params': {
                    'clustering': {
                        'resolutions': res,
                        'n_neighbors': '15',
                    }
                },
                'rationale': f'调整 clustering resolution 为 {res}',
            })
            idx += 1

        # n_neighbors sweep (跳过已在 resolution 中包含的)
        for nn in n_neighbors_options:
            if idx >= max_candidates:
                break
            candidates.append({
                'name': f"clustering_res_1.0_nn_{nn}",
                'modules': ['clustering'],
                'params': {
                    'clustering': {
                        'resolutions': '1.0',
                        'n_neighbors': nn,
                    }
                },
                'rationale': f'调整 n_neighbors 为 {nn}',
            })
            idx += 1

    return candidates[:max_candidates]


def _run_parameter_sweep(args, project_id):
    """
    [需确认] 执行参数搜索。
    输入：
        {
            "goal_id": "...",
            "base_checkpoint": "...",
            "candidates": [...],
            "evaluator": "sc_cluster_signature"
        }
    """
    goal_id = args.get('goal_id', '')
    base_checkpoint = args.get('base_checkpoint', '')
    candidates = args.get('candidates', [])

    if not goal_id:
        return {'error': '缺少 goal_id'}

    if not candidates:
        return {'error': '候选列表为空'}

    # 验证 base_checkpoint 安全性
    err = _validate_project_path(base_checkpoint, project_id)
    if err:
        return {'error': f'base_checkpoint: {err}'}

    from modules.agent_orchestrator import run_parameter_sweep as do_sweep

    return do_sweep(
        goal_id=goal_id,
        candidates=candidates,
        base_checkpoint=base_checkpoint,
        project_id=project_id,
    )


# ============================================================
# 阶段 6：目标 Agent 编排工具
# ============================================================

def _start_goal_agent(args, project_id):
    """
    [需确认] 启动目标驱动的 agent 分析会话。
    输入：
        {
            "goal_type": "target_cluster_refinement",
            "target_cell_type": "microglia",
            "positive_markers": [...],
            "negative_markers": [...],
            "user_requirement": "找到最接近小胶质细胞的分群",
            "max_candidate_runs": 6
        }
    """
    goal_type = args.get('goal_type', 'target_cluster_refinement')
    target_cell_type = args.get('target_cell_type', '')
    positive_markers = args.get('positive_markers', None)
    negative_markers = args.get('negative_markers', None)
    user_requirement = args.get('user_requirement', '')
    max_candidate_runs = min(int(args.get('max_candidate_runs', 6)), 12)

    if not target_cell_type:
        return {'error': '缺少 target_cell_type 参数'}

    from modules.agent_orchestrator import start_goal_agent as do_start

    return do_start(
        project_id=project_id,
        goal_type=goal_type,
        target_cell_type=target_cell_type,
        positive_markers=positive_markers if positive_markers else None,
        negative_markers=negative_markers if negative_markers else None,
        user_requirement=user_requirement,
        max_candidate_runs=max_candidate_runs,
    )


def _continue_goal_agent(args, project_id):
    """
    [需确认] 继续目标分析会话，处理用户反馈。
    输入：
        {
            "session_id": "...",
            "instruction": "继续细分这个cluster" | "采用候选B" | "不满意，试试更高resolution" | "停止"
        }
    """
    session_id = args.get('session_id', '')
    instruction = args.get('instruction', '')

    if not session_id:
        return {'error': '缺少 session_id'}
    if not instruction:
        return {'error': '缺少 instruction'}

    from modules.agent_orchestrator import continue_goal_agent as do_continue

    return do_continue(session_id, instruction, project_id=project_id)
