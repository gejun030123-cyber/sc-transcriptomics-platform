# modules/ai_tools.py
"""AI 工具函数 — 平台分析功能的接口"""
import os
import json
from config import Config


def execute_tool(name, args, project_id):
    """执行指定的工具函数"""
    if name == "run_analysis":
        return _run_analysis(args, project_id)
    elif name == "get_project_status":
        return _get_project_status(project_id)
    elif name == "get_task_results":
        return _get_task_results(args.get("task_id", ""))
    elif name == "list_modules":
        return _list_modules(args.get("pipeline_type", "all"))
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
        # 使用统一 context resolver（优先 accepted branch > 主线任务 > intermediate > uploads）
        ctx = resolve_current_adata_path(project_id, allow_upload=True)
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
    if os.path.isdir(uploads_dir):
        uploaded = [f for f in os.listdir(uploads_dir)
                    if f.endswith(('.h5ad', '.csv', '.txt', '.xlsx'))]

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


def _validate_analysis_params(module_name, params):
    """校验 AI 传入的分析参数，移除未知键，返回 (cleaned_params, error_msg)。"""
    from modules.schemas import PARAM_SCHEMAS

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

    return cleaned, None


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
    [只读] 检查指定 AnnData 文件的详细结构。
    输入：{"adata_path": "..."}
    """
    adata_path = args.get("adata_path", "")
    if not adata_path:
        ctx = resolve_current_adata_path(project_id, allow_upload=True)
        if ctx['path']:
            adata_path = ctx['path']
        else:
            return {"error": "未找到当前 h5ad 文件，请先上传数据或指定路径"}
    err = _validate_project_path(adata_path, project_id)
    if err:
        return {"error": err}

    try:
        import scanpy as sc
        adata = sc.read(adata_path, backed='r')

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
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "obs_columns": list(adata.obs.columns)[:50],
            "var_columns": list(adata.var.columns)[:20],
            "layers": layers,
            "obsm": obsm_keys,
            "cluster_keys": cluster_keys,
            "embedding_keys": obsm_keys,
            "sample_preview": preview,
        }
    except Exception as e:
        return {"error": f"无法读取 AnnData: {e}"}


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
    """[只读] 列出所有内置细胞类型 marker 定义."""
    from modules.cell_markers import list_builtin_cell_types
    return {"cell_types": list_builtin_cell_types()}


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
