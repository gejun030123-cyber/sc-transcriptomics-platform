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
    if not input_path:
        # 自动查找最新的中间文件
        intermediate_dir = os.path.join(project_dir, 'intermediate')
        if os.path.isdir(intermediate_dir):
            h5ad_files = sorted(
                [f for f in os.listdir(intermediate_dir) if f.endswith('.h5ad')],
                key=lambda f: os.path.getmtime(os.path.join(intermediate_dir, f)),
                reverse=True
            )
            if h5ad_files:
                input_path = os.path.join(intermediate_dir, h5ad_files[0])
        if not input_path:
            # 尝试上传目录
            uploads_dir = os.path.join(project_dir, 'uploads')
            if os.path.isdir(uploads_dir):
                upload_files = [f for f in os.listdir(uploads_dir)
                                if f.endswith(('.h5ad', '.csv', '.txt', '.xlsx'))]
                if upload_files:
                    input_path = os.path.join(uploads_dir, upload_files[0])

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
        "message": f"分析任务已提交：{module_name}，任务 ID: {task.id}"
    }


def _get_project_status(project_id):
    """获取项目状态"""
    from models import AnalysisTask

    project_dir = Config.project_dir(project_id)
    if not os.path.isdir(project_dir):
        return {"error": "项目不存在"}

    tasks = AnalysisTask.get_by_project(project_id)

    # 已上传文件
    uploads_dir = os.path.join(project_dir, 'uploads')
    uploaded = []
    if os.path.isdir(uploads_dir):
        uploaded = [f for f in os.listdir(uploads_dir)
                    if f.endswith(('.h5ad', '.csv', '.txt', '.xlsx'))]

    # 已完成任务
    completed_tasks = []
    for t in tasks:
        if t.status == 'completed':
            completed_tasks.append({
                "id": t.id,
                "module": t.module_name,
                "finished_at": str(t.finished_at) if t.finished_at else None,
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
