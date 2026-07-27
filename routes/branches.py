# routes/branches.py
"""分析候选分支 API — 创建、运行、评分、采纳"""
import os
import json
import logging
from flask import Blueprint, request, jsonify
from config import Config
from routes.auth import require_ai_token

logger = logging.getLogger(__name__)

branches_bp = Blueprint('branches', __name__)


# ============================================================
# Agent Jobs API
# ============================================================

@branches_bp.route('/api/projects/<pid>/agent/jobs/<job_id>', methods=['GET'])
@require_ai_token
def get_agent_job(pid, job_id):
    """获取 agent job 状态。"""
    from models import AgentJob
    job = AgentJob.get_by_id(job_id)
    if not job:
        return jsonify({"error": "Job 不存在"}), 404
    if job.project_id != pid:
        return jsonify({"error": "Job 不属于该项目"}), 403
    return jsonify(job.to_dict())


@branches_bp.route('/api/projects/<pid>/agent/jobs', methods=['GET'])
@require_ai_token
def list_agent_jobs(pid):
    """列出项目的所有 agent jobs。"""
    from models import AgentJob
    jobs = AgentJob.get_by_project(pid)
    return jsonify({"jobs": [j.to_dict() for j in jobs]})


# ============================================================
# Current Context API
# ============================================================

@branches_bp.route('/api/projects/<pid>/current-context', methods=['GET'])
def get_current_context(pid):
    """获取项目当前分析基线."""
    from modules.ai_tools import resolve_current_adata_path
    from models import AnalysisBranch

    ctx = resolve_current_adata_path(pid, allow_upload=True)
    result = {
        "project_id": pid,
        "current_adata_path": ctx['path'],
        "source": ctx['source'],
        "can_continue_analysis": ctx['path'] is not None,
        "warnings": ctx.get('warnings', []),
    }

    if ctx.get('accepted_branch_id'):
        result['accepted_branch_id'] = ctx['accepted_branch_id']
        branch = AnalysisBranch.get_by_id(ctx['accepted_branch_id'])
        if branch:
            result['branch_name'] = branch.branch_name

    return jsonify(result)


# ============================================================
# Agent Sessions API
# ============================================================

@branches_bp.route('/api/projects/<pid>/agent/sessions', methods=['GET'])
@require_ai_token
def list_agent_sessions(pid):
    """列出项目的所有 agent 会话."""
    from models import AgentSession
    sessions = AgentSession.get_by_project(pid)
    return jsonify({"sessions": [s.to_dict() for s in sessions]})


@branches_bp.route('/api/projects/<pid>/agent/sessions', methods=['POST'])
@require_ai_token
def create_agent_session(pid):
    """创建 agent 会话."""
    from models import AgentSession
    data = request.get_json(silent=True) or {}
    session = AgentSession(
        project_id=pid,
        status='active',
        goal_summary=data.get('goal_summary', ''),
    )
    session.save()
    return jsonify(session.to_dict())


@branches_bp.route('/api/projects/<pid>/agent/sessions/<sid>', methods=['GET'])
@require_ai_token
def get_agent_session(pid, sid):
    """获取 agent 会话详情（含 goals 和 steps）。需校验归属项目."""
    from models import AgentSession, AgentGoal, AgentStep
    session = AgentSession.get_by_id(sid)
    if not session:
        return jsonify({"error": "会话不存在"}), 404
    if session.project_id != pid:
        return jsonify({"error": "会话不属于该项目"}), 403
    goals = AgentGoal.get_by_session(sid)
    steps = AgentStep.get_by_session(sid)
    return jsonify({
        "session": session.to_dict(),
        "goals": [g.to_dict() for g in goals],
        "steps": [s.to_dict() for s in steps],
    })


@branches_bp.route('/api/projects/<pid>/branches/<branch_id>/delete', methods=['POST'])
@require_ai_token
def delete_branch(pid, branch_id):
    """软删除候选分支（仅当未被采纳时）。使用 deleted=1 标记，不物理删除。"""
    from models import AnalysisBranch
    branch = AnalysisBranch.get_by_id(branch_id)
    if not branch:
        return jsonify({"error": "分支不存在"}), 404
    if branch.project_id != pid:
        return jsonify({"error": "分支不属于该项目"}), 403
    if branch.accepted:
        return jsonify({"error": "已采纳的分支不能删除"}), 400
    branch.soft_delete()
    return jsonify({"status": "deleted", "branch_id": branch_id})


def _validate_branch_path(path, project_id):
    """验证分支路径安全."""
    from config import Config as Cfg
    is_valid, err = Cfg._validate_path(path, project_id)
    if not is_valid:
        return err
    return None


@branches_bp.route('/api/projects/<pid>/branches', methods=['GET'])
@require_ai_token
def list_branches(pid):
    """列出项目的所有候选分支."""
    from models import AnalysisBranch
    branches = AnalysisBranch.get_by_project(pid)
    return jsonify({"branches": [b.to_dict() for b in branches]})


@branches_bp.route('/api/projects/<pid>/branches/compare', methods=['GET'])
@require_ai_token
def compare_branches(pid):
    """比较项目或目标下的候选分支评分和参数差异."""
    from models import AnalysisBranch, CandidateScore
    from modules.reporting.branch_report import build_branch_comparison

    goal_id = request.args.get('goal_id', '').strip()
    if goal_id:
        branches = [b for b in AnalysisBranch.get_by_goal(goal_id) if b.project_id == pid and not b.deleted]
    else:
        branches = AnalysisBranch.get_by_project(pid)

    scores_by_branch = {
        branch.id: CandidateScore.get_by_branch(branch.id)
        for branch in branches
    }
    return jsonify(build_branch_comparison(pid, branches, scores_by_branch))


@branches_bp.route('/api/projects/<pid>/branches', methods=['POST'])
@require_ai_token
def create_branch(pid):
    """创建候选分支."""
    data = request.get_json(silent=True) or {}
    parent_adata_path = data.get('parent_adata_path', '')
    branch_name = data.get('branch_name', '')
    purpose = data.get('purpose', '')
    session_id = data.get('session_id', None)
    goal_id = data.get('goal_id', None)

    if not parent_adata_path:
        return jsonify({"error": "缺少 parent_adata_path"}), 400

    err = _validate_branch_path(parent_adata_path, pid)
    if err:
        return jsonify({"error": err}), 400

    if not branch_name:
        import uuid
        branch_name = f"branch_{uuid.uuid4().hex[:8]}"

    # 确保 branches 目录存在
    branches_dir = Config.branches_dir(pid)
    os.makedirs(branches_dir, exist_ok=True)

    from models import AnalysisBranch
    branch = AnalysisBranch(
        project_id=pid,
        session_id=session_id,
        goal_id=goal_id,
        parent_adata_path=parent_adata_path,
        branch_name=branch_name,
        purpose=purpose,
        params_json='{}',
        status='pending',
    )
    branch.save()

    # 创建分支专属目录
    branch_dir = Config.branch_dir(pid, branch.id)
    os.makedirs(branch_dir, exist_ok=True)

    return jsonify({
        "id": branch.id,
        "branch_name": branch_name,
        "message": "候选分支已创建",
    })


@branches_bp.route('/api/projects/<pid>/branches/<branch_id>', methods=['GET'])
@require_ai_token
def get_branch(pid, branch_id):
    """获取分支详情（需校验归属项目）."""
    from models import AnalysisBranch, CandidateScore
    branch = AnalysisBranch.get_by_id(branch_id)
    if not branch:
        return jsonify({"error": "分支不存在"}), 404
    if branch.project_id != pid:
        return jsonify({"error": "分支不属于该项目"}), 403

    scores = CandidateScore.get_by_branch(branch_id)
    return jsonify({
        "branch": branch.to_dict(),
        "scores": [s.to_dict() for s in scores],
    })


@branches_bp.route('/api/projects/<pid>/branches/<branch_id>/run', methods=['POST'])
@require_ai_token
def run_branch(pid, branch_id):
    """在分支上运行分析模块."""
    from models import AnalysisBranch, AnalysisTask
    from modules import MODULE_REGISTRY, validate_pipeline_order
    from modules.ai_tools import _validate_analysis_params

    branch = AnalysisBranch.get_by_id(branch_id)
    if not branch:
        return jsonify({"error": "分支不存在"}), 404

    # 校验分支归属
    if branch.project_id != pid:
        return jsonify({"error": "分支不属于该项目"}), 403

    data = request.get_json(silent=True) or {}
    modules_list = data.get('modules', [])
    params_by_module = data.get('params', {})

    if not modules_list:
        return jsonify({"error": "缺少 modules 参数"}), 400

    # 校验模块
    for mod in modules_list:
        if mod not in MODULE_REGISTRY:
            return jsonify({"error": f"未知模块: {mod}"}), 400

    # 校验 pipeline order
    is_valid, errors = validate_pipeline_order(modules_list)
    if not is_valid:
        return jsonify({"error": f"模块顺序不满足依赖: {'; '.join(errors)}"}), 400

    # 校验并清洗参数（清洗后的参数用于执行）
    cleaned_params_by_module = {}
    for mod in modules_list:
        raw_params = params_by_module.get(mod, {})
        if raw_params:
            cleaned, err = _validate_analysis_params(mod, raw_params)
            if err:
                return jsonify({"error": f"模块 {mod} 参数校验失败: {err}"}), 400
            cleaned_params_by_module[mod] = cleaned
        else:
            cleaned_params_by_module[mod] = {}

    # 检查并标记 running
    if not branch.mark_running():
        return jsonify({"error": f"分支 {branch_id} 状态不是 pending，无法启动"}), 400

    # 在后台执行
    from worker import _executor, register_task_outputs
    import traceback

    branch_dir = Config.branch_dir(branch.project_id, branch.id)
    os.makedirs(branch_dir, exist_ok=True)

    def _run_branch_modules():
        current_task = None
        try:
            current_input = branch.parent_adata_path
            for mod_name in modules_list:
                params = cleaned_params_by_module.get(mod_name, {})
                cls = MODULE_REGISTRY[mod_name]

                task = AnalysisTask(
                    project_id=branch.project_id,
                    module_name=mod_name,
                    params_json=json.dumps(params, ensure_ascii=False),
                    branch_id=branch.id,
                )
                task.save()
                task.mark_running()
                current_task = task

                task_progress_log = []

                def progress_cb(pct, message, _task=task, _module=mod_name):
                    from datetime import datetime
                    now = datetime.now().strftime('%H:%M:%S')
                    task_progress_log.append({'time': now, 'pct': pct, 'msg': f'[{_module}] {message}'})
                    _task.update_progress(pct, message, json.dumps(task_progress_log, ensure_ascii=False))

                module = cls(project_dir=branch_dir, params=params, progress_callback=progress_cb)
                result = module.run(current_input)

                output_adata = result.get('output_adata')
                if not output_adata:
                    raise ValueError(f"模块 {mod_name} 未返回 output_adata")

                summary = result.get('summary', {})
                result_error = result.get('error') if isinstance(result, dict) else None
                summary_error = summary.get('error') if isinstance(summary, dict) else None
                if result_error or summary_error:
                    summary_json = json.dumps(summary, ensure_ascii=False, default=str)
                    error_message = str(result_error or summary_error)
                    task.mark_failed(error_message, summary_json)
                    raise ValueError(f"模块 {mod_name}：{error_message}")

                register_task_outputs(task, branch.project_id, branch_dir, result)
                task.mark_completed(output_adata, json.dumps(summary, ensure_ascii=False, default=str))
                current_input = output_adata
                current_task = None

            branch.mark_completed(current_input)
        except Exception as e:
            tb = traceback.format_exc()
            logger.error(f"[Branch] Branch {branch_id} failed:\n{tb}")
            if current_task is not None:
                try:
                    current_task.mark_failed(tb)
                except Exception as db_err:
                    logger.warning(f"[Branch] Failed to mark task {current_task.id} failed: {db_err}")
            try:
                branch.mark_failed(tb)
            except Exception:
                pass

    _executor.submit(_run_branch_modules)

    return jsonify({
        "branch_id": branch.id,
        "status": "running",
        "message": f"分支 {branch.branch_name} 开始运行: {modules_list}",
    })


@branches_bp.route('/api/projects/<pid>/branches/<branch_id>/accept', methods=['POST'])
@require_ai_token
def accept_branch(pid, branch_id):
    """采纳候选分支."""
    from models import AnalysisBranch, CandidateScore

    branch = AnalysisBranch.get_by_id(branch_id)
    if not branch:
        return jsonify({"error": "分支不存在"}), 404
    if branch.project_id != pid:
        return jsonify({"error": "分支不属于该项目"}), 403

    if branch.status != 'completed':
        return jsonify({"error": f"分支状态为 {branch.status}，只有 completed 分支可以采纳"}), 400

    data = request.get_json(silent=True) or {}
    confirmed = data.get('confirm', False)

    if not confirmed:
        return jsonify({
            "status": "needs_confirmation",
            "branch_id": branch.id,
            "branch_name": branch.branch_name,
            "message": "请显式确认采纳此分支（设置 confirm=true）",
        }), 200

    scores = CandidateScore.get_by_branch(branch.id)
    from modules.reporting.branch_report import build_branch_comparison, write_acceptance_report
    branch.accept()
    branch = AnalysisBranch.get_by_id(branch_id)
    project_branches = AnalysisBranch.get_by_project(pid)
    scores_by_branch = {
        b.id: CandidateScore.get_by_branch(b.id)
        for b in project_branches
    }
    comparison = build_branch_comparison(pid, project_branches, scores_by_branch)
    report_paths = write_acceptance_report(
        project_dir=Config.project_dir(pid),
        branch=branch,
        scores=scores,
        comparison=comparison,
    )

    return jsonify({
        "status": "accepted",
        "branch_id": branch.id,
        "branch_name": branch.branch_name,
        "output_adata_path": branch.output_adata_path,
        "accepted_report": report_paths,
        "message": f"已采纳分支 {branch.branch_name}，输出文件: {branch.output_adata_path}",
    })


@branches_bp.route('/api/projects/<pid>/branches/<branch_id>/score', methods=['POST'])
@require_ai_token
def score_branch(pid, branch_id):
    """对分支进行细胞类型签名评分."""
    from models import AnalysisBranch, CandidateScore
    from modules.evaluators.sc_cluster import score_cell_type_signature

    branch = AnalysisBranch.get_by_id(branch_id)
    if not branch:
        return jsonify({"error": "分支不存在"}), 404
    if branch.project_id != pid:
        return jsonify({"error": "分支不属于该项目"}), 403

    if branch.status != 'completed' or not branch.output_adata_path:
        return jsonify({"error": "分支未完成，无法评分"}), 400

    data = request.get_json(silent=True) or {}
    target_cell_type = data.get('target_cell_type', '')
    cluster_key = data.get('cluster_key', 'leiden')
    positive_markers = data.get('positive_markers', None)
    negative_markers = data.get('negative_markers', None)
    evaluator_name = data.get('evaluator_name', 'sc_cluster_signature')

    if not target_cell_type:
        return jsonify({"error": "缺少 target_cell_type"}), 400

    result = score_cell_type_signature(
        adata_path=branch.output_adata_path,
        cluster_key=cluster_key,
        target_cell_type=target_cell_type,
        positive_markers=positive_markers,
        negative_markers=negative_markers,
    )

    if result.get('error'):
        return jsonify({"error": result['error']}), 400

    best = result.get('best_cluster')
    total_score = result['cluster_scores'][0]['total_score'] if result.get('cluster_scores') else None

    # 保存评分
    score = CandidateScore(
        branch_id=branch_id,
        project_id=branch.project_id,
        evaluator_name=evaluator_name,
        target_label=f"{target_cell_type}:{best}" if best else target_cell_type,
        score_json=json.dumps(result, ensure_ascii=False),
        total_score=total_score,
        recommendation=f"最佳 cluster: {best}, 置信度: {result.get('confidence', 'N/A')}" if best else "未找到匹配 cluster",
    )
    score.save()

    return jsonify({
        "branch_id": branch_id,
        "score": score.to_dict(),
        "detail": result,
    })
