# modules/agent_orchestrator.py
"""Agent 编排器 — 目标驱动的分析流程控制

核心流程:
1. parse_goal → 解析用户目标
2. inspect_analysis_state → 了解当前状态
3. inspect current clusters → 检查分群
4. score current result → 评分当前结果
5. if good enough: 展示证据, 用户确认
   else: propose parameter sweep
6. 用户确认 sweep
7. run candidate branches → 逐步执行
8. evaluate candidates → 自动评分
9. recommend best branch → 推荐最佳分支
10. 等待用户 accept/refine/stop
"""

import os
import json
import logging
from datetime import datetime, timezone

from models import (
    gen_id, AgentSession, AgentGoal, AgentStep, AnalysisBranch, CandidateScore
)
from config import Config

logger = logging.getLogger(__name__)

# 参数搜索空间（阶段 5 首版）
SWEEP_PARAMS_CLUSTERING = {
    'resolutions': ['0.6', '0.8', '1.0', '1.2', '1.5'],
    'n_neighbors': ['10', '15', '30'],
}

SWEEP_PARAMS_DIMRED = {
    'n_comps': ['30', '50', '80'],
}

SWEEP_PARAMS_HVG = {
    'n_top_genes': ['2000', '3000', '5000'],
}


def _now():
    return datetime.now(timezone.utc).isoformat()


def start_goal_agent(project_id, goal_type, target_cell_type,
                     positive_markers=None, negative_markers=None,
                     user_requirement='', max_candidate_runs=6):
    """
    启动一个目标分析会话。

    参数：
        project_id: str
        goal_type: str, 如 'target_cluster_refinement'
        target_cell_type: str, 如 'microglia'
        positive_markers: list[str]|None
        negative_markers: list[str]|None
        user_requirement: str, 用户自然语言需求
        max_candidate_runs: int, 最大候选数

    返回：
        dict: {session_id, goal_id, status, summary, proposed_actions, needs_confirmation}
    """
    # 1. 创建 session
    session = AgentSession(
        project_id=project_id,
        status='active',
        goal_summary=f"{goal_type}: {target_cell_type} — {user_requirement[:80]}",
    )
    session.save()

    # 2. 创建 goal
    goal_json = json.dumps({
        'goal_type': goal_type,
        'target_cell_type': target_cell_type,
        'positive_markers': positive_markers or [],
        'negative_markers': negative_markers or [],
        'user_requirement': user_requirement,
        'max_candidate_runs': max_candidate_runs,
    }, ensure_ascii=False)

    goal = AgentGoal(
        session_id=session.id,
        project_id=project_id,
        goal_type=goal_type,
        target_entity=target_cell_type,
        goal_json=goal_json,
        status='draft',
    )
    goal.save()

    # 3. 记录 step 0: 目标创建
    step = AgentStep(
        session_id=session.id,
        goal_id=goal.id,
        step_index=0,
        step_type='goal_created',
        thought_summary=f"用户目标：{goal_type}，靶细胞类型：{target_cell_type}",
        action_name='start_goal_agent',
        action_args_json=json.dumps({
            'goal_type': goal_type,
            'target_cell_type': target_cell_type,
            'max_candidate_runs': max_candidate_runs,
        }, ensure_ascii=False),
        result_json=json.dumps({'session_id': session.id, 'goal_id': goal.id}),
        status='completed',
    )
    step.save()

    # 4. 执行初始检查
    from modules.ai_tools import _find_latest_adata, _inspect_analysis_state

    state = _inspect_analysis_state(project_id)
    latest_adata = state.get('latest_adata')

    if not latest_adata:
        return {
            'session_id': session.id,
            'goal_id': goal.id,
            'status': 'blocked',
            'summary': '项目中没有可用的 h5ad 文件。请先上传数据并运行基础分析流程（至少到 clustering）。',
            'proposed_actions': [],
            'needs_confirmation': False,
        }

    # 5. 评分当前结果
    from modules.evaluators.sc_cluster import score_cell_type_signature

    cluster_key = 'leiden'
    if state.get('available_cluster_keys'):
        cluster_key = state['available_cluster_keys'][0]

    score_result = score_cell_type_signature(
        adata_path=latest_adata,
        cluster_key=cluster_key,
        target_cell_type=target_cell_type,
        positive_markers=positive_markers,
        negative_markers=negative_markers,
    )

    # 6. 记录 step 1: 初始评分
    step1 = AgentStep(
        session_id=session.id,
        goal_id=goal.id,
        step_index=1,
        step_type='initial_assessment',
        thought_summary=f"检查当前 {cluster_key} 分群中 {target_cell_type} 签名评分",
        action_name='score_cell_type_signature',
        action_args_json=json.dumps({
            'adata_path': latest_adata,
            'cluster_key': cluster_key,
            'target_cell_type': target_cell_type,
        }, ensure_ascii=False),
        result_json=json.dumps({
            'best_cluster': score_result.get('best_cluster'),
            'confidence': score_result.get('confidence'),
            'top_score': score_result['cluster_scores'][0]['total_score'] if score_result.get('cluster_scores') else None,
        }, ensure_ascii=False),
        status='completed',
    )
    step1.save()

    # 7. 判断是否需要 sweep
    confidence = score_result.get('confidence', 'none')
    top_score = score_result['cluster_scores'][0]['total_score'] if score_result.get('cluster_scores') else 0

    if confidence in ('high',) and top_score >= 0.75:
        # 当前结果已满足目标
        session.status = 'completed'
        session.save()
        goal.status = 'achieved'
        goal.save()

        return {
            'session_id': session.id,
            'goal_id': goal.id,
            'status': 'goal_achieved',
            'summary': (
                f"当前 {cluster_key} 分群中 cluster {score_result['best_cluster']} "
                f"已满足 {target_cell_type} 细胞类型要求（评分: {top_score:.2f}, 置信度: {confidence}）。"
            ),
            'proposed_actions': [],
            'needs_confirmation': False,
            'current_score': score_result,
        }

    # 需要参数搜索
    proposals = _generate_sweep_candidates(
        goal, state, score_result, max_candidate_runs
    )

    session.goal_summary = f"{goal_type}: {target_cell_type} — 需要参数搜索优化"
    session.save()

    step2 = AgentStep(
        session_id=session.id,
        goal_id=goal.id,
        step_index=2,
        step_type='sweep_proposed',
        thought_summary=f"当前评分 ({top_score:.2f}) 未达目标，提议 {len(proposals)} 个参数候选",
        action_name='propose_parameter_sweep',
        action_args_json=json.dumps({
            'current_confidence': confidence,
            'current_top_score': top_score,
            'max_candidates': max_candidate_runs,
        }, ensure_ascii=False),
        result_json=json.dumps({
            'n_candidates': len(proposals),
            'candidates': proposals,
        }, ensure_ascii=False),
        status='completed',
    )
    step2.save()

    return {
        'session_id': session.id,
        'goal_id': goal.id,
        'status': 'needs_confirmation',
        'summary': (
            f"当前 {target_cell_type} 评分 {top_score:.2f}（置信度: {confidence}）。"
            f"建议进行参数搜索（{len(proposals)} 个候选）以优化分群结果。"
        ),
        'current_score': score_result,
        'proposed_actions': [
            {
                'action': 'run_parameter_sweep',
                'description': f'对 {target_cell_type} 目标运行 {len(proposals)} 个候选参数组合',
                'candidates': proposals,
                'needs_confirmation': True,
            }
        ],
        'needs_confirmation': True,
    }


def _generate_sweep_candidates(goal, state, current_score, max_candidates=6):
    """
    生成参数搜索候选列表。
    策略：
    - 如果已有 X_pca 和 X_umap，优先只 sweep clustering
    - 如果分群结构明显不稳定，加入 dimred 参数
    - 默认最多 max_candidates 个候选
    """
    candidates = []
    goal_json = json.loads(goal.goal_json) if goal.goal_json else {}
    target = goal_json.get('target_cell_type', '')

    available_embeddings = state.get('available_embeddings', [])
    has_embeddings = bool(available_embeddings)

    # 优先只 sweep clustering（最快）
    resolutions = ['0.6', '0.8', '1.0', '1.2', '1.5']
    n_neighbors_opts = ['10', '15', '30']

    idx = 0
    for res in resolutions:
        if idx >= max_candidates:
            break
        # resolution sweep
        candidates.append({
            'name': f"clustering_res_{res}",
            'modules': ['clustering'],
            'params': {
                'clustering': {
                    'resolutions': res,
                    'n_neighbors': '15',
                }
            },
            'rationale': f'尝试 resolution={res} 来调整分群粒度',
        })
        idx += 1

    # 如果有 budget，添加 neighbor 变化
    for nn in n_neighbors_opts:
        if idx >= max_candidates:
            break
        if nn == '15':
            continue  # 已在上面作为默认值
        candidates.append({
            'name': f"clustering_nn_{nn}",
            'modules': ['clustering'],
            'params': {
                'clustering': {
                    'resolutions': '1.0',
                    'n_neighbors': nn,
                }
            },
            'rationale': f'尝试 n_neighbors={nn} 来调整局部/全局结构',
        })
        idx += 1

    return candidates[:max_candidates]


def _validate_sweep_inputs(goal_id, candidates, base_checkpoint, project_id):
    """校验 sweep 输入合法性。纯逻辑校验在前，I/O 校验在后。返回 error dict 或 None。"""
    from modules.ai_tools import _validate_project_path

    goal = AgentGoal.get_by_id(goal_id)
    if not goal:
        return {'error': 'Goal 不存在'}
    if goal.project_id != project_id:
        return {'error': f'Goal {goal_id} 不属于项目 {project_id}'}

    max_allowed = 12
    if len(candidates) > max_allowed:
        return {'error': f'候选数量 ({len(candidates)}) 超过硬上限 ({max_allowed})'}

    # 路径校验放在纯逻辑校验之后（I/O 操作最后）
    err = _validate_project_path(base_checkpoint, project_id)
    if err:
        return {'error': f'base_checkpoint: {err}'}

    return None


def _validate_and_clean_candidates(candidates):
    """校验并清洗候选参数。返回清洁后的 candidates 列表或 error dict。"""
    from modules.ai_tools import _validate_analysis_params
    from modules import MODULE_REGISTRY, validate_pipeline_order

    ALLOWED_MODULES_FIRST_ROUND = {'clustering', 'dimred', 'hvg', 'batch_correct'}
    validated = []

    for i, candidate in enumerate(candidates):
        modules_list = candidate.get('modules', ['clustering'])
        params = candidate.get('params', {})

        for mod in modules_list:
            if mod not in MODULE_REGISTRY:
                return {'error': f'候选 {i}: 未知模块 {mod}'}
            if mod not in ALLOWED_MODULES_FIRST_ROUND:
                return {'error': f'候选 {i}: 模块 {mod} 不在第一轮允许范围（{ALLOWED_MODULES_FIRST_ROUND}）'}

        is_valid, errors = validate_pipeline_order(modules_list)
        if not is_valid:
            return {'error': f'候选 {i}: pipeline order 不满足依赖: {"; ".join(errors)}'}

        cleaned_params = {}
        for mod in modules_list:
            raw = params.get(mod, {})
            if raw:
                cleaned, err_clean = _validate_analysis_params(mod, raw)
                if err_clean:
                    return {'error': f'候选 {i} 模块 {mod}: 参数校验失败: {err_clean}'}
                cleaned_params[mod] = cleaned
            else:
                cleaned_params[mod] = {}

        validated.append({
            'name': candidate['name'],
            'modules': modules_list,
            'params': cleaned_params,
            'rationale': candidate.get('rationale', ''),
        })

    return validated


def run_parameter_sweep(goal_id, candidates, base_checkpoint, project_id):
    """
    执行参数搜索（异步）：创建 job 后立即返回，不阻塞。
    使用 modules/agent_jobs.py 中的后台执行器。
    """
    from modules.agent_jobs import submit_sweep_job
    return submit_sweep_job(goal_id, candidates, base_checkpoint, project_id)




def continue_goal_agent(session_id, user_instruction, project_id=None):
    """
    继续目标分析会话，处理用户新指令。
    支持指令：继续细分 / 不满意 / 采用分支 / 停止 / 增加约束

    参数：
        project_id: 用于校验 session 归属，防止跨项目访问

    返回：
        dict: {status, summary, proposed_actions}
    """
    session = AgentSession.get_by_id(session_id)
    if not session:
        return {'error': 'Session 不存在', 'status': 'error'}

    # 校验 session 归属于当前 project
    if project_id and session.project_id != project_id:
        return {'error': f'Session {session_id} 不属于项目 {project_id}', 'status': 'error'}

    goals = AgentGoal.get_by_session(session_id)
    if not goals:
        return {'error': 'Session 中没有 goal', 'status': 'error'}

    goal = goals[0]  # 取最新 goal
    step_idx = AgentStep.next_index(session_id)
    instruction_lower = user_instruction.lower()

    # 解析用户指令
    if '采用' in user_instruction or 'accept' in instruction_lower:
        # 提取分支标识
        branches = AnalysisBranch.get_by_goal(goal.id)
        completed = [b for b in branches if b.status == 'completed']

        if not completed:
            return {
                'session_id': session_id,
                'goal_id': goal.id,
                'status': 'blocked',
                'summary': '没有已完成的候选分支可以采纳。',
                'proposed_actions': [],
                'needs_confirmation': False,
            }

        # 创建 step
        step = AgentStep(
            session_id=session_id,
            goal_id=goal.id,
            step_index=step_idx,
            step_type='user_instruction',
            thought_summary=f"用户指令: {user_instruction}",
            action_name='continue_goal_agent',
            action_args_json=json.dumps({'instruction': user_instruction}, ensure_ascii=False),
            result_json=json.dumps({'available_branches': [b.id for b in completed]}),
            status='completed',
        )
        step.save()

        return {
            'session_id': session_id,
            'goal_id': goal.id,
            'status': 'needs_confirmation',
            'summary': f'找到 {len(completed)} 个已完成的候选分支。请选择要采纳的分支。',
            'available_branches': [{
                'branch_id': b.id,
                'branch_name': b.branch_name,
                'purpose': b.purpose,
                'status': b.status,
            } for b in completed],
            'proposed_actions': [{
                'action': 'accept_branch',
                'description': '采纳选定的候选分支',
                'needs_confirmation': True,
            }],
            'needs_confirmation': True,
        }

    elif '停止' in user_instruction or 'stop' in instruction_lower:
        session.status = 'completed'
        session.updated_at = _now()
        session.save()

        step = AgentStep(
            session_id=session_id,
            goal_id=goal.id,
            step_index=step_idx,
            step_type='user_instruction',
            thought_summary=f"用户停止目标分析",
            action_name='continue_goal_agent',
            action_args_json=json.dumps({'instruction': user_instruction}),
            result_json=json.dumps({'status': 'stopped'}),
            status='completed',
        )
        step.save()

        return {
            'session_id': session_id,
            'goal_id': goal.id,
            'status': 'stopped',
            'summary': '目标分析已停止。',
            'proposed_actions': [],
            'needs_confirmation': False,
        }

    elif '继续' in user_instruction or '细分' in user_instruction or '不满意' in user_instruction or '更高' in user_instruction or '必须' in user_instruction:
        # 更新约束
        try:
            constraints = json.loads(goal.constraints_json) if goal.constraints_json else {}
        except Exception:
            constraints = {}

        if '必须' in user_instruction:
            constraints['hard_requirement'] = user_instruction
        if '不用' in user_instruction or '不要' in user_instruction:
            constraints['exclude'] = constraints.get('exclude', []) + [user_instruction]

        goal.constraints_json = json.dumps(constraints, ensure_ascii=False)
        goal.updated_at = _now()
        goal.save()

        step = AgentStep(
            session_id=session_id,
            goal_id=goal.id,
            step_index=step_idx,
            step_type='user_instruction',
            thought_summary=f"用户新增约束: {user_instruction}",
            action_name='continue_goal_agent',
            action_args_json=json.dumps({'instruction': user_instruction}, ensure_ascii=False),
            result_json=json.dumps({'updated_constraints': constraints}),
            status='completed',
        )
        step.save()

        return {
            'session_id': session_id,
            'goal_id': goal.id,
            'status': 'constraints_updated',
            'summary': f'已更新约束条件。可以运行新的参数搜索以找到更匹配的分群。',
            'updated_constraints': constraints,
            'proposed_actions': [{
                'action': 'run_parameter_sweep',
                'description': '基于新约束重新运行参数搜索',
                'needs_confirmation': True,
            }],
            'needs_confirmation': True,
        }

    else:
        # 未识别的指令
        step = AgentStep(
            session_id=session_id,
            goal_id=goal.id,
            step_index=step_idx,
            step_type='user_instruction',
            thought_summary=f"未识别指令: {user_instruction}",
            action_name='continue_goal_agent',
            action_args_json=json.dumps({'instruction': user_instruction}),
            result_json=json.dumps({'error': 'unrecognized_instruction'}),
            status='completed',
        )
        step.save()

        return {
            'session_id': session_id,
            'goal_id': goal.id,
            'status': 'unknown_instruction',
            'summary': f'无法解析指令: "{user_instruction}"。支持的指令: "采用候选X", "继续细分", "不满意，试试更高resolution", "停止"。',
            'proposed_actions': [],
            'needs_confirmation': False,
        }
