# modules/agent_jobs.py
"""Agent Job 执行器 — 异步后台任务（sweep、branch run 等）"""
import os
import json
import traceback
import logging
from models import gen_id, AgentJob, AnalysisBranch, AnalysisTask, CandidateScore, AgentStep

logger = logging.getLogger(__name__)


def submit_sweep_job(goal_id, candidates, base_checkpoint, project_id):
    """
    提交参数搜索后台任务。立即返回 job_id，不阻塞调用方。

    返回：
        dict: {status: 'submitted', job_id, goal_id, n_candidates, poll_url}
    """
    from modules.agent_orchestrator import _validate_sweep_inputs

    # 执行期校验（与同步版相同）
    error = _validate_sweep_inputs(goal_id, candidates, base_checkpoint, project_id)
    if error:
        return error

    from modules.agent_orchestrator import _validate_and_clean_candidates
    validated = _validate_and_clean_candidates(candidates)
    if isinstance(validated, dict) and 'error' in validated:
        return validated

    # 读取 goal 获取 session_id
    from models import AgentGoal
    goal = AgentGoal.get_by_id(goal_id)
    session_id = goal.session_id if goal else None

    # 创建 job
    job = AgentJob(
        project_id=project_id,
        session_id=session_id,
        goal_id=goal_id,
        job_type='parameter_sweep',
        status='pending',
        params_json=json.dumps({
            'goal_id': goal_id,
            'n_candidates': len(validated),
            'base_checkpoint': base_checkpoint,
        }, ensure_ascii=False),
    )
    job.save()

    # 提交后台执行
    from worker import _executor
    _executor.submit(_run_sweep_background, job, goal_id, validated, base_checkpoint, project_id)

    return {
        'status': 'submitted',
        'job_id': job.id,
        'goal_id': goal_id,
        'n_candidates': len(validated),
        'poll_url': f'/api/projects/{project_id}/agent/jobs/{job.id}',
        'message': f'参数搜索已提交，共 {len(validated)} 个候选',
    }


def _run_sweep_background(job, goal_id, candidates, base_checkpoint, project_id):
    """后台执行参数搜索（在 worker 线程池中运行）。顶层兜底确保任何异常都让 job 进入 failed。"""
    from modules import MODULE_REGISTRY

    if not job.mark_running():
        logger.error(f"[AgentJob] Job {job.id} not in pending state")
        return

    results = []
    total = len(candidates)
    errors = []

    try:
        _run_sweep_loop(job, goal_id, candidates, base_checkpoint, project_id, results, errors, total)
    except Exception as e:
        tb = traceback.format_exc()
        logger.error(f"[AgentJob] Job {job.id} crashed:\n{tb}")
        errors.append({'candidate': '__job__', 'error': str(e)})
        result_json = json.dumps({
            'n_total': total, 'n_completed': len([r for r in results if r.get('status') == 'completed']),
            'n_failed': len(errors), 'all_results': results, 'errors': errors,
        }, ensure_ascii=False)
        try:
            # 先持久化 partial result，便于 UI 展示部分完成/失败候选
            job.update_progress(job.progress, f'崩溃: {str(e)[:80]}', result_json)
            job.mark_failed(tb)
        except Exception:
            pass


def _run_sweep_loop(job, goal_id, candidates, base_checkpoint, project_id, results, errors, total):
    """逐候选执行 sweep（在顶层 try/except 内）。"""
    from modules import MODULE_REGISTRY
    from worker import register_task_outputs

    for i, candidate in enumerate(candidates):
        branch_name = candidate['name']
        modules_list = candidate['modules']
        params = candidate['params']

        pct = int((i / total) * 100)
        job.update_progress(
            pct,
            f'运行候选 {i+1}/{total}: {branch_name}',
            json.dumps({'results': results, 'errors': errors}, ensure_ascii=False),
        )

        # 创建 branch
        branch = AnalysisBranch(
            project_id=project_id,
            session_id=job.session_id,
            goal_id=goal_id,
            parent_adata_path=base_checkpoint,
            branch_name=branch_name,
            purpose=candidate.get('rationale', ''),
            params_json=json.dumps(params, ensure_ascii=False),
            status='pending',
        )
        branch.save()

        # 写入 step
        step = AgentStep(
            session_id=job.session_id or '',
            goal_id=goal_id,
            step_index=i,
            step_type='candidate_started',
            thought_summary=f'Job {job.id}: 启动候选 {branch_name}',
            action_name='run_branch_candidate',
            action_args_json=json.dumps({'branch_id': branch.id, 'modules': modules_list}),
            status='pending',
        )
        step.save()

        # 运行模块
        from config import Config
        branch_dir = Config.branch_dir(project_id, branch.id)
        os.makedirs(branch_dir, exist_ok=True)

        if not branch.mark_running():
            errors.append({'candidate': branch_name, 'error': '无法启动 branch（非 pending 状态）'})
            step.status = 'failed'
            step.result_json = json.dumps({'error': 'branch not pending'})
            step.save()
            results.append({'branch_id': branch.id, 'branch_name': branch_name, 'status': 'failed', 'error': '无法启动 branch'})
            continue

        try:
            current_input = base_checkpoint
            current_task = None
            for mod_name in modules_list:
                cls = MODULE_REGISTRY.get(mod_name)
                if not cls:
                    raise ValueError(f'未知模块: {mod_name}')

                task = AnalysisTask(
                    project_id=project_id,
                    module_name=mod_name,
                    params_json=json.dumps(params.get(mod_name, {}), ensure_ascii=False),
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

                module = cls(project_dir=branch_dir, params=params.get(mod_name, {}), progress_callback=progress_cb)
                result = module.run(current_input)

                output_adata = result.get('output_adata')
                if not output_adata:
                    raise ValueError(f'模块 {mod_name} 未返回 output_adata')

                register_task_outputs(task, project_id, branch_dir, result)
                task.mark_completed(output_adata, json.dumps(result.get('summary', {}), ensure_ascii=False))
                current_input = output_adata
                current_task = None

            branch.mark_completed(current_input)

            # 自动评分
            from modules.evaluators.sc_cluster import score_cell_type_signature
            from models import AgentGoal
            goal = AgentGoal.get_by_id(goal_id)
            target_cell_type = 'microglia'
            if goal:
                try:
                    goal_data = json.loads(goal.goal_json or '{}')
                    target_cell_type = goal_data.get('target_cell_type', goal.target_entity or 'microglia')
                except Exception:
                    pass

            score_result = score_cell_type_signature(
                adata_path=current_input,
                cluster_key='leiden',
                target_cell_type=target_cell_type,
            )

            top_score = None
            if not score_result.get('error') and score_result.get('cluster_scores'):
                top = score_result['cluster_scores'][0]
                top_score = top['total_score']
                cs = CandidateScore(
                    branch_id=branch.id,
                    project_id=project_id,
                    evaluator_name='sc_cluster_signature',
                    target_label=f"{target_cell_type}:{top['cluster']}",
                    score_json=json.dumps(score_result, ensure_ascii=False),
                    total_score=top_score,
                    recommendation=f"最佳 cluster: {top['cluster']}, 评分: {top_score:.3f}",
                )
                cs.save()

            step.status = 'completed'
            step.result_json = json.dumps({'branch_id': branch.id, 'status': 'completed', 'score': top_score})
            step.save()

            results.append({
                'branch_id': branch.id, 'branch_name': branch_name,
                'status': 'completed', 'score': top_score,
                'best_cluster': score_result.get('best_cluster'),
                'confidence': score_result.get('confidence'),
            })

        except Exception as e:
            tb = traceback.format_exc()
            logger.error(f"[AgentJob] Candidate {branch_name} failed:\n{tb}")
            if 'current_task' in locals() and current_task is not None:
                try:
                    current_task.mark_failed(tb)
                except Exception as db_err:
                    logger.warning(f"[AgentJob] Failed to mark task {current_task.id} failed: {db_err}")
            try:
                branch.mark_failed(tb)
            except Exception:
                pass
            step.status = 'failed'
            step.result_json = json.dumps({'error': str(e)})
            step.save()

            errors.append({'candidate': branch_name, 'error': str(e)})
            results.append({'branch_id': branch.id, 'branch_name': branch_name, 'status': 'failed', 'error': str(e)})

    # 全部完成
    completed = [r for r in results if r.get('status') == 'completed' and r.get('score') is not None]
    best = None
    if completed:
        completed.sort(key=lambda x: x.get('score', 0), reverse=True)
        best = completed[0]

    final_result = json.dumps({
        'n_total': len(results),
        'n_completed': len(completed),
        'n_failed': len(errors),
        'best_candidate': best,
        'all_results': results,
        'errors': errors,
    }, ensure_ascii=False)

    if errors and not completed:
        job.mark_failed(f'所有候选均失败: {errors[0].get("error", "unknown")}')
    else:
        job.mark_completed(final_result)
