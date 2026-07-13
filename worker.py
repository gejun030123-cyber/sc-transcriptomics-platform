import concurrent.futures
import traceback
import json
import logging
from datetime import datetime
from database import get_conn
from models import gen_id, AnalysisTask, ResultFile, PipelineRun

logger = logging.getLogger(__name__)

_executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
_active_futures = {}

def active_count():
    return sum(1 for f in _active_futures.values() if not f.done())


def register_task_outputs(task, project_id, project_dir, result, pipeline_run_id=None):
    """Persist result files, converting legacy Plotly payloads to static figures."""
    created_result_files = []
    for rf in (result.get('result_files') or []):
        try:
            if rf.get('file_type') == 'plotly_json' and rf.get('file_path'):
                from modules.reporting.static_rendering import render_plotly_json
                static_files = render_plotly_json(
                    rf['file_path'], label=rf.get('label') or rf.get('category') or '分析图'
                )
                for static_file in static_files:
                    static_file['category'] = rf.get('category', 'plot')
                    created_result_files.append(ResultFile.create(
                        task_id=task.id, project_id=project_id,
                        file_type=static_file['file_type'],
                        category=static_file['category'],
                        label=static_file['label'],
                        file_path=static_file['file_path']
                    ))
                # The JSON is an implementation detail after conversion and must
                # not remain available as an interactive result artifact.
                try:
                    import os
                    os.remove(rf['file_path'])
                except OSError:
                    pass
                continue
            created_result_files.append(ResultFile.create(
                task_id=task.id, project_id=project_id,
                file_type=rf.get('file_type', ''),
                category=rf.get('category', ''),
                label=rf.get('label', ''),
                file_path=rf.get('file_path', '')
            ))
        except Exception as e:
            logger.warning(f"[Worker] Skipping result_files insert: {rf.get('file_path', '')} ({e})")

    try:
        from modules.reporting.result_manifest import write_task_manifest
        manifest_rf = write_task_manifest(
            project_dir=project_dir,
            task=task,
            result=result,
            result_files=created_result_files or result.get('result_files') or [],
            pipeline_run_id=pipeline_run_id,
        )
        created_result_files.append(ResultFile.create(
            task_id=task.id, project_id=project_id,
            file_type=manifest_rf.get('file_type', ''),
            category=manifest_rf.get('category', ''),
            label=manifest_rf.get('label', ''),
            file_path=manifest_rf.get('file_path', '')
        ))
    except Exception as e:
        logger.warning(f"[Worker] Failed to write task manifest for {task.id}: {e}")

    return created_result_files


def submit_task(task_id, project_id, module_name, params, project_dir, input_path):
    if task_id in _active_futures:
        logger.warning(f"[Worker] Task {task_id} already submitted, skipping")
        return False
    future = _executor.submit(
        _run_task, task_id, project_id, module_name,
        params, project_dir, input_path
    )
    _active_futures[task_id] = future
    return True

def _run_task(task_id, project_id, module_name, params, project_dir, input_path):
    task = None
    try:
        task = AnalysisTask.get_by_id(task_id)
        if not task:
            logger.warning(f"[Worker] Task {task_id} not found")
            return

        proj_conn = get_conn()
        try:
            proj_row = proj_conn.execute("SELECT id FROM projects WHERE id=?", (project_id,)).fetchone()
            logger.info(f"[Worker] task_id={task_id} exists=True, project_id={project_id} exists={proj_row is not None}")
        finally:
            proj_conn.close()

        if not task.mark_running():
            logger.warning(f"[Worker] Task {task_id} not in pending state, skipping")
            return

        progress_log = []

        def progress_cb(pct, message):
            now = datetime.now().strftime('%H:%M:%S')
            progress_log.append({'time': now, 'pct': pct, 'msg': message})
            task.update_progress(pct, message, json.dumps(progress_log, ensure_ascii=False))

        from modules import MODULE_REGISTRY
        cls = MODULE_REGISTRY.get(module_name)
        if not cls:
            raise ValueError(f"Unknown module: {module_name}")

        module = cls(project_dir=project_dir, params=params, progress_callback=progress_cb)

        # Best-effort input validation before running
        try:
            adata = module.load_adata(input_path)
            validation_error = module.validate_input(adata)
            if validation_error:
                raise ValueError(f"输入验证失败: {validation_error}")
        except FileNotFoundError:
            logger.debug(f"[Worker] Input file not found for validation: {input_path}")
        except ImportError as ie:
            logger.debug(f"[Worker] Module dependency missing for validation: {ie}")
        except (OSError, ValueError) as e:
            # 非 h5ad 文件（CSV/TSV/Excel）无法通过 load_adata 加载，
            # 但模块内部的 run() 可能使用 read_expression_matrix() 正确处理
            logger.debug(f"[Worker] Input validation skipped for {input_path}: {e}")

        result = module.run(input_path)

        register_task_outputs(task, project_id, project_dir, result)

        task.mark_completed(
            result.get('output_adata'),
            json.dumps(result.get('summary', {}))
        )

    except Exception as e:
        logger.error(f"[Worker] Task {task_id} failed:\n{traceback.format_exc()}")
        if task is None:
            return
        try:
            task.mark_failed(traceback.format_exc())
        except Exception as db_err:
            logger.error(f"[Worker] Failed to mark task {task_id} as failed: {db_err}")
    finally:
        _active_futures.pop(task_id, None)


def submit_pipeline_run(run_id, project_id, modules, params_by_module, project_dir, input_path):
    """提交 pipeline run 到线程池执行。"""
    if run_id in _active_futures:
        logger.warning(f"[Pipeline] Run {run_id} already submitted, skipping")
        return False
    future = _executor.submit(
        _run_pipeline_run, run_id, project_id, modules,
        params_by_module, project_dir, input_path
    )
    _active_futures[run_id] = future
    return True


def _write_pipeline_artifacts(run_id, project_dir):
    try:
        from models import PipelineRun, AnalysisTask, ResultFile
        from modules.reporting.pipeline_report import write_pipeline_report
        from modules.schemas import MODULE_DISPLAY_MAP
        run = PipelineRun.get_by_id(run_id)
        if not run:
            return
        try:
            task_ids = json.loads(run.task_ids_json) if run.task_ids_json else []
        except (TypeError, ValueError, json.JSONDecodeError):
            task_ids = []
        tasks = []
        for tid in task_ids:
            task = AnalysisTask.get_by_id(tid)
            if task:
                tasks.append(task)
        files_by_task = {task.id: ResultFile.get_by_task(task.id) for task in tasks}
        write_pipeline_report(project_dir, run, tasks, files_by_task, MODULE_DISPLAY_MAP)
    except Exception as e:
        logger.warning(f"[Pipeline] Failed to write pipeline artifacts for {run_id}: {e}")


def _run_pipeline_run(run_id, project_id, modules, params_by_module, project_dir, input_path):
    """同步执行 pipeline run：按顺序运行多个模块，串联 input/output。"""
    pipeline_run = None
    try:
        pipeline_run = PipelineRun.get_by_id(run_id)
        if not pipeline_run:
            logger.warning(f"[Pipeline] Run {run_id} not found")
            return

        if not pipeline_run.mark_running():
            logger.warning(f"[Pipeline] Run {run_id} not in pending state, skipping")
            return

        from modules import MODULE_REGISTRY

        current_input = input_path
        task_ids = []
        total_steps = len(modules)
        progress_log = []

        for i, module_name in enumerate(modules):
            # 更新 pipeline 进度
            pct = int((i / total_steps) * 100)
            pipeline_run.update_progress(pct, module_name, json.dumps(progress_log, ensure_ascii=False))

            # 创建子任务
            task = AnalysisTask(
                project_id=project_id,
                module_name=module_name,
                params_json=json.dumps(params_by_module.get(module_name, {}), ensure_ascii=False)
            )
            task.save()
            task_ids.append(task.id)

            # 更新 pipeline 的 task_ids
            pipeline_run.update_task_ids(task_ids)

            if not task.mark_running():
                logger.warning(f"[Pipeline] Task {task.id} for {module_name} not in pending state")
                pipeline_run.mark_failed(f"模块 {module_name} 的任务无法启动")
                _write_pipeline_artifacts(run_id, project_dir)
                return

            # 实例化模块
            cls = MODULE_REGISTRY.get(module_name)
            if not cls:
                task.mark_failed(f"未知模块: {module_name}")
                pipeline_run.mark_failed(f"未知模块: {module_name}")
                _write_pipeline_artifacts(run_id, project_dir)
                return

            progress_log_entry = {'step': i + 1, 'module': module_name, 'status': 'running'}

            def progress_cb(pct, message, _task_id=task.id, _module=module_name):
                now = datetime.now().strftime('%H:%M:%S')
                progress_log.append({'time': now, 'pct': pct, 'msg': f'[{_module}] {message}'})
                task.update_progress(pct, message, json.dumps(progress_log, ensure_ascii=False))

            module = cls(project_dir=project_dir, params=params_by_module.get(module_name, {}),
                         progress_callback=progress_cb)

            # 执行模块
            try:
                result = module.run(current_input)
            except Exception as e:
                tb = traceback.format_exc()
                task.mark_failed(tb)
                progress_log_entry['status'] = 'failed'
                progress_log_entry['error'] = str(e)
                progress_log.append(progress_log_entry)
                pipeline_run.mark_failed(f"模块 {module_name} 执行失败:\n{tb}")
                _write_pipeline_artifacts(run_id, project_dir)
                return

            # 检查输出
            output_adata = result.get('output_adata')
            if not output_adata:
                error_msg = f"模块 {module_name} 未返回 output_adata"
                task.mark_failed(error_msg)
                pipeline_run.mark_failed(error_msg)
                _write_pipeline_artifacts(run_id, project_dir)
                return

            register_task_outputs(task, project_id, project_dir, result, pipeline_run_id=run_id)

            # 标记任务完成
            task.mark_completed(output_adata, json.dumps(result.get('summary', {}), ensure_ascii=False))

            # 更新进度日志
            progress_log_entry['status'] = 'completed'
            progress_log.append(progress_log_entry)

            # 链式传递
            current_input = output_adata

        # 全部完成
        pipeline_run.mark_completed()
        _write_pipeline_artifacts(run_id, project_dir)

    except Exception as e:
        logger.error(f"[Pipeline] Run {run_id} failed:\n{traceback.format_exc()}")
        if pipeline_run:
            pipeline_run.mark_failed(traceback.format_exc())
            _write_pipeline_artifacts(run_id, project_dir)
    finally:
        _active_futures.pop(run_id, None)
