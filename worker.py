import concurrent.futures
import traceback
import json
import logging
from datetime import datetime
from database import get_conn
from models import gen_id, AnalysisTask, ResultFile

logger = logging.getLogger(__name__)

_executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
_active_futures = {}

def active_count():
    return sum(1 for f in _active_futures.values() if not f.done())

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

        result = module.run(input_path)

        for rf in (result.get('result_files') or []):
            try:
                ResultFile.create(
                    task_id=task_id, project_id=project_id,
                    file_type=rf.get('file_type', ''),
                    category=rf.get('category', ''),
                    label=rf.get('label', ''),
                    file_path=rf.get('file_path', '')
                )
            except Exception as e:
                logger.warning(f"[Worker] Skipping result_files insert: {rf.get('file_path', '')} ({e})")

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
