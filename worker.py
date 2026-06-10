import threading
import concurrent.futures
import traceback
import sys
import json
import sqlite3
from datetime import datetime
from config import Config
from models import gen_id

_executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
_active_futures = {}

def active_count():
    return sum(1 for f in _active_futures.values() if not f.done())

def submit_task(task_id, project_id, module_name, params, project_dir, input_path):
    future = _executor.submit(
        _run_task, task_id, project_id, module_name,
        params, project_dir, input_path
    )
    _active_futures[task_id] = future

def _run_task(task_id, project_id, module_name, params, project_dir, input_path):
    # Worker uses its own DB connection to avoid conflicts with Flask thread
    db = sqlite3.connect(Config.DB_PATH)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA foreign_keys=ON")
    try:
        db.execute(
            "UPDATE analysis_tasks SET status='running', started_at=CURRENT_TIMESTAMP WHERE id=?",
            (task_id,)
        )
        db.commit()

        progress_log = []

        def progress_cb(pct, message):
            now = datetime.now().strftime('%H:%M:%S')
            progress_log.append({'time': now, 'pct': pct, 'msg': message})
            db.execute(
                "UPDATE analysis_tasks SET progress=?, progress_message=?, log_text=? WHERE id=?",
                (pct, message, json.dumps(progress_log, ensure_ascii=False), task_id)
            )
            db.commit()

        from modules import MODULE_REGISTRY
        cls = MODULE_REGISTRY.get(module_name)
        if not cls:
            raise ValueError(f"Unknown module: {module_name}")

        module = cls(project_dir=project_dir, params=params, progress_callback=progress_cb)
        result = module.run(input_path)

        for rf in result.get('result_files', []):
            db.execute(
                "INSERT INTO result_files (id, task_id, project_id, file_type, category, label, file_path) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (gen_id(), task_id, project_id,
                 rf.get('file_type', ''), rf.get('category', ''), rf.get('label', ''), rf.get('file_path', ''))
            )

        db.execute(
            "UPDATE analysis_tasks SET status='completed', progress=100, "
            "progress_message='已完成', finished_at=CURRENT_TIMESTAMP, "
            "output_adata_path=?, result_json=? WHERE id=?",
            (result.get('output_adata'), json.dumps(result.get('summary', {})), task_id)
        )
        db.commit()

    except Exception as e:
        print(f"[Worker] Task {task_id} failed:\n{traceback.format_exc()}", file=sys.stderr)
        db.execute(
            "UPDATE analysis_tasks SET status='failed', error_traceback=?, "
            "progress_message='失败', finished_at=CURRENT_TIMESTAMP WHERE id=?",
            (traceback.format_exc(), task_id)
        )
        db.commit()
    finally:
        db.close()
        _active_futures.pop(task_id, None)
