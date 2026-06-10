import sqlite3
import os
from config import Config

_db = None

def get_db():
    global _db
    if _db is None:
        _db = sqlite3.connect(Config.DB_PATH, check_same_thread=False)
        _db.row_factory = sqlite3.Row
        _db.execute("PRAGMA journal_mode=WAL")
        _db.execute("PRAGMA foreign_keys=ON")
    return _db

def init_db():
    db = get_db()
    db.executescript("""
        CREATE TABLE IF NOT EXISTS projects (
            id TEXT PRIMARY KEY,
            name TEXT NOT NULL,
            description TEXT DEFAULT '',
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            status TEXT DEFAULT 'empty',
            metadata_json TEXT DEFAULT '{}'
        );

        CREATE TABLE IF NOT EXISTS analysis_tasks (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            module_name TEXT NOT NULL,
            status TEXT DEFAULT 'pending',
            progress INTEGER DEFAULT 0,
            progress_message TEXT DEFAULT '',
            params_json TEXT DEFAULT '{}',
            result_json TEXT DEFAULT '{}',
            error_traceback TEXT,
            started_at DATETIME,
            finished_at DATETIME,
            output_adata_path TEXT,
            log_text TEXT DEFAULT '',
            FOREIGN KEY (project_id) REFERENCES projects(id)
        );

        CREATE TABLE IF NOT EXISTS result_files (
            id TEXT PRIMARY KEY,
            task_id TEXT NOT NULL REFERENCES analysis_tasks(id) ON DELETE CASCADE,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            file_type TEXT NOT NULL,
            category TEXT NOT NULL,
            label TEXT NOT NULL,
            file_path TEXT NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (task_id) REFERENCES analysis_tasks(id),
            FOREIGN KEY (project_id) REFERENCES projects(id)
        );
    """)
    db.execute("CREATE INDEX IF NOT EXISTS idx_tasks_project ON analysis_tasks(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_tasks_status ON analysis_tasks(status)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_result_files_task ON result_files(task_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_result_files_project ON result_files(project_id)")
    db.commit()
