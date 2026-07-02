import sqlite3
from config import Config


def get_conn():
    """创建新的 SQLite 连接，设置 PRAGMA。每次调用返回新连接。"""
    conn = sqlite3.connect(Config.DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def init_db():
    db = get_conn()
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

        CREATE TABLE IF NOT EXISTS pipeline_runs (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL,
            name TEXT NOT NULL,
            analysis_type TEXT DEFAULT '',
            status TEXT DEFAULT 'pending',
            current_module TEXT DEFAULT '',
            progress INTEGER DEFAULT 0,
            input_path TEXT DEFAULT '',
            modules_json TEXT DEFAULT '[]',
            params_json TEXT DEFAULT '{}',
            task_ids_json TEXT DEFAULT '[]',
            error_traceback TEXT,
            started_at DATETIME,
            finished_at DATETIME,
            log_text TEXT DEFAULT '',
            FOREIGN KEY (project_id) REFERENCES projects(id) ON DELETE CASCADE
        );
    """)
    db.execute("CREATE INDEX IF NOT EXISTS idx_tasks_project ON analysis_tasks(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_tasks_status ON analysis_tasks(status)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_result_files_task ON result_files(task_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_result_files_project ON result_files(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_pipeline_runs_project ON pipeline_runs(project_id)")
    db.commit()
    db.close()
