import uuid
import json
import sqlite3
from datetime import datetime
from config import Config


def gen_id():
    return str(uuid.uuid4())[:8]


def _get_conn():
    conn = sqlite3.connect(Config.DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


class Project:
    def __init__(self, id=None, name='', description='', created_at=None, updated_at=None, status='empty', metadata_json='{}'):
        self.id = id or gen_id()
        self.name = name
        self.description = description
        self.created_at = created_at or datetime.now().isoformat()
        self.updated_at = updated_at or datetime.now().isoformat()
        self.status = status
        self.metadata_json = metadata_json

    def save(self):
        conn = _get_conn()
        conn.execute(
            "INSERT OR REPLACE INTO projects (id, name, description, created_at, updated_at, status, metadata_json) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (self.id, self.name, self.description, self.created_at, self.updated_at, self.status, self.metadata_json)
        )
        conn.commit()
        conn.close()

    def delete(self):
        conn = _get_conn()
        conn.execute("DELETE FROM projects WHERE id=?", (self.id,))
        conn.commit()
        conn.close()

    def to_dict(self):
        return {
            'id': self.id, 'name': self.name, 'description': self.description,
            'created_at': self.created_at, 'updated_at': self.updated_at,
            'status': self.status, 'metadata_json': self.metadata_json
        }

    @classmethod
    def get_by_id(cls, pid):
        conn = _get_conn()
        row = conn.execute("SELECT * FROM projects WHERE id=?", (pid,)).fetchone()
        conn.close()
        if row:
            return cls(**dict(row))
        return None

    @classmethod
    def get_all(cls):
        conn = _get_conn()
        rows = conn.execute("SELECT * FROM projects ORDER BY created_at DESC").fetchall()
        conn.close()
        return [cls(**dict(r)) for r in rows]

    def get_tasks(self):
        return AnalysisTask.get_by_project(self.id)

    def get_latest_adata_path(self):
        conn = _get_conn()
        row = conn.execute(
            "SELECT output_adata_path FROM analysis_tasks WHERE project_id=? AND status='completed' "
            "AND output_adata_path IS NOT NULL ORDER BY finished_at DESC LIMIT 1",
            (self.id,)
        ).fetchone()
        conn.close()
        return row['output_adata_path'] if row else None


class AnalysisTask:
    def __init__(self, id=None, project_id='', module_name='', status='pending',
                 progress=0, progress_message='', params_json='{}', result_json='{}',
                 error_traceback=None, started_at=None, finished_at=None,
                 output_adata_path=None, log_text=''):
        self.id = id or gen_id()
        self.project_id = project_id
        self.module_name = module_name
        self.status = status
        self.progress = progress
        self.progress_message = progress_message
        self.params_json = params_json
        self.result_json = result_json
        self.error_traceback = error_traceback
        self.started_at = started_at
        self.finished_at = finished_at
        self.output_adata_path = output_adata_path
        self.log_text = log_text

    def save(self):
        conn = _get_conn()
        conn.execute(
            "INSERT OR REPLACE INTO analysis_tasks "
            "(id, project_id, module_name, status, progress, progress_message, "
            "params_json, result_json, error_traceback, started_at, finished_at, "
            "output_adata_path, log_text) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (self.id, self.project_id, self.module_name, self.status, self.progress,
             self.progress_message, self.params_json, self.result_json,
             self.error_traceback, self.started_at, self.finished_at,
             self.output_adata_path, self.log_text)
        )
        conn.commit()
        conn.close()

    def to_dict(self):
        result = json.loads(self.result_json) if self.result_json else {}
        params = json.loads(self.params_json) if self.params_json else {}
        return {
            'id': self.id, 'project_id': self.project_id, 'module_name': self.module_name,
            'status': self.status, 'progress': self.progress,
            'progress_message': self.progress_message, 'params': params,
            'result': result, 'error_traceback': self.error_traceback,
            'started_at': self.started_at, 'finished_at': self.finished_at,
            'output_adata_path': self.output_adata_path
        }

    @classmethod
    def get_by_id(cls, tid):
        conn = _get_conn()
        row = conn.execute("SELECT * FROM analysis_tasks WHERE id=?", (tid,)).fetchone()
        conn.close()
        if row:
            return cls(**dict(row))
        return None

    @classmethod
    def get_by_project(cls, project_id):
        conn = _get_conn()
        rows = conn.execute(
            "SELECT * FROM analysis_tasks WHERE project_id=? ORDER BY started_at DESC",
            (project_id,)
        ).fetchall()
        conn.close()
        return [cls(**dict(r)) for r in rows]


class ResultFile:
    def __init__(self, id=None, task_id='', project_id='', file_type='',
                 category='', label='', file_path='', created_at=None):
        self.id = id or gen_id()
        self.task_id = task_id
        self.project_id = project_id
        self.file_type = file_type
        self.category = category
        self.label = label
        self.file_path = file_path
        self.created_at = created_at or datetime.now().isoformat()

    def save(self):
        conn = _get_conn()
        conn.execute(
            "INSERT INTO result_files (id, task_id, project_id, file_type, category, label, file_path, created_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (self.id, self.task_id, self.project_id, self.file_type, self.category,
             self.label, self.file_path, self.created_at)
        )
        conn.commit()
        conn.close()

    def to_dict(self):
        return {
            'id': self.id, 'task_id': self.task_id, 'project_id': self.project_id,
            'file_type': self.file_type, 'category': self.category,
            'label': self.label, 'file_path': self.file_path
        }

    @classmethod
    def get_by_task(cls, task_id):
        conn = _get_conn()
        rows = conn.execute("SELECT * FROM result_files WHERE task_id=?", (task_id,)).fetchall()
        conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def get_by_project(cls, project_id):
        conn = _get_conn()
        rows = conn.execute("SELECT * FROM result_files WHERE project_id=?", (project_id,)).fetchall()
        conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def get_by_id(cls, fid):
        conn = _get_conn()
        row = conn.execute("SELECT * FROM result_files WHERE id=?", (fid,)).fetchone()
        conn.close()
        if row:
            return cls(**dict(row))
        return None
