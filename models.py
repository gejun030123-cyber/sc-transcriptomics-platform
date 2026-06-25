import uuid
import json
from datetime import datetime, timezone
from database import get_conn


def gen_id():
    return str(uuid.uuid4())[:12]


class Project:
    def __init__(self, id=None, name='', description='', created_at=None, updated_at=None, status='empty', metadata_json='{}'):
        self.id = id or gen_id()
        self.name = name
        self.description = description
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()
        self.updated_at = updated_at or datetime.now(timezone.utc).isoformat()
        self.status = status
        self.metadata_json = metadata_json

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO projects (id, name, description, created_at, updated_at, status, metadata_json) "
                "VALUES (?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT(id) DO UPDATE SET name=excluded.name, description=excluded.description, "
                "created_at=excluded.created_at, updated_at=excluded.updated_at, "
                "status=excluded.status, metadata_json=excluded.metadata_json",
                (self.id, self.name, self.description, self.created_at, self.updated_at, self.status, self.metadata_json)
            )
            conn.commit()
        finally:
            conn.close()

    def delete(self):
        conn = get_conn()
        try:
            conn.execute("DELETE FROM projects WHERE id=?", (self.id,))
            conn.commit()
        finally:
            conn.close()

    def to_dict(self):
        return {
            'id': self.id, 'name': self.name, 'description': self.description,
            'created_at': self.created_at, 'updated_at': self.updated_at,
            'status': self.status, 'metadata_json': self.metadata_json
        }

    @classmethod
    def get_by_id(cls, pid):
        conn = get_conn()
        try:
            row = conn.execute("SELECT * FROM projects WHERE id=?", (pid,)).fetchone()
        finally:
            conn.close()
        if row:
            return cls(**dict(row))
        return None

    @classmethod
    def get_all(cls):
        conn = get_conn()
        try:
            rows = conn.execute("SELECT * FROM projects ORDER BY created_at DESC").fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]

    def get_tasks(self):
        return AnalysisTask.get_by_project(self.id)

    def get_latest_adata_path(self):
        conn = get_conn()
        try:
            row = conn.execute(
                "SELECT output_adata_path FROM analysis_tasks WHERE project_id=? AND status='completed' "
                "AND output_adata_path IS NOT NULL ORDER BY finished_at DESC LIMIT 1",
                (self.id,)
            ).fetchone()
        finally:
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
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO analysis_tasks "
                "(id, project_id, module_name, status, progress, progress_message, "
                "params_json, result_json, error_traceback, started_at, finished_at, "
                "output_adata_path, log_text) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET project_id=excluded.project_id, "
                "module_name=excluded.module_name, status=excluded.status, "
                "progress=excluded.progress, progress_message=excluded.progress_message, "
                "params_json=excluded.params_json, result_json=excluded.result_json, "
                "error_traceback=excluded.error_traceback, started_at=excluded.started_at, "
                "finished_at=excluded.finished_at, output_adata_path=excluded.output_adata_path, "
                "log_text=excluded.log_text",
                (self.id, self.project_id, self.module_name, self.status, self.progress,
                 self.progress_message, self.params_json, self.result_json,
                 self.error_traceback, self.started_at, self.finished_at,
                 self.output_adata_path, self.log_text)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_running(self):
        conn = get_conn()
        try:
            cursor = conn.execute(
                "UPDATE analysis_tasks SET status='running', started_at=CURRENT_TIMESTAMP "
                "WHERE id=? AND status='pending'",
                (self.id,)
            )
            conn.commit()
            return cursor.rowcount > 0
        finally:
            conn.close()

    def update_progress(self, pct, message, log_text=None):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE analysis_tasks SET progress=?, progress_message=?, log_text=? WHERE id=?",
                (pct, message, log_text, self.id)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_completed(self, output_adata, result_json):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE analysis_tasks SET status='completed', progress=100, "
                "progress_message='已完成', finished_at=CURRENT_TIMESTAMP, "
                "output_adata_path=?, result_json=? WHERE id=? AND status='running'",
                (output_adata, result_json, self.id)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_failed(self, error_traceback):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE analysis_tasks SET status='failed', error_traceback=?, "
                "progress_message='失败', finished_at=CURRENT_TIMESTAMP "
                "WHERE id=? AND status IN ('pending', 'running')",
                (error_traceback, self.id)
            )
            conn.commit()
        finally:
            conn.close()

    def to_dict(self):
        try:
            result = json.loads(self.result_json) if self.result_json else {}
        except (json.JSONDecodeError, ValueError):
            result = {}
        try:
            params = json.loads(self.params_json) if self.params_json else {}
        except (json.JSONDecodeError, ValueError):
            params = {}
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
        conn = get_conn()
        try:
            row = conn.execute("SELECT * FROM analysis_tasks WHERE id=?", (tid,)).fetchone()
        finally:
            conn.close()
        if row:
            return cls(**dict(row))
        return None

    @classmethod
    def get_by_project(cls, project_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM analysis_tasks WHERE project_id=? ORDER BY started_at DESC",
                (project_id,)
            ).fetchall()
        finally:
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
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO result_files (id, task_id, project_id, file_type, category, label, file_path, created_at) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT(id) DO UPDATE SET task_id=excluded.task_id, "
                "project_id=excluded.project_id, file_type=excluded.file_type, "
                "category=excluded.category, label=excluded.label, "
                "file_path=excluded.file_path, created_at=excluded.created_at",
                (self.id, self.task_id, self.project_id, self.file_type, self.category,
                 self.label, self.file_path, self.created_at)
            )
            conn.commit()
        finally:
            conn.close()

    @classmethod
    def create(cls, task_id, project_id, file_type, category, label, file_path):
        rf = cls(task_id=task_id, project_id=project_id,
                 file_type=file_type, category=category, label=label, file_path=file_path)
        rf.save()
        return rf

    def to_dict(self):
        return {
            'id': self.id, 'task_id': self.task_id, 'project_id': self.project_id,
            'file_type': self.file_type, 'category': self.category,
            'label': self.label, 'file_path': self.file_path
        }

    @classmethod
    def get_by_task(cls, task_id):
        conn = get_conn()
        try:
            rows = conn.execute("SELECT * FROM result_files WHERE task_id=?", (task_id,)).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def get_by_project(cls, project_id):
        conn = get_conn()
        try:
            rows = conn.execute("SELECT * FROM result_files WHERE project_id=?", (project_id,)).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def get_by_id(cls, fid):
        conn = get_conn()
        try:
            row = conn.execute("SELECT * FROM result_files WHERE id=?", (fid,)).fetchone()
        finally:
            conn.close()
        if row:
            return cls(**dict(row))
        return None
