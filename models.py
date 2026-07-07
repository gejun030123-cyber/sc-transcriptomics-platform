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
        """获取项目当前最新有效 h5ad 路径。

        优先级：
        1. 如果已采纳 branch，返回该 branch 的 output_adata_path
        2. 否则返回最新主线任务（branch_id IS NULL）的输出
        """
        conn = get_conn()
        try:
            # 检查是否有已采纳 branch
            meta_row = conn.execute(
                "SELECT metadata_json FROM projects WHERE id=?", (self.id,)
            ).fetchone()
            if meta_row:
                try:
                    meta = json.loads(meta_row['metadata_json'] or '{}')
                except Exception:
                    meta = {}
                accepted_path = meta.get('current_adata_path', '')
                accepted_branch_id = meta.get('accepted_branch_id', '')
                if accepted_path and accepted_branch_id:
                    # 校验 branch 仍存在且为 accepted 状态
                    branch_row = conn.execute(
                        "SELECT id FROM analysis_branches WHERE id=? AND accepted=1 AND deleted=0",
                        (accepted_branch_id,)
                    ).fetchone()
                    if branch_row:
                        import os
                        if os.path.isfile(accepted_path):
                            return accepted_path

            # 回退：最新主线任务
            row = conn.execute(
                "SELECT output_adata_path FROM analysis_tasks "
                "WHERE project_id=? AND status='completed' "
                "AND branch_id IS NULL "
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
                 output_adata_path=None, log_text='', branch_id=None):
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
        self.branch_id = branch_id

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO analysis_tasks "
                "(id, project_id, module_name, status, progress, progress_message, "
                "params_json, result_json, error_traceback, started_at, finished_at, "
                "output_adata_path, log_text, branch_id) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET project_id=excluded.project_id, "
                "module_name=excluded.module_name, status=excluded.status, "
                "progress=excluded.progress, progress_message=excluded.progress_message, "
                "params_json=excluded.params_json, result_json=excluded.result_json, "
                "error_traceback=excluded.error_traceback, started_at=excluded.started_at, "
                "finished_at=excluded.finished_at, output_adata_path=excluded.output_adata_path, "
                "log_text=excluded.log_text, branch_id=excluded.branch_id",
                (self.id, self.project_id, self.module_name, self.status, self.progress,
                 self.progress_message, self.params_json, self.result_json,
                 self.error_traceback, self.started_at, self.finished_at,
                 self.output_adata_path, self.log_text, self.branch_id)
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
            'output_adata_path': self.output_adata_path,
            'branch_id': self.branch_id,
            'is_candidate': bool(self.branch_id),
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
    def get_by_project(cls, project_id, include_branches=False):
        conn = get_conn()
        try:
            if include_branches:
                rows = conn.execute(
                    "SELECT * FROM analysis_tasks WHERE project_id=? ORDER BY started_at DESC",
                    (project_id,)
                ).fetchall()
            else:
                rows = conn.execute(
                    "SELECT * FROM analysis_tasks WHERE project_id=? AND branch_id IS NULL ORDER BY started_at DESC",
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


class PipelineRun:
    def __init__(self, id=None, project_id='', name='', analysis_type='',
                 status='pending', current_module='', progress=0, input_path='',
                 modules_json='[]', params_json='{}', task_ids_json='[]',
                 error_traceback=None, started_at=None, finished_at=None,
                 log_text=''):
        self.id = id or gen_id()
        self.project_id = project_id
        self.name = name
        self.analysis_type = analysis_type
        self.status = status
        self.current_module = current_module
        self.progress = progress
        self.input_path = input_path
        self.modules_json = modules_json
        self.params_json = params_json
        self.task_ids_json = task_ids_json
        self.error_traceback = error_traceback
        self.started_at = started_at
        self.finished_at = finished_at
        self.log_text = log_text

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO pipeline_runs "
                "(id, project_id, name, analysis_type, status, current_module, progress, "
                "input_path, modules_json, params_json, task_ids_json, error_traceback, "
                "started_at, finished_at, log_text) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET project_id=excluded.project_id, "
                "name=excluded.name, analysis_type=excluded.analysis_type, "
                "status=excluded.status, current_module=excluded.current_module, "
                "progress=excluded.progress, input_path=excluded.input_path, "
                "modules_json=excluded.modules_json, params_json=excluded.params_json, "
                "task_ids_json=excluded.task_ids_json, error_traceback=excluded.error_traceback, "
                "started_at=excluded.started_at, finished_at=excluded.finished_at, "
                "log_text=excluded.log_text",
                (self.id, self.project_id, self.name, self.analysis_type, self.status,
                 self.current_module, self.progress, self.input_path, self.modules_json,
                 self.params_json, self.task_ids_json, self.error_traceback,
                 self.started_at, self.finished_at, self.log_text)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_running(self):
        conn = get_conn()
        try:
            cursor = conn.execute(
                "UPDATE pipeline_runs SET status='running', started_at=CURRENT_TIMESTAMP "
                "WHERE id=? AND status='pending'",
                (self.id,)
            )
            conn.commit()
            return cursor.rowcount > 0
        finally:
            conn.close()

    def update_progress(self, pct, current_module='', log_text=None):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE pipeline_runs SET progress=?, current_module=?, log_text=? WHERE id=?",
                (pct, current_module, log_text, self.id)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_completed(self):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE pipeline_runs SET status='completed', progress=100, "
                "finished_at=CURRENT_TIMESTAMP WHERE id=? AND status='running'",
                (self.id,)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_failed(self, error_traceback):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE pipeline_runs SET status='failed', error_traceback=?, "
                "finished_at=CURRENT_TIMESTAMP WHERE id=? AND status IN ('pending', 'running')",
                (error_traceback, self.id)
            )
            conn.commit()
        finally:
            conn.close()

    def to_dict(self):
        try:
            modules = json.loads(self.modules_json) if self.modules_json else []
        except (json.JSONDecodeError, ValueError):
            modules = []
        try:
            params = json.loads(self.params_json) if self.params_json else {}
        except (json.JSONDecodeError, ValueError):
            params = {}
        try:
            task_ids = json.loads(self.task_ids_json) if self.task_ids_json else []
        except (json.JSONDecodeError, ValueError):
            task_ids = []
        return {
            'id': self.id, 'project_id': self.project_id, 'name': self.name,
            'analysis_type': self.analysis_type, 'status': self.status,
            'current_module': self.current_module, 'progress': self.progress,
            'input_path': self.input_path, 'modules': modules, 'params': params,
            'task_ids': task_ids, 'error_traceback': self.error_traceback,
            'started_at': self.started_at, 'finished_at': self.finished_at
        }

    @classmethod
    def get_by_id(cls, rid):
        conn = get_conn()
        try:
            row = conn.execute("SELECT * FROM pipeline_runs WHERE id=?", (rid,)).fetchone()
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
                "SELECT * FROM pipeline_runs WHERE project_id=? ORDER BY started_at DESC",
                (project_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]


# ============================================================
# Agent 系统模型（阶段 1）
# ============================================================

class AgentSession:
    """一次用户与 agent 的目标分析会话"""
    def __init__(self, id=None, project_id='', status='active',
                 goal_summary='', created_at=None, updated_at=None):
        self.id = id or gen_id()
        self.project_id = project_id
        self.status = status
        self.goal_summary = goal_summary
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()
        self.updated_at = updated_at or datetime.now(timezone.utc).isoformat()

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO agent_sessions (id, project_id, status, goal_summary, "
                "created_at, updated_at) VALUES (?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET status=excluded.status, "
                "goal_summary=excluded.goal_summary, updated_at=excluded.updated_at",
                (self.id, self.project_id, self.status, self.goal_summary,
                 self.created_at, self.updated_at)
            )
            conn.commit()
        finally:
            conn.close()

    def to_dict(self):
        return {
            'id': self.id, 'project_id': self.project_id, 'status': self.status,
            'goal_summary': self.goal_summary,
            'created_at': self.created_at, 'updated_at': self.updated_at
        }

    @classmethod
    def get_by_id(cls, sid):
        conn = get_conn()
        try:
            row = conn.execute("SELECT * FROM agent_sessions WHERE id=?", (sid,)).fetchone()
        finally:
            conn.close()
        return cls(**dict(row)) if row else None

    @classmethod
    def get_by_project(cls, project_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM agent_sessions WHERE project_id=? ORDER BY created_at DESC",
                (project_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]


class AgentGoal:
    """结构化目标，例如"寻找 microglia-like cluster"."""
    def __init__(self, id=None, session_id='', project_id='',
                 goal_type='', target_entity='', goal_json='{}',
                 constraints_json='{}', status='draft',
                 created_at=None, updated_at=None):
        self.id = id or gen_id()
        self.session_id = session_id
        self.project_id = project_id
        self.goal_type = goal_type
        self.target_entity = target_entity
        self.goal_json = goal_json
        self.constraints_json = constraints_json
        self.status = status
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()
        self.updated_at = updated_at or datetime.now(timezone.utc).isoformat()

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO agent_goals (id, session_id, project_id, goal_type, "
                "target_entity, goal_json, constraints_json, status, created_at, updated_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET "
                "goal_type=excluded.goal_type, target_entity=excluded.target_entity, "
                "goal_json=excluded.goal_json, constraints_json=excluded.constraints_json, "
                "status=excluded.status, updated_at=excluded.updated_at",
                (self.id, self.session_id, self.project_id, self.goal_type,
                 self.target_entity, self.goal_json, self.constraints_json,
                 self.status, self.created_at, self.updated_at)
            )
            conn.commit()
        finally:
            conn.close()

    def to_dict(self):
        try:
            goal = json.loads(self.goal_json) if self.goal_json else {}
        except (json.JSONDecodeError, ValueError):
            goal = {}
        try:
            constraints = json.loads(self.constraints_json) if self.constraints_json else {}
        except (json.JSONDecodeError, ValueError):
            constraints = {}
        return {
            'id': self.id, 'session_id': self.session_id, 'project_id': self.project_id,
            'goal_type': self.goal_type, 'target_entity': self.target_entity,
            'goal': goal, 'constraints': constraints,
            'status': self.status, 'created_at': self.created_at, 'updated_at': self.updated_at
        }

    @classmethod
    def get_by_id(cls, gid):
        conn = get_conn()
        try:
            row = conn.execute("SELECT * FROM agent_goals WHERE id=?", (gid,)).fetchone()
        finally:
            conn.close()
        return cls(**dict(row)) if row else None

    @classmethod
    def get_by_session(cls, session_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM agent_goals WHERE session_id=? ORDER BY created_at DESC",
                (session_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]


class AgentStep:
    """记录 agent 每一步的观察、假设、动作和结果."""
    def __init__(self, id=None, session_id='', goal_id=None,
                 step_index=0, step_type='', thought_summary='',
                 action_name='', action_args_json='{}',
                 result_json='{}', status='completed', created_at=None):
        self.id = id or gen_id()
        self.session_id = session_id
        self.goal_id = goal_id
        self.step_index = step_index
        self.step_type = step_type
        self.thought_summary = thought_summary
        self.action_name = action_name
        self.action_args_json = action_args_json
        self.result_json = result_json
        self.status = status
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO agent_steps (id, session_id, goal_id, step_index, "
                "step_type, thought_summary, action_name, action_args_json, "
                "result_json, status, created_at) VALUES (?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET "
                "step_index=excluded.step_index, step_type=excluded.step_type, "
                "thought_summary=excluded.thought_summary, action_name=excluded.action_name, "
                "action_args_json=excluded.action_args_json, result_json=excluded.result_json, "
                "status=excluded.status",
                (self.id, self.session_id, self.goal_id, self.step_index,
                 self.step_type, self.thought_summary, self.action_name,
                 self.action_args_json, self.result_json, self.status, self.created_at)
            )
            conn.commit()
        finally:
            conn.close()

    def to_dict(self):
        try:
            args = json.loads(self.action_args_json) if self.action_args_json else {}
        except (json.JSONDecodeError, ValueError):
            args = {}
        try:
            result = json.loads(self.result_json) if self.result_json else {}
        except (json.JSONDecodeError, ValueError):
            result = {}
        return {
            'id': self.id, 'session_id': self.session_id, 'goal_id': self.goal_id,
            'step_index': self.step_index, 'step_type': self.step_type,
            'thought_summary': self.thought_summary, 'action_name': self.action_name,
            'action_args': args, 'result': result,
            'status': self.status, 'created_at': self.created_at
        }

    @classmethod
    def get_by_session(cls, session_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM agent_steps WHERE session_id=? ORDER BY step_index ASC",
                (session_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def get_by_goal(cls, goal_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM agent_steps WHERE goal_id=? ORDER BY step_index ASC",
                (goal_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def next_index(cls, session_id):
        """获取 session 的下一个 step_index."""
        conn = get_conn()
        try:
            row = conn.execute(
                "SELECT COALESCE(MAX(step_index), -1) + 1 AS next_idx "
                "FROM agent_steps WHERE session_id=?", (session_id,)
            ).fetchone()
        finally:
            conn.close()
        return row['next_idx'] if row else 0


class AnalysisBranch:
    """候选分析分支，避免覆盖主线结果。软删除通过 deleted=1 标记。"""
    def __init__(self, id=None, project_id='', session_id=None, goal_id=None,
                 parent_task_id=None, parent_adata_path='', branch_name='',
                 purpose='', params_json='{}', output_adata_path=None,
                 status='pending', accepted=0, deleted=0, error_traceback=None,
                 created_at=None, finished_at=None):
        self.id = id or gen_id()
        self.project_id = project_id
        self.session_id = session_id
        self.goal_id = goal_id
        self.parent_task_id = parent_task_id
        self.parent_adata_path = parent_adata_path
        self.branch_name = branch_name
        self.purpose = purpose
        self.params_json = params_json
        self.output_adata_path = output_adata_path
        self.status = status
        self.accepted = accepted
        self.deleted = deleted
        self.error_traceback = error_traceback
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()
        self.finished_at = finished_at

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO analysis_branches (id, project_id, session_id, goal_id, "
                "parent_task_id, parent_adata_path, branch_name, purpose, params_json, "
                "output_adata_path, status, accepted, deleted, error_traceback, created_at, finished_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET "
                "status=excluded.status, accepted=excluded.accepted, "
                "deleted=excluded.deleted, "
                "output_adata_path=excluded.output_adata_path, "
                "error_traceback=excluded.error_traceback, "
                "finished_at=excluded.finished_at",
                (self.id, self.project_id, self.session_id, self.goal_id,
                 self.parent_task_id, self.parent_adata_path, self.branch_name,
                 self.purpose, self.params_json, self.output_adata_path,
                 self.status, self.accepted, self.deleted, self.error_traceback,
                 self.created_at, self.finished_at)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_running(self):
        conn = get_conn()
        try:
            cursor = conn.execute(
                "UPDATE analysis_branches SET status='running' "
                "WHERE id=? AND status='pending'", (self.id,)
            )
            conn.commit()
            return cursor.rowcount > 0
        finally:
            conn.close()

    def mark_completed(self, output_adata_path):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE analysis_branches SET status='completed', "
                "output_adata_path=?, finished_at=CURRENT_TIMESTAMP "
                "WHERE id=? AND status='running'",
                (output_adata_path, self.id)
            )
            conn.commit()
            self.status = 'completed'
            self.output_adata_path = output_adata_path
        finally:
            conn.close()

    def mark_failed(self, error_traceback):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE analysis_branches SET status='failed', error_traceback=?, "
                "finished_at=CURRENT_TIMESTAMP WHERE id=? AND status IN ('pending', 'running')",
                (error_traceback, self.id)
            )
            conn.commit()
        finally:
            conn.close()

    def accept(self):
        """采纳此 branch，同时更新项目 metadata 记录 accepted branch."""
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE analysis_branches SET accepted=1 WHERE id=?", (self.id,)
            )
            # 更新项目 metadata 记录当前采纳的 branch
            if self.output_adata_path:
                import json as _json
                proj_row = conn.execute(
                    "SELECT metadata_json FROM projects WHERE id=?", (self.project_id,)
                ).fetchone()
                if proj_row:
                    try:
                        meta = _json.loads(proj_row['metadata_json'] or '{}')
                    except Exception:
                        meta = {}
                    meta['accepted_branch_id'] = self.id
                    meta['current_adata_path'] = self.output_adata_path
                    conn.execute(
                        "UPDATE projects SET metadata_json=? WHERE id=?",
                        (_json.dumps(meta, ensure_ascii=False), self.project_id)
                    )
            conn.commit()
            self.accepted = 1
        finally:
            conn.close()

    def soft_delete(self):
        """软删除此 branch（设置 deleted=1）。不会触发 FK 级联删除任务。"""
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE analysis_branches SET deleted=1 WHERE id=?", (self.id,)
            )
            conn.commit()
            self.deleted = 1
        finally:
            conn.close()

    def to_dict(self):
        try:
            params = json.loads(self.params_json) if self.params_json else {}
        except (json.JSONDecodeError, ValueError):
            params = {}
        return {
            'id': self.id, 'project_id': self.project_id,
            'session_id': self.session_id, 'goal_id': self.goal_id,
            'parent_task_id': self.parent_task_id, 'parent_adata_path': self.parent_adata_path,
            'branch_name': self.branch_name, 'purpose': self.purpose,
            'params': params, 'output_adata_path': self.output_adata_path,
            'status': self.status, 'accepted': self.accepted, 'deleted': self.deleted,
            'error_traceback': self.error_traceback,
            'created_at': self.created_at, 'finished_at': self.finished_at
        }

    @classmethod
    def get_by_id(cls, bid, include_deleted=False):
        conn = get_conn()
        try:
            if include_deleted:
                row = conn.execute("SELECT * FROM analysis_branches WHERE id=?", (bid,)).fetchone()
            else:
                row = conn.execute(
                    "SELECT * FROM analysis_branches WHERE id=? AND deleted=0", (bid,)
                ).fetchone()
        finally:
            conn.close()
        return cls(**dict(row)) if row else None

    @classmethod
    def get_by_project(cls, project_id, include_deleted=False):
        conn = get_conn()
        try:
            if include_deleted:
                rows = conn.execute(
                    "SELECT * FROM analysis_branches WHERE project_id=? ORDER BY created_at DESC",
                    (project_id,)
                ).fetchall()
            else:
                rows = conn.execute(
                    "SELECT * FROM analysis_branches WHERE project_id=? AND deleted=0 ORDER BY created_at DESC",
                    (project_id,)
                ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def get_by_goal(cls, goal_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM analysis_branches WHERE goal_id=? ORDER BY created_at DESC",
                (goal_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]


class CandidateScore:
    """候选结果评分."""
    def __init__(self, id=None, branch_id='', project_id='',
                 evaluator_name='', target_label='', score_json='{}',
                 total_score=None, recommendation='', created_at=None):
        self.id = id or gen_id()
        self.branch_id = branch_id
        self.project_id = project_id
        self.evaluator_name = evaluator_name
        self.target_label = target_label
        self.score_json = score_json
        self.total_score = total_score
        self.recommendation = recommendation
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO candidate_scores (id, branch_id, project_id, "
                "evaluator_name, target_label, score_json, total_score, "
                "recommendation, created_at) VALUES (?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET "
                "total_score=excluded.total_score, recommendation=excluded.recommendation",
                (self.id, self.branch_id, self.project_id, self.evaluator_name,
                 self.target_label, self.score_json, self.total_score,
                 self.recommendation, self.created_at)
            )
            conn.commit()
        finally:
            conn.close()

    def to_dict(self):
        try:
            score = json.loads(self.score_json) if self.score_json else {}
        except (json.JSONDecodeError, ValueError):
            score = {}
        return {
            'id': self.id, 'branch_id': self.branch_id, 'project_id': self.project_id,
            'evaluator_name': self.evaluator_name, 'target_label': self.target_label,
            'score': score, 'total_score': self.total_score,
            'recommendation': self.recommendation, 'created_at': self.created_at
        }

    @classmethod
    def get_by_branch(cls, branch_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM candidate_scores WHERE branch_id=? ORDER BY total_score DESC",
                (branch_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]

    @classmethod
    def get_by_project(cls, project_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM candidate_scores WHERE project_id=? ORDER BY total_score DESC",
                (project_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]


class AgentJob:
    """异步 Agent 任务（参数 sweep、branch run 等）。"""
    def __init__(self, id=None, project_id='', session_id=None, goal_id=None,
                 job_type='', status='pending', progress=0, current_step='',
                 params_json='{}', result_json='{}', error_traceback=None,
                 created_at=None, started_at=None, finished_at=None):
        self.id = id or gen_id()
        self.project_id = project_id
        self.session_id = session_id
        self.goal_id = goal_id
        self.job_type = job_type
        self.status = status
        self.progress = progress
        self.current_step = current_step
        self.params_json = params_json
        self.result_json = result_json
        self.error_traceback = error_traceback
        self.created_at = created_at or datetime.now(timezone.utc).isoformat()
        self.started_at = started_at
        self.finished_at = finished_at

    def save(self):
        conn = get_conn()
        try:
            conn.execute(
                "INSERT INTO agent_jobs (id, project_id, session_id, goal_id, "
                "job_type, status, progress, current_step, params_json, result_json, "
                "error_traceback, created_at, started_at, finished_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET "
                "status=excluded.status, progress=excluded.progress, "
                "current_step=excluded.current_step, result_json=excluded.result_json, "
                "error_traceback=excluded.error_traceback, finished_at=excluded.finished_at",
                (self.id, self.project_id, self.session_id, self.goal_id,
                 self.job_type, self.status, self.progress, self.current_step,
                 self.params_json, self.result_json, self.error_traceback,
                 self.created_at, self.started_at, self.finished_at)
            )
            conn.commit()
        finally:
            conn.close()

    def mark_running(self):
        conn = get_conn()
        try:
            cursor = conn.execute(
                "UPDATE agent_jobs SET status='running', started_at=CURRENT_TIMESTAMP "
                "WHERE id=? AND status='pending'", (self.id,)
            )
            conn.commit()
            self.status = 'running'
            return cursor.rowcount > 0
        finally:
            conn.close()

    def update_progress(self, pct, current_step, result_json=None):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE agent_jobs SET progress=?, current_step=?, result_json=? WHERE id=?",
                (pct, current_step, result_json, self.id)
            )
            conn.commit()
            self.progress = pct
            self.current_step = current_step
        finally:
            conn.close()

    def mark_completed(self, result_json):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE agent_jobs SET status='completed', progress=100, "
                "result_json=?, finished_at=CURRENT_TIMESTAMP "
                "WHERE id=? AND status='running'",
                (result_json, self.id)
            )
            conn.commit()
            self.status = 'completed'
        finally:
            conn.close()

    def mark_failed(self, error_traceback):
        conn = get_conn()
        try:
            conn.execute(
                "UPDATE agent_jobs SET status='failed', error_traceback=?, "
                "finished_at=CURRENT_TIMESTAMP WHERE id=? AND status IN ('pending', 'running')",
                (error_traceback, self.id)
            )
            conn.commit()
            self.status = 'failed'
        finally:
            conn.close()

    def to_dict(self):
        try:
            params = json.loads(self.params_json) if self.params_json else {}
        except Exception:
            params = {}
        try:
            result = json.loads(self.result_json) if self.result_json else {}
        except Exception:
            result = {}
        return {
            'id': self.id, 'project_id': self.project_id,
            'session_id': self.session_id, 'goal_id': self.goal_id,
            'job_type': self.job_type, 'status': self.status,
            'progress': self.progress, 'current_step': self.current_step,
            'params': params, 'result': result,
            'error_traceback': self.error_traceback,
            'created_at': self.created_at, 'started_at': self.started_at,
            'finished_at': self.finished_at,
        }

    @classmethod
    def get_by_id(cls, jid):
        conn = get_conn()
        try:
            row = conn.execute("SELECT * FROM agent_jobs WHERE id=?", (jid,)).fetchone()
        finally:
            conn.close()
        return cls(**dict(row)) if row else None

    @classmethod
    def get_by_project(cls, project_id):
        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT * FROM agent_jobs WHERE project_id=? ORDER BY created_at DESC",
                (project_id,)
            ).fetchall()
        finally:
            conn.close()
        return [cls(**dict(r)) for r in rows]
