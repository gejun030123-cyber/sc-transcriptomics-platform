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
            branch_id TEXT,
            FOREIGN KEY (project_id) REFERENCES projects(id),
            FOREIGN KEY (branch_id) REFERENCES analysis_branches(id) ON DELETE SET NULL
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

        -- Platform-level overrides entered from the AI API settings page.
        -- Values are never returned to the browser in clear text; the API key
        -- is only read by the server-side adapter at request time.
        CREATE TABLE IF NOT EXISTS platform_settings (
            key TEXT PRIMARY KEY,
            value TEXT DEFAULT '',
            updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
        );
    """)
    db.execute("CREATE INDEX IF NOT EXISTS idx_tasks_project ON analysis_tasks(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_tasks_status ON analysis_tasks(status)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_result_files_task ON result_files(task_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_result_files_project ON result_files(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_pipeline_runs_project ON pipeline_runs(project_id)")

    # Generic assay/workflow assets.  These tables are deliberately additive:
    # existing AnalysisTask/PipelineRun records remain the compatibility view
    # for single-cell and Bulk RNA modules.
    db.executescript("""
        CREATE TABLE IF NOT EXISTS data_assets (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            sample_id TEXT DEFAULT '',
            assay_type TEXT NOT NULL DEFAULT '',
            role TEXT DEFAULT '',
            artifact_kind TEXT NOT NULL,
            format TEXT DEFAULT '',
            file_path TEXT NOT NULL,
            size_bytes INTEGER DEFAULT 0,
            checksum TEXT DEFAULT '',
            reference_build TEXT DEFAULT '',
            metadata_json TEXT DEFAULT '{}',
            parent_asset_ids_json TEXT DEFAULT '[]',
            created_by_run_id TEXT DEFAULT '',
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS sample_manifests (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            version INTEGER NOT NULL,
            content_json TEXT NOT NULL,
            checksum TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'draft',
            validation_json TEXT DEFAULT '{}',
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            UNIQUE(project_id, version),
            UNIQUE(project_id, checksum)
        );

        CREATE TABLE IF NOT EXISTS reference_assets (
            id TEXT PRIMARY KEY,
            species TEXT DEFAULT '',
            assembly TEXT NOT NULL,
            bundle_version TEXT NOT NULL,
            asset_type TEXT NOT NULL,
            file_path TEXT NOT NULL,
            checksum TEXT DEFAULT '',
            source TEXT DEFAULT '',
            license_note TEXT DEFAULT '',
            metadata_json TEXT DEFAULT '{}',
            status TEXT NOT NULL DEFAULT 'registered',
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS workflow_runs (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            workflow_key TEXT NOT NULL,
            manifest_id TEXT NOT NULL REFERENCES sample_manifests(id),
            status TEXT NOT NULL DEFAULT 'pending',
            external_run_id TEXT DEFAULT '',
            parameter_signature TEXT DEFAULT '',
            launch_json TEXT DEFAULT '{}',
            provenance_json TEXT DEFAULT '{}',
            error_text TEXT DEFAULT '',
            executor TEXT DEFAULT 'nextflow',
            workflow_release TEXT DEFAULT '',
            profile TEXT DEFAULT '',
            run_dir TEXT DEFAULT '',
            pid INTEGER,
            process_group_id INTEGER,
            stdout_path TEXT DEFAULT '',
            stderr_path TEXT DEFAULT '',
            exit_code INTEGER,
            cancel_requested_at DATETIME,
            started_at DATETIME,
            finished_at DATETIME,
            updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS workflow_artifacts (
            id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL REFERENCES workflow_runs(id) ON DELETE CASCADE,
            asset_id TEXT REFERENCES data_assets(id) ON DELETE SET NULL,
            artifact_kind TEXT NOT NULL,
            file_path TEXT NOT NULL,
            checksum TEXT DEFAULT '',
            required INTEGER NOT NULL DEFAULT 0,
            metadata_json TEXT DEFAULT '{}',
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS wes_sra_jobs (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            sample_id TEXT NOT NULL,
            patient_id TEXT NOT NULL,
            role TEXT NOT NULL,
            capture_kit_id TEXT DEFAULT '',
            source_path TEXT NOT NULL,
            source_asset_id TEXT DEFAULT '',
            output_dir TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'pending',
            progress INTEGER NOT NULL DEFAULT 0,
            progress_message TEXT DEFAULT '',
            input_size_bytes INTEGER NOT NULL DEFAULT 0,
            output_paths_json TEXT DEFAULT '{}',
            output_asset_ids_json TEXT DEFAULT '[]',
            log_text TEXT DEFAULT '',
            error_text TEXT DEFAULT '',
            started_at DATETIME,
            finished_at DATETIME,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        );
    """)
    db.execute("CREATE INDEX IF NOT EXISTS idx_assets_project ON data_assets(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_manifests_project ON sample_manifests(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_workflow_runs_project ON workflow_runs(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_workflow_artifacts_run ON workflow_artifacts(run_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_wes_sra_jobs_project ON wes_sra_jobs(project_id)")

    wes_sra_columns = {
        row["name"] for row in db.execute("PRAGMA table_info(wes_sra_jobs)").fetchall()
    }
    if "capture_kit_id" not in wes_sra_columns:
        db.execute("ALTER TABLE wes_sra_jobs ADD COLUMN capture_kit_id TEXT DEFAULT ''")

    # Existing lab databases may have been created before the lightweight WES
    # runtime fields were introduced.  Keep this additive and migration-free:
    # only fixed, trusted column names are used here.
    workflow_run_columns = {
        row["name"] for row in db.execute("PRAGMA table_info(workflow_runs)").fetchall()
    }
    workflow_run_additions = {
        "executor": "TEXT DEFAULT 'nextflow'",
        "workflow_release": "TEXT DEFAULT ''",
        "profile": "TEXT DEFAULT ''",
        "run_dir": "TEXT DEFAULT ''",
        "pid": "INTEGER",
        "process_group_id": "INTEGER",
        "stdout_path": "TEXT DEFAULT ''",
        "stderr_path": "TEXT DEFAULT ''",
        "exit_code": "INTEGER",
        "cancel_requested_at": "DATETIME",
        "updated_at": "DATETIME DEFAULT CURRENT_TIMESTAMP",
    }
    for column, definition in workflow_run_additions.items():
        if column not in workflow_run_columns:
            db.execute(f"ALTER TABLE workflow_runs ADD COLUMN {column} {definition}")

    # 图形美化工作台：原始图片、编辑版本与分析任务解耦，确保编辑永不覆盖
    # 原始分析输出。
    db.executescript("""
        CREATE TABLE IF NOT EXISTS figure_assets (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            label TEXT NOT NULL,
            file_type TEXT NOT NULL,
            file_path TEXT NOT NULL,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS figure_versions (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            source_kind TEXT NOT NULL,
            source_id TEXT NOT NULL,
            edit_mode TEXT NOT NULL,
            label TEXT NOT NULL,
            style_json TEXT DEFAULT '{}',
            png_path TEXT DEFAULT '',
            svg_path TEXT DEFAULT '',
            created_at TEXT DEFAULT CURRENT_TIMESTAMP
        );
    """)
    db.execute("CREATE INDEX IF NOT EXISTS idx_figure_assets_project ON figure_assets(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_figure_versions_project ON figure_versions(project_id)")

    # === Agent 系统表（阶段 1）===
    db.executescript("""
        CREATE TABLE IF NOT EXISTS agent_sessions (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            status TEXT DEFAULT 'active',
            goal_summary TEXT DEFAULT '',
            created_at TEXT DEFAULT CURRENT_TIMESTAMP,
            updated_at TEXT DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS agent_goals (
            id TEXT PRIMARY KEY,
            session_id TEXT NOT NULL REFERENCES agent_sessions(id) ON DELETE CASCADE,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            goal_type TEXT NOT NULL,
            target_entity TEXT DEFAULT '',
            goal_json TEXT NOT NULL,
            constraints_json TEXT DEFAULT '{}',
            status TEXT DEFAULT 'draft',
            created_at TEXT DEFAULT CURRENT_TIMESTAMP,
            updated_at TEXT DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS agent_steps (
            id TEXT PRIMARY KEY,
            session_id TEXT NOT NULL REFERENCES agent_sessions(id) ON DELETE CASCADE,
            goal_id TEXT REFERENCES agent_goals(id) ON DELETE SET NULL,
            step_index INTEGER DEFAULT 0,
            step_type TEXT DEFAULT '',
            thought_summary TEXT DEFAULT '',
            action_name TEXT DEFAULT '',
            action_args_json TEXT DEFAULT '{}',
            result_json TEXT DEFAULT '{}',
            status TEXT DEFAULT 'completed',
            created_at TEXT DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS analysis_branches (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            session_id TEXT REFERENCES agent_sessions(id) ON DELETE SET NULL,
            goal_id TEXT REFERENCES agent_goals(id) ON DELETE SET NULL,
            parent_task_id TEXT REFERENCES analysis_tasks(id) ON DELETE SET NULL,
            parent_adata_path TEXT NOT NULL,
            branch_name TEXT DEFAULT '',
            purpose TEXT DEFAULT '',
            params_json TEXT DEFAULT '{}',
            output_adata_path TEXT,
            status TEXT DEFAULT 'pending',
            accepted INTEGER DEFAULT 0,
            deleted INTEGER DEFAULT 0,
            error_traceback TEXT,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP,
            finished_at TEXT
        );

        CREATE TABLE IF NOT EXISTS candidate_scores (
            id TEXT PRIMARY KEY,
            branch_id TEXT NOT NULL REFERENCES analysis_branches(id) ON DELETE CASCADE,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            evaluator_name TEXT NOT NULL,
            target_label TEXT DEFAULT '',
            score_json TEXT NOT NULL,
            total_score REAL,
            recommendation TEXT DEFAULT '',
            created_at TEXT DEFAULT CURRENT_TIMESTAMP
        );
    """)
    db.execute("CREATE INDEX IF NOT EXISTS idx_agent_sessions_project ON agent_sessions(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_agent_goals_session ON agent_goals(session_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_agent_goals_project ON agent_goals(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_agent_steps_session ON agent_steps(session_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_agent_steps_goal ON agent_steps(goal_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_analysis_branches_project ON analysis_branches(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_analysis_branches_session ON analysis_branches(session_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_candidate_scores_branch ON candidate_scores(branch_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_candidate_scores_project ON candidate_scores(project_id)")

    # 迁移：为已有 analysis_tasks 表添加 branch_id 列
    try:
        db.execute("ALTER TABLE analysis_tasks ADD COLUMN branch_id TEXT REFERENCES analysis_branches(id) ON DELETE SET NULL")
    except Exception:
        pass  # 列已存在
    # 迁移：为已有 analysis_branches 表添加 deleted 列（软删除）
    try:
        db.execute("ALTER TABLE analysis_branches ADD COLUMN deleted INTEGER DEFAULT 0")
    except Exception:
        pass  # 列已存在

    db.execute("CREATE INDEX IF NOT EXISTS idx_tasks_branch ON analysis_tasks(branch_id)")

    # === Agent Jobs 表（Phase D：异步任务） ===
    db.executescript("""
        CREATE TABLE IF NOT EXISTS agent_jobs (
            id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
            session_id TEXT REFERENCES agent_sessions(id) ON DELETE SET NULL,
            goal_id TEXT REFERENCES agent_goals(id) ON DELETE SET NULL,
            job_type TEXT NOT NULL,
            status TEXT DEFAULT 'pending',
            progress INTEGER DEFAULT 0,
            current_step TEXT DEFAULT '',
            params_json TEXT DEFAULT '{}',
            result_json TEXT DEFAULT '{}',
            error_traceback TEXT,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP,
            started_at TEXT,
            finished_at TEXT
        );
    """)
    db.execute("CREATE INDEX IF NOT EXISTS idx_agent_jobs_project ON agent_jobs(project_id)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_agent_jobs_session ON agent_jobs(session_id)")

    db.commit()
    db.close()
