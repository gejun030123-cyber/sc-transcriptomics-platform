"""Constrained local Nextflow executor for the internal WES platform.

The web layer never accepts a command string. A :class:`LaunchSpec` is built
from a registry workflow, a server-rendered samplesheet and administrator
configuration; execution is disabled unless ``WES_EXECUTOR_ENABLED`` is
explicitly enabled.
"""

import json
import os
import signal
import shutil
import subprocess
import time
from datetime import datetime, timezone
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from config import Config

from .contracts import WorkflowSpec
from .runs import transition_workflow_run


class WorkflowNotConfiguredError(RuntimeError):
    pass


@dataclass(frozen=True)
class LaunchSpec:
    workflow: WorkflowSpec
    project_id: str
    run_id: str
    manifest_id: str
    command: List[str]
    work_dir: str
    results_dir: str
    samplesheet_path: str = ""
    profile: str = "docker"
    stdout_path: str = ""
    stderr_path: str = ""
    intervals_path: str = ""
    input_type: str = ""
    params_file_path: str = ""
    reference_bundle_id: str = ""
    capture_bed_id: str = ""
    parameters: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "workflow": self.workflow.key,
            "release": self.workflow.release,
            "project_id": self.project_id,
            "run_id": self.run_id,
            "manifest_id": self.manifest_id,
            "command": list(self.command),
            "work_dir": self.work_dir,
            "results_dir": self.results_dir,
            "samplesheet_path": self.samplesheet_path,
            "profile": self.profile,
            "stdout_path": self.stdout_path,
            "stderr_path": self.stderr_path,
            "intervals_path": self.intervals_path,
            "input_type": self.input_type,
            "params_file_path": self.params_file_path,
            "reference_bundle_id": self.reference_bundle_id,
            "capture_bed_id": self.capture_bed_id,
            "parameters": dict(self.parameters),
        }


# The platform is intentionally small and single-process. Keeping handles
# here lets the run status endpoint collect an exit code without a separate
# queue/worker service.
_PROCESS_HANDLES: Dict[str, subprocess.Popen] = {}
_PROCESS_LOG_HANDLES: Dict[str, tuple[Any, Any]] = {}


class NextflowExecutor:
    """Build and optionally execute a pinned nf-core/sarek command."""

    _WORKFLOW_STEPS = {
        "wes_germline": "mapping",
        "wes_somatic": "mapping",
        "wes_annotate_only": "annotate",
    }
    _WORKFLOW_TOOLS = {
        "wes_germline": "haplotypecaller",
        "wes_somatic": "mutect2",
        "wes_annotate_only": "vep",
    }

    def __init__(self, *, nextflow_bin: Optional[str] = None,
                 pipeline: Optional[str] = None, profile: Optional[str] = None,
                 enabled: Optional[bool] = None,
                 enforce_reference_catalog: Optional[bool] = None):
        self.nextflow_bin = nextflow_bin or Config.WES_NEXTFLOW_BIN
        self.pipeline = pipeline or Config.WES_NEXTFLOW_PIPELINE
        self.profile = profile or Config.WES_NEXTFLOW_PROFILE
        self.enabled = Config.WES_EXECUTOR_ENABLED if enabled is None else bool(enabled)
        self.enforce_reference_catalog = (
            Config.WES_REQUIRE_VALIDATED_REFERENCES
            if enforce_reference_catalog is None else bool(enforce_reference_catalog)
        )

    def prepare(self, workflow: WorkflowSpec, project_id: str, run_id: str,
                manifest_id: str, work_dir: str, results_dir: str,
                samplesheet_path: str = "", profile: Optional[str] = None,
                stdout_path: str = "", stderr_path: str = "",
                resume: bool = False, intervals_path: str = "",
                params_file_path: str = "", input_type: str = "",
                reference_bundle_id: str = "", capture_bed_id: str = "") -> LaunchSpec:
        selected_profile = profile or self.profile
        input_path = samplesheet_path or manifest_id
        normalized_input_type = str(input_type or "").strip().lower()
        step = self._WORKFLOW_STEPS.get(workflow.key)
        if workflow.key in {"wes_germline", "wes_somatic"} and normalized_input_type in {"bam", "cram"}:
            step = "variant_calling"
        tools = self._WORKFLOW_TOOLS.get(workflow.key)
        if not step or not tools:
            raise WorkflowNotConfiguredError(f"WES workflow 尚未定义 Sarek step/tools: {workflow.key}")
        if not Config.WES_NEXTFLOW_GENOME:
            raise WorkflowNotConfiguredError("尚未配置经验证的 WES_NEXTFLOW_GENOME")
        if workflow.key in {"wes_germline", "wes_somatic"} and not intervals_path:
            raise WorkflowNotConfiguredError("WES calling run 必须提供 capture BED")
        resolved_params_path = os.path.abspath(
            params_file_path or os.path.join(
                os.path.dirname(samplesheet_path) if samplesheet_path else os.path.dirname(work_dir),
                "parameters.json",
            )
        )
        parameters: Dict[str, Any] = {
            "input": os.path.abspath(input_path),
            "outdir": os.path.abspath(results_dir),
            "genome": Config.WES_NEXTFLOW_GENOME,
            "wes": True,
            "step": step,
            "tools": tools,
        }
        if intervals_path:
            parameters["intervals"] = os.path.abspath(intervals_path)
        if Config.WES_NEXTFLOW_IGENOMES_BASE:
            parameters["igenomes_base"] = Config.validate_wes_reference_path(
                Config.WES_NEXTFLOW_IGENOMES_BASE, require_directory=True
            )
        if workflow.key == "wes_annotate_only" and Config.WES_NEXTFLOW_VEP_CACHE:
            parameters["vep_cache"] = Config.validate_wes_reference_path(
                Config.WES_NEXTFLOW_VEP_CACHE, require_directory=True
            )
        if workflow.key == "wes_somatic":
            if Config.WES_NEXTFLOW_PON:
                parameters["pon"] = Config.validate_wes_reference_path(Config.WES_NEXTFLOW_PON)
            if Config.WES_NEXTFLOW_GERMLINE_RESOURCE:
                parameters["germline_resource"] = Config.validate_wes_reference_path(
                    Config.WES_NEXTFLOW_GERMLINE_RESOURCE
                )
        command = [
            self.nextflow_bin, "run", self.pipeline, "-r", workflow.release,
            "-profile", selected_profile,
        ]
        if resume:
            command.append("-resume")
        command.extend(["-params-file", resolved_params_path])
        return LaunchSpec(
            workflow=workflow,
            project_id=project_id,
            run_id=run_id,
            manifest_id=manifest_id,
            command=command,
            work_dir=os.path.abspath(work_dir),
            results_dir=os.path.abspath(results_dir),
            samplesheet_path=os.path.abspath(samplesheet_path) if samplesheet_path else "",
            profile=selected_profile,
            stdout_path=os.path.abspath(stdout_path) if stdout_path else "",
            stderr_path=os.path.abspath(stderr_path) if stderr_path else "",
            intervals_path=os.path.abspath(intervals_path) if intervals_path else "",
            input_type=normalized_input_type,
            params_file_path=resolved_params_path,
            reference_bundle_id=str(reference_bundle_id or "").strip(),
            capture_bed_id=str(capture_bed_id or "").strip(),
            parameters=parameters,
        )

    def _resolve_binary(self) -> str:
        candidate = str(self.nextflow_bin or "").strip()
        if not candidate:
            raise WorkflowNotConfiguredError("未配置 Nextflow 可执行文件")
        resolved = os.path.abspath(candidate) if os.path.sep in candidate else shutil.which(candidate)
        if not resolved or not os.path.isfile(resolved) or not os.access(resolved, os.X_OK):
            raise WorkflowNotConfiguredError(f"Nextflow 可执行文件不可用: {candidate}")
        return resolved

    def launch(self, spec: LaunchSpec, *, project_id: Optional[str] = None,
               expected_statuses: tuple[str, ...] = ("prepared",)):
        if not self.enabled:
            raise WorkflowNotConfiguredError(
                "WES_EXECUTOR_ENABLED 未开启；当前只允许生成并审阅 LaunchSpec"
            )
        if self.enforce_reference_catalog:
            from .references import launch_reference_readiness
            readiness = launch_reference_readiness(
                reference_bundle_id=spec.reference_bundle_id,
                capture_bed_id=spec.capture_bed_id,
                capture_bed_path=spec.intervals_path,
                workflow_key=spec.workflow.key,
            )
            if not readiness["valid"]:
                raise WorkflowNotConfiguredError(
                    "WES reference/capture 启动门禁未通过: " + "; ".join(readiness["errors"])
                )
            if not Config.WES_NEXTFLOW_IGENOMES_BASE:
                raise WorkflowNotConfiguredError("尚未配置本地 WES_NEXTFLOW_IGENOMES_BASE")
            if spec.workflow.key == "wes_annotate_only" and not Config.WES_NEXTFLOW_VEP_CACHE:
                raise WorkflowNotConfiguredError("注释 workflow 尚未配置本地 WES_NEXTFLOW_VEP_CACHE")
            if spec.workflow.key == "wes_somatic":
                if not Config.WES_NEXTFLOW_PON:
                    raise WorkflowNotConfiguredError("somatic workflow 尚未配置本地 WES_NEXTFLOW_PON")
                if not Config.WES_NEXTFLOW_GERMLINE_RESOURCE:
                    raise WorkflowNotConfiguredError(
                        "somatic workflow 尚未配置本地 WES_NEXTFLOW_GERMLINE_RESOURCE"
                    )
        binary = self._resolve_binary()
        if not spec.command or spec.command[0] != self.nextflow_bin:
            raise WorkflowNotConfiguredError("LaunchSpec 与当前 Nextflow 配置不一致")
        if not spec.samplesheet_path or not os.path.isfile(spec.samplesheet_path):
            raise WorkflowNotConfiguredError("LaunchSpec 缺少已生成的 samplesheet")
        if spec.workflow.key != "wes_annotate_only":
            if (not spec.intervals_path or not os.path.isfile(spec.intervals_path) or
                    os.path.islink(spec.intervals_path)):
                raise WorkflowNotConfiguredError("WES calling run 缺少有效且非符号链接的 capture BED")
        run_root = os.path.realpath(os.path.dirname(spec.work_dir))
        sample_path = os.path.realpath(spec.samplesheet_path)
        if (os.path.islink(spec.samplesheet_path) or
                not sample_path.startswith(os.path.join(run_root, "launch") + os.sep)):
            raise WorkflowNotConfiguredError("LaunchSpec samplesheet 不在当前 run 的 launch 目录内")
        params_path = os.path.realpath(spec.params_file_path)
        if (not spec.params_file_path or not os.path.isfile(params_path) or
                os.path.islink(spec.params_file_path) or
                not params_path.startswith(os.path.join(run_root, "launch") + os.sep)):
            raise WorkflowNotConfiguredError("LaunchSpec parameters.json 不在当前 run 的 launch 目录内")
        try:
            with open(params_path, encoding="utf-8") as handle:
                stored_parameters = json.load(handle)
        except (OSError, ValueError) as exc:
            raise WorkflowNotConfiguredError("LaunchSpec parameters.json 无法读取") from exc
        if stored_parameters != spec.parameters:
            raise WorkflowNotConfiguredError("LaunchSpec parameters.json 与固定运行参数不一致")
        for directory in (spec.work_dir, spec.results_dir):
            real_directory = os.path.realpath(directory)
            if not real_directory.startswith(run_root + os.sep):
                raise WorkflowNotConfiguredError("LaunchSpec 工作目录越过 run 目录边界")
        disk_parent = os.path.dirname(spec.results_dir)
        os.makedirs(disk_parent, exist_ok=True)
        free_gb = shutil.disk_usage(disk_parent).free / (1024 ** 3)
        if free_gb < max(0.0, float(Config.WES_MIN_FREE_GB)):
            raise WorkflowNotConfiguredError(
                f"WES 输出卷剩余空间 {free_gb:.1f} GB，低于启动门槛 {Config.WES_MIN_FREE_GB:.1f} GB"
            )
        os.makedirs(spec.work_dir, exist_ok=True)
        os.makedirs(spec.results_dir, exist_ok=True)
        stdout_path = spec.stdout_path or os.path.join(spec.work_dir, "..", "logs", "stdout.log")
        stderr_path = spec.stderr_path or os.path.join(spec.work_dir, "..", "logs", "stderr.log")
        stdout_path = os.path.abspath(stdout_path)
        stderr_path = os.path.abspath(stderr_path)
        os.makedirs(os.path.dirname(stdout_path), exist_ok=True)
        os.makedirs(os.path.dirname(stderr_path), exist_ok=True)
        out_handle = open(stdout_path, "ab")
        err_handle = open(stderr_path, "ab")
        try:
            command = list(spec.command)
            command[0] = binary
            process = subprocess.Popen(
                command, cwd=spec.work_dir, stdin=subprocess.DEVNULL,
                stdout=out_handle, stderr=err_handle, shell=False,
                start_new_session=True, close_fds=True,
            )
        except Exception:
            out_handle.close()
            err_handle.close()
            raise
        _PROCESS_HANDLES[spec.run_id] = process
        _PROCESS_LOG_HANDLES[spec.run_id] = (out_handle, err_handle)
        pgid = os.getpgid(process.pid) if os.name == "posix" else process.pid
        run = transition_workflow_run(
            spec.run_id, project_id or spec.project_id, "running",
            expected_statuses=expected_statuses, pid=process.pid,
            process_group_id=pgid, stdout_path=stdout_path,
            stderr_path=stderr_path,
            started_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        )
        return {"process": process, "pid": process.pid, "process_group_id": pgid,
                "stdout_path": stdout_path, "stderr_path": stderr_path, "run": run}

    @staticmethod
    def poll(run_id: str, project_id: str, *, collect_artifacts: bool = True):
        """Collect a local process exit and update the persisted lifecycle."""
        process = _PROCESS_HANDLES.get(run_id)
        if process is not None:
            code = process.poll()
            if code is None:
                from .runs import get_workflow_run
                current = get_workflow_run(run_id, project_id)
                if current and current.get("status") == "cancel_requested":
                    requested = current.get("cancel_requested_at") or ""
                    try:
                        age = (datetime.now(timezone.utc) -
                               datetime.fromisoformat(str(requested).replace("Z", "+00:00"))).total_seconds()
                    except (TypeError, ValueError):
                        age = 0
                    if age >= max(1, int(Config.WES_CANCEL_GRACE_SECONDS)):
                        try:
                            if os.name == "posix" and current.get("process_group_id"):
                                os.killpg(int(current["process_group_id"]), signal.SIGKILL)
                            else:
                                os.kill(int(process.pid), signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                return current or transition_workflow_run(
                    run_id, project_id, "running", expected_statuses=("running",)
                )
            _PROCESS_HANDLES.pop(run_id, None)
            handles = _PROCESS_LOG_HANDLES.pop(run_id, ())
            for handle in handles:
                try:
                    handle.close()
                except Exception:
                    pass
            from .runs import get_workflow_run
            current = get_workflow_run(run_id, project_id)
            target = "cancelled" if current and current.get("status") == "cancel_requested" else (
                "completed" if code == 0 else "failed"
            )
            completed_run = transition_workflow_run(
                run_id, project_id, target,
                expected_statuses=("running", "cancel_requested"), exit_code=code,
                finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            )
            if target == "completed" and collect_artifacts:
                from .artifacts import collect_workflow_artifacts
                try:
                    collection = collect_workflow_artifacts(run_id, project_id, strict=True)
                except (OSError, ValueError) as exc:
                    transition_workflow_run(
                        run_id, project_id, "failed", expected_statuses=("completed",),
                        error_text=f"WES 工件收集失败: {exc}",
                    )
                    collection = {"valid": False}
                if not collection["valid"]:
                    from .runs import get_workflow_run
                    return get_workflow_run(run_id, project_id)
            return completed_run
        from .runs import get_workflow_run
        run = get_workflow_run(run_id, project_id)
        if run and run.get("status") == "running":
            return transition_workflow_run(
                run_id, project_id, "interrupted", expected_statuses=("running",),
                error_text="平台进程重启后无法继续跟踪 Nextflow 子进程",
            )
        return run

    @staticmethod
    def cancel(run_id: str, project_id: str):
        from .runs import get_workflow_run
        run = get_workflow_run(run_id, project_id)
        if not run:
            return None
        if run.get("status") not in {"prepared", "running"}:
            return run
        pid = run.get("pid")
        pgid = run.get("process_group_id") or pid
        if not pid:
            return transition_workflow_run(
                run_id, project_id, "cancelled", expected_statuses=("prepared",),
                cancel_requested_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            )
        try:
            if os.name == "posix" and pgid:
                os.killpg(int(pgid), signal.SIGTERM)
            else:
                os.kill(int(pid), signal.SIGTERM)
        except ProcessLookupError:
            pass
        except PermissionError as exc:
            return transition_workflow_run(
                run_id, project_id, "failed", expected_statuses=("running",),
                error_text=f"取消 Nextflow 进程失败: {exc}",
            )
        return transition_workflow_run(
            run_id, project_id, "cancel_requested",
            expected_statuses=("running", "prepared"),
            cancel_requested_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        )
