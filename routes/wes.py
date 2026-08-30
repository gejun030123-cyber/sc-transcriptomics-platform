"""Small project-scoped WES dashboard and safe artifact downloads."""

import hashlib
import os

from flask import Blueprint, abort, flash, jsonify, redirect, render_template, request, send_file, url_for

from config import Config
from models import Project
from modules.workflows.artifacts import get_workflow_artifact, list_workflow_artifacts
from modules.workflows.capture_uploads import store_capture_bed_upload
from modules.workflows.references import list_capture_kit_profiles, list_reference_assets
from modules.workflows.registry import list_workflows
from modules.workflows.runs import list_workflow_runs
from modules.workflows.sra import list_sra_jobs


wes_bp = Blueprint("wes", __name__)


@wes_bp.route("/<pid>/wes")
def index(pid):
    project = Project.get_by_id(pid)
    if not project:
        flash("项目未找到", "danger")
        return redirect(url_for("main.index"))
    runs = list_workflow_runs(pid)
    for run in runs:
        run["artifacts"] = list_workflow_artifacts(run["id"], pid)
    settings = {
        "executor_enabled": Config.WES_EXECUTOR_ENABLED,
        "profile": Config.WES_NEXTFLOW_PROFILE,
        "genome": Config.WES_NEXTFLOW_GENOME,
        "igenomes_configured": bool(Config.WES_NEXTFLOW_IGENOMES_BASE),
        "vep_configured": bool(Config.WES_NEXTFLOW_VEP_CACHE),
        "pon_configured": bool(Config.WES_NEXTFLOW_PON),
        "germline_resource_configured": bool(Config.WES_NEXTFLOW_GERMLINE_RESOURCE),
    }
    return render_template(
        "wes.html", project=project, workflows=list_workflows("wes"), runs=runs,
        capture_kits=list_capture_kit_profiles(),
        references=list_reference_assets(), sra_jobs=list_sra_jobs(pid), settings=settings,
    )


@wes_bp.route("/<pid>/wes/capture-kits/upload", methods=["POST"])
def upload_capture_kit(pid):
    """Upload one project-provided calling BED and register it as test_only."""
    project = Project.get_by_id(pid)
    if not project:
        abort(404)
    uploaded = request.files.get("calling_bed") or request.files.get("bed") or request.files.get("file")
    try:
        result = store_capture_bed_upload(
            project_id=pid,
            file_storage=uploaded,
            capture_kit_id=request.form.get("capture_kit_id", ""),
            name=request.form.get("name", ""),
            version=request.form.get("version", ""),
            assembly=request.form.get("assembly", "GRCh38"),
        )
    except (OSError, ValueError) as exc:
        if request.accept_mimetypes.best == "application/json":
            return jsonify({"error": str(exc)}), 400
        flash(f"BED 上传失败：{exc}", "danger")
        return redirect(url_for("wes.index", pid=pid))
    profile = result["profile"]
    payload = {
        "message": "BED 已上传并登记为 test_only，需管理员审核后才能用于生产分析",
        "profile": profile,
        "checksum": result["checksum"],
        "size_bytes": result["size_bytes"],
    }
    if request.accept_mimetypes.best == "application/json":
        return jsonify(payload), 201
    flash(payload["message"], "success")
    return redirect(url_for("wes.index", pid=pid))


@wes_bp.route("/<pid>/wes/artifacts/<artifact_id>")
def artifact(pid, artifact_id):
    if not Project.get_by_id(pid):
        abort(404)
    item = get_workflow_artifact(artifact_id, pid)
    if not item:
        abort(404)
    path = os.path.abspath(str(item.get("file_path") or ""))
    run = next((candidate for candidate in list_workflow_runs(pid)
                if candidate["id"] == item["run_id"]), None)
    run_root = os.path.realpath(str((run or {}).get("run_dir") or ""))
    real_path = os.path.realpath(path)
    if (not run_root or not real_path.startswith(run_root + os.sep) or
            os.path.islink(path) or not os.path.isfile(real_path)):
        abort(403)
    digest = hashlib.sha256()
    with open(real_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != item.get("checksum"):
        abort(409, description="工件 checksum 已变化，拒绝下载")
    inline = item.get("artifact_kind") in {"multiqc_report", "pipeline_provenance"}
    return send_file(real_path, as_attachment=not inline)
