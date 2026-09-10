"""Small project-scoped WES dashboard and safe artifact downloads."""

import hashlib
import os

from flask import Blueprint, abort, flash, jsonify, redirect, render_template, request, send_file, url_for

from config import Config
from models import Project
from modules.workflows.artifacts import (
    get_workflow_artifact,
    is_visual_artifact,
    list_workflow_artifacts,
)
from modules.workflows.capture_uploads import store_capture_bed_upload
from modules.workflows.references import list_capture_kit_profiles, list_reference_assets
from modules.workflows.registry import list_workflows
from modules.workflows.runs import get_workflow_run, list_workflow_runs
from modules.workflows.sra import list_sra_jobs


wes_bp = Blueprint("wes", __name__)


_GERMLINE_REFERENCE_REQUIREMENTS = {
    "fasta", "fai", "dict", "dbsnp", "dbsnp_tbi",
    "germline_resource", "germline_resource_tbi",
}


def _reference_bundle_choices(references):
    """Build safe, exact-ID choices for the WES form without reading files."""
    grouped = {}
    for asset in references:
        bundle_id = str(asset.get("bundle_version") or "")
        if not bundle_id or bundle_id.startswith("capture:"):
            continue
        record = grouped.setdefault(bundle_id, {
            "id": bundle_id,
            "assemblies": set(),
            "validated_types": set(),
        })
        record["assemblies"].add(str(asset.get("assembly") or ""))
        if asset.get("status") == "validated":
            record["validated_types"].add(str(asset.get("asset_type") or ""))
    choices = []
    for record in grouped.values():
        record["assembly"] = next(iter(record["assemblies"]), "") if len(record["assemblies"]) == 1 else "mixed"
        record["production_allowed"] = _GERMLINE_REFERENCE_REQUIREMENTS.issubset(
            record["validated_types"]
        )
        choices.append({
            "id": record["id"],
            "assembly": record["assembly"],
            "production_allowed": record["production_allowed"],
        })
    return sorted(choices, key=lambda item: (not item["production_allowed"], item["assembly"], item["id"]))


@wes_bp.route("/<pid>/wes")
def index(pid):
    project = Project.get_by_id(pid)
    if not project:
        flash("项目未找到", "danger")
        return redirect(url_for("main.index"))
    runs = list_workflow_runs(pid)
    # Poll within the serving process so its in-memory subprocess handles are
    # available.  This lets a page refresh turn a finished Nextflow process
    # into its persisted terminal state without introducing a second worker.
    if any(run.get("status") in {"running", "cancel_requested"} for run in runs):
        from modules.workflows.nextflow import NextflowExecutor
        for run in runs:
            if run.get("status") in {"running", "cancel_requested"}:
                NextflowExecutor.poll(run["id"], pid)
        runs = list_workflow_runs(pid)
    from modules.workflows.stages import next_stage, stage_from_launch
    for run in runs:
        stage = stage_from_launch(run.get("launch") or {})
        run["stage"] = stage
        run["next_stage"] = (
            next_stage(stage["key"], stage["analysis_workflow_key"])
            if stage else None
        )
        run["artifacts"] = list_workflow_artifacts(run["id"], pid)
    has_active_runs = any(
        run.get("status") in {"running", "cancel_requested"} for run in runs
    )
    settings = {
        "executor_enabled": Config.WES_EXECUTOR_ENABLED,
        "profile": Config.WES_NEXTFLOW_PROFILE,
        "genome": Config.WES_NEXTFLOW_GENOME,
        "igenomes_configured": bool(Config.WES_NEXTFLOW_IGENOMES_BASE),
        "vep_configured": bool(Config.WES_NEXTFLOW_VEP_CACHE),
        "pon_configured": bool(Config.WES_NEXTFLOW_PON),
        "germline_resource_configured": bool(Config.WES_NEXTFLOW_GERMLINE_RESOURCE),
    }
    references = list_reference_assets()
    from modules.workflows.guided import WES_GUIDED_STAGES
    return render_template(
        "wes.html", project=project, workflows=list_workflows("wes"), runs=runs,
        capture_kits=list_capture_kit_profiles(),
        references=references, reference_bundles=_reference_bundle_choices(references),
        guided_stages=WES_GUIDED_STAGES, sra_jobs=list_sra_jobs(pid), settings=settings,
        has_active_runs=has_active_runs,
    )


@wes_bp.route("/<pid>/wes/runs/<run_id>")
def run_detail(pid, run_id):
    """Show one stage's own inputs, status and result artifacts."""
    project = Project.get_by_id(pid)
    if not project:
        abort(404)
    run = get_workflow_run(run_id, pid)
    if not run:
        abort(404)
    if run.get("status") in {"running", "cancel_requested"}:
        from modules.workflows.nextflow import NextflowExecutor
        run = NextflowExecutor.poll(run_id, pid)
    from modules.workflows.stages import next_stage, stage_from_launch
    stage = stage_from_launch(run.get("launch") or {})
    artifacts = list_workflow_artifacts(run_id, pid)
    # All reports and visual outputs remain project- and run-scoped through
    # their registered artifact IDs; filesystem paths never reach the page.
    visual_artifacts = [item for item in artifacts if is_visual_artifact(item)]
    return render_template(
        "wes_run.html", project=project, run=run, stage=stage,
        next_stage=(next_stage(stage["key"], stage["analysis_workflow_key"])
                    if stage else None),
        artifacts=artifacts, visual_artifacts=visual_artifacts,
        is_active=run.get("status") in {"running", "cancel_requested"},
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
    inline = item.get("artifact_kind") == "pipeline_provenance" or is_visual_artifact(item)
    response = send_file(real_path, as_attachment=not inline)
    if is_visual_artifact(item):
        # Generated reports and vector figures are useful to render in the
        # platform, but must not inherit the platform's origin or session when
        # opened directly in another tab.
        response.headers["Content-Security-Policy"] = "sandbox allow-scripts"
        response.headers["X-Content-Type-Options"] = "nosniff"
    return response
