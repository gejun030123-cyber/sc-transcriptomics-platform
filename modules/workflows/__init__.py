"""Workflow contracts and preflight helpers for non-AnnData assays.

The existing ``MODULE_REGISTRY`` remains the compatibility registry for the
single-cell and Bulk RNA modules.  WES workflows are deliberately registered
separately because they operate on sample manifests and heterogeneous file
artifacts rather than a single ``.h5ad`` path.
"""

from .contracts import AnalysisResult, PreflightResult, WorkflowSpec
from .preflight import validate_manifest
from .registry import WORKFLOW_REGISTRY, WES_WORKFLOW_REGISTRY, get_workflow, list_workflows
from .references import list_reference_assets
from .runs import (
    create_workflow_run, get_workflow_run, list_workflow_runs,
    transition_workflow_run,
)
from .sarek import render_samplesheet, write_launch_bundle

__all__ = [
    "AnalysisResult",
    "PreflightResult",
    "WorkflowSpec",
    "WORKFLOW_REGISTRY",
    "WES_WORKFLOW_REGISTRY",
    "get_workflow",
    "list_workflows",
    "validate_manifest",
    "list_reference_assets",
    "create_workflow_run",
    "get_workflow_run",
    "list_workflow_runs",
    "transition_workflow_run",
    "render_samplesheet",
    "write_launch_bundle",
]
