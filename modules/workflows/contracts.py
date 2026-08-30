"""Stable contracts shared by workflow adapters and API layers."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


@dataclass(frozen=True)
class WorkflowSpec:
    key: str
    display_name: str
    description: str
    assay_type: str = "wes"
    modes: Tuple[str, ...] = ()
    input_types: Tuple[str, ...] = ()
    executor: str = "nextflow"
    release: str = "3.10.0"
    status: str = "planned"
    confirmation_required: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "key": self.key,
            "display_name": self.display_name,
            "description": self.description,
            "assay_type": self.assay_type,
            "modes": list(self.modes),
            "input_types": list(self.input_types),
            "executor": self.executor,
            "release": self.release,
            "status": self.status,
            "confirmation_required": self.confirmation_required,
        }


@dataclass
class PreflightResult:
    valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    normalized_manifest: Optional[Dict[str, Any]] = None
    checks: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "valid": self.valid,
            "errors": list(self.errors),
            "warnings": list(self.warnings),
            "manifest": self.normalized_manifest,
            "checks": dict(self.checks),
        }


@dataclass
class AnalysisResult:
    """Generic result envelope; external executors may add artifact metadata."""

    status: str
    primary_output: Optional[Dict[str, Any]] = None
    artifacts: List[Dict[str, Any]] = field(default_factory=list)
    metrics: Dict[str, Any] = field(default_factory=dict)
    summary: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    provenance: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status,
            "primary_output": self.primary_output,
            "artifacts": list(self.artifacts),
            "metrics": dict(self.metrics),
            "summary": dict(self.summary),
            "warnings": list(self.warnings),
            "provenance": dict(self.provenance),
        }
