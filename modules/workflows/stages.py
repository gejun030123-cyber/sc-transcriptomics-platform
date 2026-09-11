"""Small, explicit WES stage contracts built on nf-core/sarek outputs.

Sarek's ``--step`` is a starting point, not an end point.  These contracts
therefore use the pipeline's published CSV hand-off files between separate,
reviewable local runs instead of pretending that a single Sarek invocation is
several independently completed analyses.
"""

from __future__ import annotations

from typing import Any, Dict, Optional


_ANALYSIS_WORKFLOWS = {"wes_germline", "wes_somatic"}


_STAGES = {
    "preprocess_mapping": {
        "key": "preprocess_mapping",
        "title": "1. FASTQ 预处理、比对与基础 QC",
        "description": "运行 fastp、比对、重复标记和 BQSR，输出可复用的 recalibrated.csv 与 MultiQC。",
        "run_workflow_key": None,
        "sarek_step": "mapping",
        "input_kind": "manifest",
        "output_csv": "recalibrated.csv",
        "next_stage": "variant_calling",
        "parameter_keys": (),
    },
    "variant_calling": {
        "key": "variant_calling",
        "title": "2. 变异调用、过滤与 VCF QC",
        "description": "以已完成的 recalibrated.csv 为输入，固定使用对应的 HaplotypeCaller 或 Mutect2。",
        "run_workflow_key": None,
        "sarek_step": "variant_calling",
        "input_kind": "recalibrated_csv",
        "output_csv": "variantcalled.csv",
        "next_stage": "annotation",
        "parameter_keys": ("filter_vcfs", "normalize_vcfs", "joint_germline"),
    },
    "annotation": {
        "key": "annotation",
        "title": "3. VCF 标准化与注释",
        "description": "以已完成的 variantcalled.csv 为输入，在已认证本地注释资源上运行 VEP。",
        "run_workflow_key": "wes_annotate_only",
        "sarek_step": "annotate",
        "input_kind": "variantcalled_csv",
        "output_csv": "",
        "next_stage": "",
        "parameter_keys": (),
    },
}


def get_stage(stage_key: str, analysis_workflow_key: str) -> Optional[Dict[str, Any]]:
    """Return one stage only for a supported germline/somatic analysis."""
    if analysis_workflow_key not in _ANALYSIS_WORKFLOWS:
        return None
    stage = _STAGES.get(str(stage_key or "").strip())
    if not stage:
        return None
    result = dict(stage)
    result["analysis_workflow_key"] = analysis_workflow_key
    result["run_workflow_key"] = result["run_workflow_key"] or analysis_workflow_key
    if result["key"] == "variant_calling":
        result["tools"] = "haplotypecaller" if analysis_workflow_key == "wes_germline" else "mutect2"
    elif result["key"] == "annotation":
        result["tools"] = "vep"
    else:
        # Empty tools deliberately means no variant caller or annotator in the
        # preprocessing run. Sarek still performs the requested mapping path.
        result["tools"] = None
    return result


def first_stage(analysis_workflow_key: str) -> Optional[Dict[str, Any]]:
    return get_stage("preprocess_mapping", analysis_workflow_key)


def next_stage(stage_key: str, analysis_workflow_key: str) -> Optional[Dict[str, Any]]:
    stage = get_stage(stage_key, analysis_workflow_key)
    if not stage or not stage.get("next_stage"):
        return None
    return get_stage(stage["next_stage"], analysis_workflow_key)


def stage_from_launch(launch: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Recover stage metadata persisted in a run's immutable launch bundle."""
    if not isinstance(launch, dict):
        return None
    return get_stage(
        str(launch.get("stage_key") or ""),
        str(launch.get("analysis_workflow_key") or launch.get("workflow") or ""),
    )
