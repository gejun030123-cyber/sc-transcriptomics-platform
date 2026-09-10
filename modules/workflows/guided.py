"""Constrained, reviewable WES stage configuration.

The browser may select from this small schema, but it never supplies a
command line, a tool name, or a reference path.  Values are normalized before
they enter a manifest and once more immediately before a LaunchSpec is made.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping


WES_GUIDED_STAGES = (
    {
        "key": "input",
        "title": "1. 输入与完整性检查",
        "description": "登记样本、参考 bundle 和 capture kit；可完整读取 FASTQ，并发现 gzip 截断文件。",
    },
    {
        "key": "preprocess",
        "title": "2. FASTQ 质控与预处理",
        "description": "选择 fastp 剪切、最短保留长度与可选固定端剪切。",
    },
    {
        "key": "mapping",
        "title": "3. 比对与覆盖度 QC",
        "description": "选择经验证的 BWA 比对器，以及是否保留对齐文件供人工复核。",
    },
    {
        "key": "calling",
        "title": "4. 变异检测与过滤",
        "description": "按胚系或肿瘤 workflow 使用固定 caller；可选择联合胚系、VCF 过滤与规范化。",
    },
    {
        "key": "review",
        "title": "5. 审阅、提交与结果",
        "description": "审阅不可变 LaunchSpec 后确认提交；运行状态、QC、VCF 与 MultiQC 在本页持续更新。",
    },
)


_BASE_DEFAULTS = {
    "trim_fastq": True,
    "trim_nextseq": False,
    "length_required": 15,
    "clip_r1": 0,
    "clip_r2": 0,
    "save_trimmed": False,
    "aligner": "bwa-mem",
    "save_mapped": False,
    "save_output_as_bam": False,
    "filter_vcfs": True,
    "normalize_vcfs": True,
}

_WORKFLOW_DEFAULTS = {
    "wes_germline": {"joint_germline": False},
    "wes_somatic": {},
    "wes_annotate_only": {},
}

# These are part of the registered workflow contract, not user-facing knobs.
# They are added to a LaunchSpec every time but cannot be overridden by a
# browser request or a hand-edited manifest.
_FIXED_BY_WORKFLOW = {
    "wes_somatic": {"only_paired_variant_calling": True},
}

_ALLOWED_BY_WORKFLOW = {
    "wes_germline": set(_BASE_DEFAULTS) | {"joint_germline"},
    "wes_somatic": set(_BASE_DEFAULTS),
    "wes_annotate_only": {"filter_vcfs", "normalize_vcfs"},
}

_BOOLEAN_KEYS = {
    "trim_fastq", "trim_nextseq", "save_trimmed", "save_mapped",
    "save_output_as_bam", "filter_vcfs", "normalize_vcfs",
    "joint_germline", "only_paired_variant_calling",
}
_INTEGER_LIMITS = {
    "length_required": (1, 300),
    "clip_r1": (0, 75),
    "clip_r2": (0, 75),
}
_ALIGNERS = {"bwa-mem", "bwa-mem2"}


def guided_defaults(workflow_key: str) -> Dict[str, Any]:
    """Return a fresh set of safe defaults for one registered WES workflow."""
    if workflow_key not in _WORKFLOW_DEFAULTS:
        raise ValueError(f"未知 WES workflow: {workflow_key}")
    defaults = dict(_BASE_DEFAULTS)
    defaults.update(_WORKFLOW_DEFAULTS[workflow_key])
    allowed = _ALLOWED_BY_WORKFLOW[workflow_key]
    result = {key: value for key, value in defaults.items() if key in allowed}
    result.update(_FIXED_BY_WORKFLOW.get(workflow_key, {}))
    return result


def _as_bool(value: Any, label: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in {0, 1}:
        return bool(value)
    raise ValueError(f"{label} 必须是布尔值")


def _as_int(value: Any, label: str, lower: int, upper: int) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} 必须是整数")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} 必须是整数") from exc
    if str(value).strip() not in {str(parsed), f"+{parsed}"} and not isinstance(value, int):
        raise ValueError(f"{label} 必须是整数")
    if not lower <= parsed <= upper:
        raise ValueError(f"{label} 必须介于 {lower} 和 {upper} 之间")
    return parsed


def normalize_guided_options(raw: Mapping[str, Any] | None, workflow_key: str,
                             *, allow_fixed: bool = False) -> Dict[str, Any]:
    """Validate the browser's limited option object and fill every default.

    This deliberately rejects unknown keys so a saved manifest cannot smuggle
    additional Sarek parameters, paths, shell fragments, or a caller change
    into a later run.
    """
    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raise ValueError("pipeline_options 必须是对象")
    defaults = guided_defaults(workflow_key)
    allowed = _ALLOWED_BY_WORKFLOW[workflow_key]
    fixed = _FIXED_BY_WORKFLOW.get(workflow_key, {})
    accepted = allowed | (set(fixed) if allow_fixed else set())
    unknown = sorted(set(raw) - accepted)
    if unknown:
        raise ValueError("pipeline_options 包含不支持的参数: " + ", ".join(unknown))
    normalized = dict(defaults)
    for key, value in raw.items():
        if key in fixed:
            if _as_bool(value, key) != fixed[key]:
                raise ValueError(f"{key} 是 {workflow_key} 的固定参数，不能修改")
            continue
        if key in _BOOLEAN_KEYS:
            normalized[key] = _as_bool(value, key)
        elif key in _INTEGER_LIMITS:
            normalized[key] = _as_int(value, key, *_INTEGER_LIMITS[key])
        elif key == "aligner":
            chosen = str(value or "").strip()
            if chosen not in _ALIGNERS:
                raise ValueError("aligner 只支持 bwa-mem 或 bwa-mem2")
            normalized[key] = chosen
    if not normalized.get("trim_fastq", False):
        if normalized.get("trim_nextseq") or normalized.get("save_trimmed") or any(
            normalized.get(key, 0) for key in ("clip_r1", "clip_r2")
        ):
            raise ValueError("关闭 fastp 时不能设置剪切、poly-G 去除或保存 trimmed FASTQ")
    if normalized.get("save_output_as_bam") and not normalized.get("save_mapped"):
        raise ValueError("save_output_as_bam 需要同时启用 save_mapped")
    return normalized
