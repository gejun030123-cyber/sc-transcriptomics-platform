"""Small, reusable experiment-design checks shown before an analysis starts.

The platform intentionally keeps metadata in ``AnnData.obs``.  This module
does not invent an experimental design or alter the input; it surfaces the
available columns, checks the selected design, and distinguishes hard input
requirements from analyses that are merely exploratory.
"""

from __future__ import annotations

from collections import Counter

from modules.io_utils import (
    infer_expression_measurement,
    infer_sample_group_candidates,
    obs_grouping_info,
)


_STATUS_RANK = {"pass": 0, "review": 1, "warning": 2, "blocked": 3}
_TECHNICAL_BATCH_NAMES = {
    "batch", "technical_batch", "sequencing_batch", "library_batch",
    "run", "lane", "sequencing_run",
}
_CONDITION_TOKENS = (
    "condition", "treatment", "group", "disease", "status", "cohort",
    "response", "phenotype", "case", "control", "stim", "drug",
)
_SAMPLE_TOKENS = ("sample", "library", "orig.ident", "specimen", "patient_sample")
_BATCH_TOKENS = ("batch", "lane", "sequencing", "library", "run")
_CELLTYPE_TOKENS = ("celltype", "cell_type", "annotation", "cluster", "leiden", "louvain")
_TIME_TOKENS = ("time", "day", "hour", "week", "stage", "minute")


def _check(name, status, message, value=None):
    item = {"name": str(name), "status": str(status), "message": str(message)}
    if value is not None:
        item["value"] = value
    return item


def _unique_strings(values, limit=80):
    result = []
    for value in values:
        text = str(value)
        if text not in result:
            result.append(text)
        if len(result) >= limit:
            break
    return result


def _column_profile(adata, column):
    values = adata.obs[column]
    non_missing = values.dropna()
    string_values = non_missing.astype(str).str.strip()
    missing = int(values.isna().sum() + string_values.eq("").sum())
    counts = Counter(string_values[string_values.ne("")].tolist())
    return {
        "column": str(column),
        "n_unique": int(values.nunique(dropna=True)),
        "n_missing": missing,
        "values": _unique_strings(counts.keys(), limit=12),
        "counts": {key: int(value) for key, value in counts.items()},
    }


def _candidate_columns(adata):
    """Return cautious role candidates without assuming their biological meaning."""
    columns = [str(column) for column in adata.obs.columns]
    grouping = []
    profiles = {}
    for column in columns:
        profile = _column_profile(adata, column)
        profiles[column] = profile
        info = obs_grouping_info(
            adata, column, max_categories=80, max_numeric_categories=30,
            require_multiple=True,
        )
        if info.get("valid"):
            grouping.append(column)

    lower = {column: column.lower() for column in columns}
    by_tokens = lambda tokens, source=grouping: [
        column for column in source
        if any(token in lower[column] for token in tokens)
    ]

    # A sample identifier may legitimately have more categories than a plotting
    # grouping.  Do not accept nearly one-cell-per-ID columns as a replicate key.
    sample_candidates = []
    for column in columns:
        n_unique = profiles[column]["n_unique"]
        if not any(token in lower[column] for token in _SAMPLE_TOKENS):
            continue
        if n_unique < 2 or n_unique > max(20, int(0.5 * max(1, adata.n_obs))):
            continue
        sample_candidates.append(column)

    from modules.sc_timecourse import _time_value
    import pandas as pd

    time_candidates = []
    for column in columns:
        profile = profiles[column]
        values = profile["values"]
        numeric = pd.api.types.is_numeric_dtype(adata.obs[column])
        looks_temporal = bool(values) and all(_time_value(value) is not None for value in values)
        if 3 <= profile["n_unique"] <= 30 and (
            numeric or looks_temporal or any(token in lower[column] for token in _TIME_TOKENS)
        ):
            time_candidates.append(column)

    def ordered(preferred, fallback):
        return _unique_strings([*preferred, *fallback], limit=80)

    condition_candidates = ordered(by_tokens(_CONDITION_TOKENS), grouping)
    batch_candidates = ordered(by_tokens(_BATCH_TOKENS), [])
    celltype_candidates = ordered(by_tokens(_CELLTYPE_TOKENS), [
        column for column in ("celltype", "annotation", "leiden") if column in columns
    ])
    sample_candidates = ordered([
        column for exact in ("sample_id", "sample", "library_id", "orig.ident")
        for column in columns if lower[column] == exact
    ], sample_candidates)

    auto_group_candidates = infer_sample_group_candidates(adata.obs.index.tolist())
    return {
        "columns": profiles,
        "grouping_candidates": grouping,
        "condition_candidates": condition_candidates,
        "sample_candidates": sample_candidates,
        "batch_candidates": batch_candidates,
        "celltype_candidates": celltype_candidates,
        "time_candidates": time_candidates,
        "auto_group_candidates": auto_group_candidates,
    }


def _selected_value(params, key, default=""):
    value = (params or {}).get(key, default)
    return str(value or "").strip()


def _as_int(value, default):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return int(default)


def _selected_groupby(module_name, params, candidates):
    if module_name == "proportion":
        return _selected_value(params, "condition_key") or _selected_value(params, "batch_key")
    if module_name in {"bulk_deg", "deg"}:
        selected = _selected_value(params, "groupby")
        if selected:
            return selected
        if module_name == "bulk_deg":
            # The actual Bulk DEG runner can safely materialize a name-based
            # mapping only when no obs grouping exists.  If obs already has a
            # candidate condition, require the caller to map it explicitly
            # instead of silently guessing which factor is the contrast.
            if not candidates.get("grouping_candidates") and candidates.get("auto_group_candidates"):
                return "_auto_group_"
        return ""
    if module_name == "sc_timecourse":
        return _selected_value(params, "condition_key")
    return _selected_value(params, "groupby") or _selected_value(params, "group_column")


def _group_counts(adata, column):
    if not column or column not in adata.obs.columns:
        return {}
    values = adata.obs[column]
    missing = values.isna() | values.astype(str).str.strip().eq("")
    return {
        str(key): int(value)
        for key, value in values.loc[~missing].astype(str).value_counts().items()
    }


def _make_contrast(group_counts, params):
    groups = sorted(group_counts)
    selected_1 = _selected_value(params, "group1")
    selected_2 = _selected_value(params, "group2")
    if selected_1 not in group_counts:
        selected_1 = groups[0] if groups else ""
    if selected_2 not in group_counts or selected_2 == selected_1:
        selected_2 = next((group for group in groups if group != selected_1), "")
    pairs = []
    for index, group_a in enumerate(groups):
        for group_b in groups[index + 1:]:
            count_a, count_b = group_counts[group_a], group_counts[group_b]
            pairs.append({
                "group1": group_a,
                "group2": group_b,
                "n_group1": int(count_a),
                "n_group2": int(count_b),
                "inference_ready": bool(count_a >= 2 and count_b >= 2),
            })
    return {
        "groups": [{"name": key, "n": int(group_counts[key])} for key in groups],
        "selected": {"group1": selected_1, "group2": selected_2},
        "pairs": pairs,
    }


def _sample_design_check(adata, sample_key, condition_key, min_samples=2):
    """Return a conservative sample-level design check for SC composition/time."""
    if not sample_key:
        return "warning", "未指定生物学样本列；仅可解释为描述性细胞组成。", None
    if sample_key not in adata.obs.columns:
        return "warning", f"样本列 '{sample_key}' 不存在；仅可解释为描述性细胞组成。", None
    profile = _column_profile(adata, sample_key)
    if profile["n_missing"]:
        return "warning", f"样本列 '{sample_key}' 有缺失值，已关闭样本级推断。", None
    per_sample_counts = list(profile["counts"].values())
    median_cells = sorted(per_sample_counts)[len(per_sample_counts) // 2] if per_sample_counts else 0
    if (profile["n_unique"] < 2
            or profile["n_unique"] > max(20, int(0.5 * max(1, adata.n_obs)))
            or median_cells < 2):
        return "warning", f"样本列 '{sample_key}' 接近一细胞一个 ID 或独立样本不足，不能作为统计重复。", None
    if sample_key.lower() in _TECHNICAL_BATCH_NAMES:
        return "warning", f"'{sample_key}' 看起来是技术 batch；未确认其为生物学样本，已关闭样本级推断。", None
    if not condition_key:
        return "review", f"已识别 {profile['n_unique']} 个样本；未指定条件列，仅输出样本级描述统计。", None
    if condition_key not in adata.obs.columns:
        return "warning", f"条件列 '{condition_key}' 不存在；无法构建样本级比较。", None

    import pandas as pd

    design = pd.DataFrame({
        "sample": adata.obs[sample_key].astype(str),
        "condition": adata.obs[condition_key].astype(str),
    })
    per_sample = design.groupby("sample", observed=True)["condition"].nunique()
    if bool((per_sample > 1).any()):
        return "warning", "同一个样本对应多个条件；当前小型检验不支持配对/交互模型，已关闭推断。", None
    condition_counts = design.drop_duplicates().groupby("condition", observed=True)["sample"].nunique()
    counts = {str(key): int(value) for key, value in condition_counts.items()}
    if len(counts) < 2:
        return "warning", "条件列只有一个有效取值，无法进行组间比较。", counts
    undersized = {key: value for key, value in counts.items() if value < int(min_samples)}
    if undersized:
        text = ", ".join(f"{key}=n{value}" for key, value in undersized.items())
        return "warning", f"以下条件的独立样本不足 {min_samples}：{text}；仅输出描述性比例。", counts
    return "pass", "样本与条件设计满足当前样本级比例比较的最小重复要求。", counts


def _missing_column_check(adata, role, column, role_candidates=None,
                          grouping_candidates=None):
    """Blocked design check for a required obs column that is absent."""
    columns = [str(item) for item in
               (grouping_candidates if grouping_candidates else adata.obs.columns)]
    shown = columns[:10]
    detail = "可用分组列：" + "、".join(shown) + ("…" if len(columns) > len(shown) else "") + "。"
    if role_candidates:
        detail += " 候选：" + "、".join(str(item) for item in role_candidates[:5]) + "。"
    detail += (" 请通过批量 10x manifest 导入写入 sample_id/condition，"
               "或在 h5ad.obs 中补充元数据后重试。")
    return _check(f"{role}列", "blocked", f"缺少{role}列 '{column}'；{detail}")


_CONDITION_RECOMMEND_NAMES = (
    "condition", "treatment", "treatment_group", "disease", "disease_status",
    "cohort", "response", "phenotype", "case", "control", "stim", "drug",
    "sample_group", "group",
)
_CONDITION_RECOMMEND_PREFIXES = (
    "condition", "treatment", "disease", "cohort", "response", "phenotype",
)
_CONDITION_EXCLUDE_TOKENS = (
    "annotation", "score", "state", "qc", "flag", "doublet", "marker",
    "passing", "candidate", "ambient", "lineage", "detected", "missing",
    "broad", "second", "top", "prediction",
)


def _token_condition_candidates(candidates):
    """Only condition-named columns are safe automatic recommendations.

    The generic candidate list also contains QC/annotation flags whose names
    embed words like "status" or "response"; those must never be auto-filled
    as an experimental condition.
    """
    safe = []
    for column in (candidates.get("condition_candidates") or []):
        lower = str(column).lower()
        if lower in _CONDITION_RECOMMEND_NAMES:
            safe.append(str(column))
            continue
        if lower.startswith(_CONDITION_RECOMMEND_PREFIXES) and not any(
            token in lower for token in _CONDITION_EXCLUDE_TOKENS
        ):
            safe.append(str(column))
    return safe


def build_design_preflight(adata, module_name, params=None, input_path=""):
    """Build JSON-safe design and contrast checks for a proposed analysis.

    ``blocked`` is deliberately narrow: it is reserved for prerequisites the
    corresponding analysis cannot run without.  A missing replicate design is
    a warning/exploratory state, not a reason to hide useful descriptive plots.
    """
    params = dict(params or {})
    checks = []
    candidates = _candidate_columns(adata)
    measurement = infer_expression_measurement(adata, input_path)
    selected_groupby = _selected_groupby(module_name, params, candidates)
    recommended = {}

    # Map roles to fields that actually exist in each small module.
    if module_name == "bulk_deg":
        suggested = selected_groupby or (candidates["condition_candidates"] or [""])[0]
        if suggested:
            recommended["groupby"] = suggested
    elif module_name == "deg":
        if not selected_groupby and candidates["grouping_candidates"]:
            recommended["groupby"] = candidates["grouping_candidates"][0]
    elif module_name == "proportion":
        if candidates["celltype_candidates"]:
            recommended["groupby"] = candidates["celltype_candidates"][0]
        if candidates["condition_candidates"]:
            recommended["condition_key"] = candidates["condition_candidates"][0]
        if candidates["sample_candidates"]:
            recommended["sample_key"] = candidates["sample_candidates"][0]
    elif module_name == "batch_correct" and candidates["batch_candidates"]:
        recommended["batch_key"] = candidates["batch_candidates"][0]
    elif module_name == "sc_timecourse":
        if candidates["time_candidates"]:
            recommended["timepoint_key"] = candidates["time_candidates"][0]
        if candidates["sample_candidates"]:
            recommended["sample_key"] = candidates["sample_candidates"][0]
        if candidates["condition_candidates"]:
            recommended["condition_key"] = candidates["condition_candidates"][0]
        if candidates["celltype_candidates"]:
            recommended["celltype_key"] = candidates["celltype_candidates"][0]

    count_only_methods = {
        "bulk_deg": {"deseq2", "edger"},
        "bulk_normalize": {"deseq2", "tmm", "cpm", "vst", "rlog"},
    }
    method = _selected_value(params, "method")
    if method in count_only_methods.get(module_name, set()):
        if measurement == "raw_counts":
            checks.append(_check("表达量尺度", "pass", f"检测到原始计数，可使用 {method}。", measurement))
        else:
            checks.append(_check(
                "表达量尺度", "blocked",
                f"当前数据为 {measurement}，不能使用需要原始计数的 {method}。",
                measurement,
            ))
    elif module_name.startswith("bulk_"):
        checks.append(_check("表达量尺度", "pass", f"检测到 {measurement}。", measurement))

    contrast = {"groups": [], "selected": {"group1": "", "group2": ""}, "pairs": []}
    needs_group = module_name in {"bulk_deg", "deg", "proportion"}
    if module_name == "proportion":
        selected_groupby = _selected_value(params, "condition_key") or _selected_value(params, "batch_key")
        if not selected_groupby:
            selected_groupby = (candidates["condition_candidates"] or [""])[0]

    if needs_group:
        if selected_groupby == "_auto_group_":
            auto_candidates = candidates.get("auto_group_candidates") or []
            if not auto_candidates:
                checks.append(_check("比较分组", "blocked", "无法从样本名推断可靠的重复分组。"))
            else:
                counts = {
                    str(key): int(value)
                    for key, value in (auto_candidates[0].get("group_sizes") or {}).items()
                }
                contrast = _make_contrast(counts, params)
                if any(value < 2 for value in counts.values()):
                    detail = ", ".join(f"{key}=n{value}" for key, value in counts.items())
                    checks.append(_check("生物学重复", "blocked", f"自动分组中每组至少需要 2 个样本；当前 {detail}。"))
                else:
                    checks.append(_check("比较分组", "pass", "已从样本名推断可复核的重复分组。", counts))
        elif not selected_groupby:
            status = "blocked" if module_name == "bulk_deg" else "warning"
            checks.append(_check("比较分组", status, "未选择可用的分类分组列。"))
        elif selected_groupby not in adata.obs.columns:
            status = "blocked" if module_name == "bulk_deg" else "warning"
            checks.append(_check("比较分组", status, f"分组列 '{selected_groupby}' 不在 adata.obs 中。"))
        else:
            info = obs_grouping_info(
                adata, selected_groupby, max_categories=80,
                max_numeric_categories=30, require_multiple=True,
            )
            if not info.get("valid"):
                status = "blocked" if module_name == "bulk_deg" else "warning"
                checks.append(_check("比较分组", status, info.get("reason", "分组列不可用")))
            else:
                counts = _group_counts(adata, selected_groupby)
                contrast = _make_contrast(counts, params)
                recommended.setdefault("groupby", selected_groupby)
                if len(counts) < 2:
                    checks.append(_check("比较分组", "blocked", "有效分组少于 2 个，无法比较。"))
                elif module_name == "bulk_deg" and any(value < 2 for value in counts.values()):
                    detail = ", ".join(f"{key}=n{value}" for key, value in counts.items())
                    checks.append(_check("生物学重复", "blocked", f"Bulk DEG 每组至少需要 2 个样本；当前 {detail}。"))
                elif any(value < 2 for value in counts.values()):
                    checks.append(_check("分组规模", "warning", "存在仅 1 个观测单位的分组；结果只能作为探索性证据。", counts))
                else:
                    checks.append(_check("比较分组", "pass", "已识别可用分组及可比较的组别。", counts))

    if module_name == "proportion":
        sample_key = _selected_value(params, "sample_key") or (candidates["sample_candidates"] or [""])[0]
        condition_key = _selected_value(params, "condition_key") or (candidates["condition_candidates"] or [""])[0]
        status, message, values = _sample_design_check(adata, sample_key, condition_key)
        checks.append(_check("样本级组成设计", status, message, values))
        if sample_key:
            recommended.setdefault("sample_key", sample_key)
        if condition_key:
            recommended.setdefault("condition_key", condition_key)

    if module_name == "sc_timecourse":
        time_key = _selected_value(params, "timepoint_key") or (candidates["time_candidates"] or [""])[0]
        sample_key = _selected_value(params, "sample_key") or (candidates["sample_candidates"] or [""])[0]
        if not time_key or time_key not in adata.obs.columns:
            checks.append(_check("时间点", "blocked", "未找到真实采样时间列；时序模块无法运行。"))
        else:
            profile = _column_profile(adata, time_key)
            if profile["n_missing"]:
                checks.append(_check("时间点", "blocked", f"时间列 '{time_key}' 有缺失值。"))
            elif profile["n_unique"] < 3:
                checks.append(_check("时间点", "blocked", "时序分析至少需要 3 个时间点。"))
            else:
                checks.append(_check("时间点", "pass", f"已识别 {profile['n_unique']} 个时间点。", time_key))
        sample_status, sample_message, sample_values = _sample_design_check(
            adata, sample_key, _selected_value(params, "condition_key"),
            min_samples=_as_int(params.get("min_replicates_per_timepoint", 2), 2),
        )
        checks.append(_check("样本重复", sample_status, sample_message, sample_values))

    if module_name == "batch_correct":
        batch_key = _selected_value(params, "batch_key") or (candidates["batch_candidates"] or [""])[0]
        if not batch_key or batch_key not in adata.obs.columns:
            checks.append(_check("批次列", "blocked", "未找到可用 batch 列，无法执行批次整合。"))
        else:
            info = obs_grouping_info(adata, batch_key, max_categories=80, max_numeric_categories=30, require_multiple=True)
            checks.append(_check(
                "批次列", "pass" if info.get("valid") else "blocked",
                "已识别可用 batch 列。" if info.get("valid") else info.get("reason", "batch 列不可用"),
                batch_key,
            ))

    if module_name == "sc_cell_deg":
        comparison_type = _selected_value(params, "comparison_type") or "condition"
        sample_key = _selected_value(params, "sample_key")
        condition_key = _selected_value(params, "condition_key")
        cluster_key = _selected_value(params, "cluster_key")
        if comparison_type == "condition":
            if not condition_key:
                checks.append(_check("条件列", "blocked", "未指定条件列；条件间细胞级比较需要一个条件列。"))
            elif condition_key not in adata.obs.columns:
                checks.append(_missing_column_check(
                    adata, "条件", condition_key,
                    role_candidates=_token_condition_candidates(candidates),
                    grouping_candidates=candidates["grouping_candidates"],
                ))
                token_condition = _token_condition_candidates(candidates)
                if token_condition:
                    recommended["condition_key"] = token_condition[0]
            else:
                profile = _column_profile(adata, condition_key)
                if profile["n_missing"]:
                    checks.append(_check("条件列", "blocked", f"条件列 '{condition_key}' 存在空值。"))
                elif profile["n_unique"] < 2:
                    checks.append(_check("条件列", "blocked", f"条件列 '{condition_key}' 只有一个取值，无法比较。"))
                else:
                    checks.append(_check("条件列", "pass", f"已识别 {profile['n_unique']} 个条件组。", condition_key))
            if condition_key:
                recommended.setdefault("condition_key", condition_key)
        else:
            if not sample_key:
                checks.append(_check("样本列", "blocked", "未指定样本列；该比较模式需要真实样本 ID 列。"))
            elif sample_key not in adata.obs.columns:
                checks.append(_missing_column_check(
                    adata, "样本", sample_key,
                    role_candidates=candidates["sample_candidates"],
                    grouping_candidates=candidates["grouping_candidates"],
                ))
                if candidates["sample_candidates"]:
                    recommended["sample_key"] = candidates["sample_candidates"][0]
            else:
                profile = _column_profile(adata, sample_key)
                if profile["n_missing"]:
                    checks.append(_check("样本列", "blocked", f"样本列 '{sample_key}' 存在空值。"))
                elif profile["n_unique"] < 2:
                    checks.append(_check("样本列", "blocked", f"样本列 '{sample_key}' 只有一个样本，无法比较。"))
                else:
                    checks.append(_check("样本列", "pass", f"已识别 {profile['n_unique']} 个样本。", sample_key))
                    if str(sample_key).lower() in _TECHNICAL_BATCH_NAMES:
                        checks.append(_check(
                            "样本列", "warning",
                            f"'{sample_key}' 看起来是技术 batch；未确认其为生物学样本，请人工核对实验设计后使用。",
                            sample_key,
                        ))
            if sample_key:
                recommended.setdefault("sample_key", sample_key)
            if comparison_type in {"within_sample_clusters", "between_samples_within_cluster"}:
                if not cluster_key:
                    checks.append(_check("聚类列", "blocked", "未指定聚类列。"))
                elif cluster_key not in adata.obs.columns:
                    checks.append(_missing_column_check(
                        adata, "聚类", cluster_key,
                        role_candidates=candidates["celltype_candidates"],
                        grouping_candidates=candidates["grouping_candidates"],
                    ))
                else:
                    checks.append(_check("聚类列", "pass", f"已使用聚类列 '{cluster_key}'。", cluster_key))
                if cluster_key:
                    recommended.setdefault("cluster_key", cluster_key)

    if module_name == "sc_pseudobulk_deg":
        sample_key = _selected_value(params, "sample_key")
        condition_key = _selected_value(params, "condition_key")
        cluster_key = _selected_value(params, "cluster_key")
        scope = _selected_value(params, "analysis_scope") or "per_cluster"
        token_condition = _token_condition_candidates(candidates)
        if not sample_key:
            checks.append(_check("样本列", "blocked", "未指定生物学样本列；pseudobulk DEG 需要真实样本 ID。"))
        elif sample_key not in adata.obs.columns:
            checks.append(_missing_column_check(
                adata, "样本", sample_key,
                role_candidates=candidates["sample_candidates"],
                grouping_candidates=candidates["grouping_candidates"],
            ))
            if candidates["sample_candidates"]:
                recommended["sample_key"] = candidates["sample_candidates"][0]
        else:
            profile = _column_profile(adata, sample_key)
            if profile["n_missing"]:
                checks.append(_check("样本列", "blocked", f"样本列 '{sample_key}' 存在空值。"))
            elif profile["n_unique"] < 2:
                checks.append(_check("样本列", "blocked", f"样本列 '{sample_key}' 只有一个样本，无法进行组间比较。"))
            else:
                checks.append(_check("样本列", "pass", f"已识别 {profile['n_unique']} 个样本。", sample_key))
                if str(sample_key).lower() in _TECHNICAL_BATCH_NAMES:
                    checks.append(_check(
                        "样本列", "warning",
                        f"'{sample_key}' 看起来是技术 batch；未确认其为生物学样本，pseudobulk 结论需人工核对实验设计。",
                        sample_key,
                    ))
        if sample_key:
            recommended.setdefault("sample_key", sample_key)

        if not condition_key:
            checks.append(_check("条件列", "blocked", "未指定条件列；pseudobulk DEG 需要每个样本一个 condition。"))
        elif condition_key not in adata.obs.columns:
            checks.append(_missing_column_check(
                adata, "条件", condition_key,
                role_candidates=token_condition,
                grouping_candidates=candidates["grouping_candidates"],
            ))
            if token_condition:
                recommended["condition_key"] = token_condition[0]
        else:
            profile = _column_profile(adata, condition_key)
            if profile["n_missing"]:
                checks.append(_check("条件列", "blocked", f"条件列 '{condition_key}' 存在空值。"))
            elif profile["n_unique"] < 2:
                checks.append(_check("条件列", "blocked", f"条件列 '{condition_key}' 只有一个取值。"))
            else:
                checks.append(_check("条件列", "pass", f"已识别 {profile['n_unique']} 个条件组。", condition_key))
        if condition_key:
            recommended.setdefault("condition_key", condition_key)

        if scope == "per_cluster":
            if not cluster_key:
                checks.append(_check("聚类列", "blocked", "未指定聚类列。"))
            elif cluster_key not in adata.obs.columns:
                checks.append(_missing_column_check(
                    adata, "聚类", cluster_key,
                    role_candidates=candidates["celltype_candidates"],
                    grouping_candidates=candidates["grouping_candidates"],
                ))
            else:
                checks.append(_check("聚类列", "pass", f"已使用聚类列 '{cluster_key}'。", cluster_key))
            if cluster_key:
                recommended.setdefault("cluster_key", cluster_key)
        batch_key = _selected_value(params, "batch_key")
        if batch_key and batch_key not in adata.obs.columns:
            checks.append(_check("批次列", "blocked", f"批次列 '{batch_key}' 不在 adata.obs 中。"))

    # A condition perfectly nested in a technical batch should not automatically
    # stop exploratory work, but the user needs to see it before integration.
    # 注意：proportion 模块的 groupby 是细胞类型列，不是条件列；只有真正的
    # condition_key（或伪 bulk/时序模块的 groupby 语义）才能参与混杂判断。
    condition_key = _selected_value(params, "condition_key")
    if not condition_key and module_name != "proportion":
        condition_key = _selected_value(params, "groupby")
    batch_key = _selected_value(params, "batch_key")
    if (
        condition_key and condition_key in adata.obs.columns
        and batch_key and batch_key in adata.obs.columns
        and condition_key != batch_key
    ):
        import pandas as pd

        n_batches = adata.obs[batch_key].nunique()
        table = pd.crosstab(adata.obs[condition_key].astype(str), adata.obs[batch_key].astype(str))
        if n_batches >= 2 and not table.empty and bool((table.gt(0).sum(axis=1) == 1).all()) and table.shape[0] > 1:
            checks.append(_check(
                "batch 与条件混杂", "warning",
                "每个条件只出现在一个 batch；批次校正可能移除真实条件差异，不能据此做因果解释。",
            ))

    # 分析范围（scope_key/scope_values）一致性检查：适用于所有带范围参数的模块
    scope_key = _selected_value(params, "scope_key")
    if scope_key:
        if scope_key not in adata.obs.columns:
            checks.append(_check("分析范围", "blocked", f"范围列 '{scope_key}' 不在 adata.obs 中。"))
        else:
            scope_values_raw = _selected_value(params, "scope_values")
            scope_values = [v.strip() for v in scope_values_raw.replace("\n", ",").split(",") if v.strip()]
            if not scope_values:
                checks.append(_check("分析范围", "warning", f"已填写范围列 '{scope_key}' 但未选择取值；将分析全部细胞。"))
            else:
                present = [v for v in scope_values if v in set(adata.obs[scope_key].astype(str))]
                missing = [v for v in scope_values if v not in present]
                if not present:
                    checks.append(_check("分析范围", "blocked", f"范围取值都不存在于 '{scope_key}'：{scope_values}"))
                elif missing:
                    checks.append(_check("分析范围", "warning", f"范围取值不存在于 '{scope_key}'：{missing}"))

    status = max((check["status"] for check in checks), key=lambda value: _STATUS_RANK.get(value, 1), default="review")
    if status == "pass":
        overall = "ready"
    elif status == "blocked":
        overall = "blocked"
    else:
        overall = "exploratory"
    blockers = [check["message"] for check in checks if check["status"] == "blocked"]
    warnings = [check["message"] for check in checks if check["status"] in {"warning", "review"}]
    return {
        "status": overall,
        "module_name": module_name,
        "measurement_type": measurement,
        "checks": checks,
        "metadata_roles": candidates,
        "selected_design": {
            "groupby": selected_groupby,
            "sample_key": _selected_value(params, "sample_key"),
            "condition_key": _selected_value(params, "condition_key"),
            "batch_key": _selected_value(params, "batch_key"),
            "timepoint_key": _selected_value(params, "timepoint_key"),
            "celltype_key": _selected_value(params, "celltype_key"),
        },
        "contrast": contrast,
        "recommended_params": recommended,
        "blockers": blockers,
        "warnings": warnings,
    }


def preflight_blockers(adata, module_name, params=None, input_path=""):
    """Return hard blockers for submission paths shared by form and AI tools."""
    return build_design_preflight(adata, module_name, params, input_path).get("blockers", [])
