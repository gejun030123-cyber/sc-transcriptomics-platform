"""Cell-level exploratory and sample-level pseudobulk enrichment.

The historical module name is kept for compatibility with existing projects,
but the implementation now has two explicit evidence levels:

* ``cell_level``: exploratory marker/condition results; cells are not treated
  as biological replicates.
* ``pseudobulk``: sample-level DEG results; ORA uses the tested-gene universe
  from each cluster and GSEA uses the complete ranked DEG table.
"""

import os
import json
import re
import ast
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from config import Config
from modules.base import BaseAnalysis
from modules.enrichment_statistics import (
    add_gsea_leading_edge_metrics,
    geneset_manifest,
    load_human_genesets,
    prepare_deg_gene_table,
    prepare_gsea_ranking,
    run_gsea_prerank,
    run_ora_full,
)
from modules.sc_batch import _available_path
from modules.sc_batch_export import batch_csv_package_dirs


GO_RESULT_COLUMNS = [
    "comparison_id", "comparison", "experimental_group", "control_group",
    "analysis_level", "statistical_status", "deg_scope", "cluster", "direction", "gene_set", "Term", "Overlap",
    "P-value", "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes",
    "n_input_genes", "n_background_genes", "method", "status", "reason",
]


LOCAL_GO_ASPECTS = {
    "BP": "GO_Biological_Process_2023",
    "CC": "GO_Cellular_Component_2023",
    "MF": "GO_Molecular_Function_2023",
}

GO_FOCUS_OVERVIEW_DATABASES = {
    "GO_Biological_Process_2023": "GO_BP",
    "GO_Cellular_Component_2023": "GO_CC",
    "GO_Molecular_Function_2023": "GO_MF",
}

# These Human snapshots are tested independently.  The mapping only provides
# stable display section names for the KEGG/Reactome/WikiPathways triptych; it
# never pools their p values or FDR families.
PATHWAY_TRIPTYCH_DATABASES = {
    "KEGG_2021_Human": "KEGG",
    "Reactome_2022": "Reactome",
    "WikiPathway_2021_Human": "WikiPathways",
}


def _safe_name(value, fallback="gene_set"):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return text or fallback


def _safe_int(value, default=0):
    parsed = pd.to_numeric(value, errors="coerce")
    return int(parsed) if pd.notna(parsed) and np.isfinite(parsed) else int(default)


def _single_text_value(frame, columns, default=""):
    """Read one provenance value and reject mixed contracts in one DEG file."""
    for column in columns:
        if column not in frame.columns:
            continue
        values = {
            str(value).strip() for value in frame[column].dropna()
            if str(value).strip()
        }
        if len(values) > 1:
            raise ValueError(f"DEG 文件的 {column} 列混有多个值")
        if values:
            return next(iter(values))
    return str(default)


def _analysis_scope_label(frame, deg_scope, cluster):
    """Return a user-facing unit name without losing legacy ``cluster`` data.

    DEG tables keep the ``cluster`` column as a stable interchange field.
    New annotation-aware tasks add ``analysis_group_label`` so an annotated
    cell type is never mislabeled as an unsupervised cluster in enrichment
    figures, audits, or result cards.
    """
    if str(deg_scope) == "all_cells":
        return "全体细胞"
    label = _single_text_value(frame, ("analysis_group_label",), "cluster")
    return f"{label} {cluster}"


def _requested_gene_sets(params):
    """Resolve one or more GO libraries without silently widening an analysis."""
    requested = str(params.get("gene_sets", "") or "").strip()
    if requested:
        libraries = [item.strip() for item in re.split(r"[,;\n]+", requested) if item.strip()]
    else:
        aspect_mode = str(params.get("go_aspects", "single") or "single")
        if aspect_mode == "BP_CC_MF":
            libraries = list(LOCAL_GO_ASPECTS.values())
        else:
            libraries = [str(params.get("gene_set", LOCAL_GO_ASPECTS["BP"]) or LOCAL_GO_ASPECTS["BP"])]
    deduplicated = []
    for library in libraries:
        if (
            library in {'.', '..'}
            or not re.fullmatch(r'[A-Za-z0-9_.-]+', library)
        ):
            raise ValueError(
                f"基因集名称 '{library}' 无效；只允许字母、数字、点、下划线和连字符。"
            )
        if library not in deduplicated:
            deduplicated.append(library)
    return deduplicated


def _requested_clusters(value):
    """Parse an explicit cluster subset while preserving exact cluster labels."""
    raw_values = [value] if isinstance(value, str) else list(value or [])
    requested = []
    for raw_value in raw_values:
        for item in re.split(r"[,;\n]+", str(raw_value or "")):
            cluster = item.strip()
            if not cluster:
                continue
            if len(cluster) > 256:
                raise ValueError("cluster 名称过长，无法作为富集运行范围")
            if cluster not in requested:
                requested.append(cluster)
    if len(requested) > 200:
        raise ValueError("一次最多选择 200 个 cluster 运行富集")
    return requested


def _local_gene_set_path(library, directory):
    """Return a platform-local GMT/TXT file for a configured GO library."""
    if library in {'.', '..'} or not re.fullmatch(r'[A-Za-z0-9_.-]+', str(library)):
        raise ValueError('本地 GO 基因集名称不是安全的文件标识符。')
    base = Path(directory).resolve()
    for suffix in (".gmt", ".txt"):
        candidate = base / f"{library}{suffix}"
        if candidate.is_symlink():
            raise ValueError('本地 GO 基因集文件不能是符号链接。')
        resolved = candidate.resolve()
        if base not in resolved.parents:
            raise ValueError('本地 GO 基因集文件必须位于选定目录内。')
        if resolved.is_file():
            return resolved
    raise FileNotFoundError(
        f"本地 GO 基因集缺失: {library}；请在 {base} 放置 {library}.gmt 或 .txt"
    )


def _resolve_local_gene_set_dir(value, project_dir):
    """Permit only the managed global cache or files inside this project."""
    managed_dir = (Path(Config.DATA_DIR) / "go_gene_sets").resolve()
    project_root = Path(project_dir).resolve()
    supplied = str(value or "").strip()
    if supplied:
        unresolved = Path(supplied)
        if not unresolved.is_absolute():
            unresolved = project_root / unresolved
        if unresolved.is_symlink():
            raise ValueError("本地 GO 基因集目录不能是符号链接")
        requested = unresolved.resolve()
    else:
        requested = managed_dir
    if not any(requested == root or root in requested.parents for root in (managed_dir, project_root)):
        raise ValueError("本地 GO 基因集目录必须位于平台缓存或当前项目目录内")
    return str(requested)


def _read_local_gene_sets(path):
    """Read Enrichr-style GMT/TXT without any network request."""
    gene_sets = {}
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term = fields[0].strip()
            genes = [gene.strip() for gene in fields[2:] if gene.strip()]
            if term and genes:
                gene_sets[term] = genes
    if not gene_sets:
        raise ValueError(f"本地 GO 基因集文件为空或格式无效: {path}")
    return gene_sets


def _run_ora(gp, selected, library, organism, execution_mode, local_gene_set_dir,
             background=None):
    """Run ORA locally when requested; Enrichr remains an explicit opt-in."""
    if execution_mode == "local":
        gene_sets = _read_local_gene_sets(_local_gene_set_path(library, local_gene_set_dir))
        kwargs = {
            "gene_list": selected, "gene_sets": gene_sets,
            "outdir": None, "no_plot": True,
        }
        if background:
            # gseapy's local ORA supports an explicit tested-gene universe.  It
            # is essential for cell-level results too: using the whole genome
            # as an implicit background overstates enrichment.
            kwargs["background"] = background
        return gp.enrich(**kwargs).results
    return gp.enrichr(
        gene_list=selected, gene_sets=[library], organism=organism,
        outdir=None, no_plot=True,
    ).results


def _latest_deg_files(deg_dir, prefix):
    """Return an unambiguous legacy CSV for each (comparison, DEG scope) pair.

    Older exports did not carry ``deg_scope``.  They remain readable by
    inferring whole-dataset versus per-cluster scope from the cluster column.
    """
    candidates = sorted(Path(deg_dir).glob(f"{prefix}_deg_*.csv"))
    latest = {}
    for path in candidates:
        # The combined pseudobulk table contains multiple comparison IDs and
        # is not a single analysis unit.  Only consume the per-comparison CSVs
        # here, otherwise the first two rows would silently define the scope.
        if path.stem.endswith("_deg_all_comparisons"):
            continue
        try:
            frame = pd.read_csv(path, nrows=2)
        except (OSError, pd.errors.EmptyDataError):
            continue
        if frame.empty or "comparison_id" not in frame.columns:
            continue
        comparison_id = str(frame["comparison_id"].iloc[0])
        if "deg_scope" in frame.columns:
            deg_scope = str(frame["deg_scope"].iloc[0])
        elif "cluster" in frame.columns and str(frame["cluster"].iloc[0]) == "All":
            deg_scope = "all_cells"
        else:
            deg_scope = "per_cluster"
        key = (comparison_id, deg_scope)
        if key in latest:
            raise ValueError(
                f"检测到多个旧版 DEG 文件对应 {comparison_id} / {deg_scope}；"
                "请在富集页面明确选择一个 DEG 任务，不能按文件时间自动猜测来源。"
            )
        latest[key] = path
    return latest


def _selected_deg_files(value, project_dir):
    """Validate explicitly selected internal DEG sources from one task."""
    if isinstance(value, str):
        values = [value] if value.strip() else []
    else:
        values = list(value or [])
    project_root = Path(project_dir).resolve()
    selected = {}
    for raw_path in values:
        supplied = Path(str(raw_path))
        if supplied.is_symlink():
            raise ValueError("所选 DEG 内部结果不能是符号链接")
        path = supplied.resolve()
        if project_root not in path.parents or not path.is_file():
            raise ValueError("所选 DEG 内部结果不存在或不属于当前项目")
        try:
            frame = pd.read_csv(path, nrows=2)
        except (OSError, pd.errors.EmptyDataError) as exc:
            raise ValueError(f"无法读取所选 DEG 结果: {path.name}") from exc
        if frame.empty or "comparison_id" not in frame.columns:
            raise ValueError(f"所选 DEG 结果缺少 comparison_id: {path.name}")
        comparison_id = str(frame["comparison_id"].iloc[0])
        deg_scope = str(frame.get("deg_scope", pd.Series(["all_cells"])).iloc[0])
        key = (comparison_id, deg_scope)
        if key in selected:
            raise ValueError(f"所选 DEG 任务含重复比较单元: {comparison_id} / {deg_scope}")
        selected[key] = path
    return selected


def _parse_focus_terms(value):
    """解析 focus_terms（列表或逗号/分号/换行分隔文本），去重并限量。

    AI 参数校验层会把 list 值字符串化（如 "['term1', 'term2']"），
    因此对每个拆分项额外剥离方括号与引号，避免残留杂质。
    """
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        items = [str(item).strip() for item in value]
    else:
        raw_value = str(value).strip()
        # AI 参数校验层会把列表转为 Python 字面量。优先安全地还原它，
        # 以免 term 本身含逗号时被错误拆分。
        if raw_value.startswith("[") and raw_value.endswith("]"):
            try:
                parsed = ast.literal_eval(raw_value)
            except (SyntaxError, ValueError):
                parsed = None
            if isinstance(parsed, (list, tuple)):
                items = [str(item).strip() for item in parsed]
            else:
                items = [item.strip() for item in re.split(r"[,;\n]+", raw_value)]
        else:
            items = [item.strip() for item in re.split(r"[,;\n]+", raw_value)]
    terms = []
    for item in items:
        cleaned = item.strip("[]'\" ").strip()
        if cleaned and cleaned not in terms:
            terms.append(cleaned)
    if len(terms) > 60:
        raise ValueError("focus_terms 最多支持 60 个通路 term")
    for term in terms:
        if len(term) > 200:
            raise ValueError("focus_terms 中单个通路名称过长")
    return terms


def _focus_term_mask(frame, focus_terms):
    """对完整 term 名称或去掉数据库 ID 后的显示名做精确匹配。"""
    if frame is None or frame.empty or "Term" not in frame.columns:
        return pd.Series(False, index=frame.index if frame is not None else [])
    terms = frame["Term"].astype(str)
    normalized = terms.str.lower().str.strip()
    display = terms.str.replace(r"\s*\([^)]*\)\s*$", "", regex=True).str.lower().str.strip()
    targets = {term.lower().strip() for term in focus_terms}
    return normalized.isin(targets) | display.isin(targets)


def _emit_focus_outputs(analysis, result, *, focus_terms, show_plots, plot_top_n,
                        enrichment_max_genes, similarity_threshold, running_term_n,
                        package_dirs, output_prefix, library, deg_scope,
                        comparison_id, comparison, method, source_level,
                        padj_cutoff, gsea_results_by_cluster, audit,
                        result_files, figure_audits, plot_warnings):
    """主题聚焦输出：全量统计与全量表保持不变，仅额外导出主题子表与聚焦图。

    聚焦图复用 figure engine 的指定通路选择（``pathway_selection=selected``），
    图中只展示用户确认的 term；FDR 仍为全库校正结果，审计中显式记录这一点。
    """
    if not focus_terms:
        return
    completed = result.loc[result["status"].eq("completed")]
    matched = completed.loc[_focus_term_mask(completed, focus_terms)].copy()
    matched_terms = sorted(matched["Term"].astype(str).unique().tolist())
    focus_audit = {
        "requested_terms": list(focus_terms),
        "matched_terms": matched_terms,
        "n_matched_rows": int(len(matched)),
        "matching_rule": "term 全名（忽略大小写和首尾空白）或去除末尾数据库 ID 后的显示名精确匹配；不使用子串匹配。",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "statistics_note": "FDR 在全库检验族上计算；聚焦输出仅为主题筛选展示，不是独立的校正检验族。",
        "status": "completed" if matched_terms else "no_match",
    }
    audit["focus"] = focus_audit
    if matched.empty:
        plot_warnings.append(
            f"{comparison} / {library}: focus_terms 未命中任何已完成通路，仅记录审计。"
        )
        return
    focus_path = _available_path(
        package_dirs["go"],
        f"{output_prefix}_focus_{_safe_name(library)}_{_safe_name(deg_scope)}_"
        f"{_safe_name(comparison_id, 'comparison')}", ".csv",
    )
    matched.to_csv(focus_path, index=False)
    focus_label = (
        f"{'Cell-level' if source_level == 'cell_level' else 'Pseudobulk'} {method} "
        f"focus terms ({library}, {deg_scope}): {comparison}"
    )
    result_files.append({
        "file_path": focus_path, "file_type": "csv", "category": "table",
        "label": focus_label,
    })
    if not show_plots:
        return
    from modules.bulk_enrichment import _render_enrichment_variants

    significant = matched["Significant"]
    if significant.dtype != bool:
        significant = significant.astype(str).str.lower().isin({"true", "1", "yes"})
    plot_frame = matched.loc[significant]
    if plot_frame.empty:
        plot_warnings.append(
            f"{comparison} / {library}: 聚焦通路中没有 FDR<{padj_cutoff:g} 的结果，"
            "不生成聚焦图；主题子表已导出。"
        )
        return
    level_title = (
        "Cell-level exploratory GO" if source_level == "cell_level"
        else f"Pseudobulk {method}"
    )
    for (cluster, direction), unit_frame in plot_frame.groupby(
        ["cluster", "direction"], sort=False, observed=True
    ):
        if unit_frame.empty:
            continue
        scope_label = _analysis_scope_label(unit_frame, deg_scope, cluster)
        cluster_suffix = "" if deg_scope == "all_cells" else f" · {scope_label}"
        direction_suffix = "" if method == "GSEA" else f" · {direction}"
        focus_params = dict(analysis.params)
        focus_params.update({
            "top_n": min(30, max(1, min(int(plot_top_n), int(unit_frame["Term"].nunique())))),
            "enrichment_plot_suite": "核心图",
            "pathway_selection": "selected",
            "target_pathways": list(focus_terms),
            "enrichment_max_genes": enrichment_max_genes,
            "similarity_threshold": similarity_threshold,
            "running_term_n": running_term_n,
        })
        stem = (
            f"{output_prefix}_focus_{_safe_name(library)}_{_safe_name(deg_scope)}_"
            f"{_safe_name(comparison_id)}_{_safe_name(cluster, 'All')}_"
            f"{_safe_name(direction, 'All')}"
        )
        semantic_warnings = ()
        if source_level == "cell_level":
            semantic_warnings = (
                "该富集基于细胞级探索性 DEG；细胞不是独立生物学重复，"
                "正式条件比较应以样本级 pseudobulk 结果验证。",
            )
        try:
            exported, audits = _render_enrichment_variants(
                analysis, unit_frame, method=method,
                title=f"{level_title} · {comparison}{cluster_suffix}{direction_suffix} · 聚焦",
                base_output_key=stem,
                label=f"{focus_label}: {scope_label} {direction}",
                params=focus_params,
                pre_res=gsea_results_by_cluster.get(str(cluster)),
                semantic_warnings=semantic_warnings,
            )
            result_files.extend(exported)
            figure_audits.extend(audits)
        except Exception as exc:
            plot_warnings.append(
                f"{comparison} / {library} / cluster {cluster} / {direction}: "
                f"聚焦绘图失败（{exc}）"
            )


_GO_PRIORITY_POLICY = (
    "每个 GO 本体按展示名额选择；先按炎症相关 term、再按脂代谢相关 term、"
    "再按其他确认主题 term、最后以全库 FDR 排名补足。"
)

# These controlled vocabulary terms intentionally reuse the platform theme
# lexicon wording.  They classify a *display priority* only: exact user
# focus_terms stay auditable separately and neither list changes ORA/FDR.
_INFLAMMATION_PRIORITY_KEYWORDS = (
    "inflamm", "cytokine", "chemokine", "interleukin", "interferon",
    "tumor necrosis factor", "nf-kappa", "toll-like", "nod-like",
    "acute phase", "inflammasome", "neutrophil", "leukocyte", "macrophage",
    "monocyte", "granulocyte", "antimicrobial", "defense response",
)
_LIPID_PRIORITY_KEYWORDS = (
    "lipid", "fatty acid", "acyl", "sterol", "cholesterol", "triglyceride",
    "phospholipid", "sphingolipid", "ceramide", "lipoprotein", "lipase",
    "glycerolipid", "eicosanoid", "prostaglandin", "leukotriene", "ketone", "ppar",
)
_GO_PRIORITY_ALLOCATION_DEFAULTS = {
    "inflammation_priority": 4,
    "lipid_metabolism_priority": 2,
    "confirmed_theme_priority": 0,
    "top_pathway_fill": 2,
}
_GO_PRIORITY_ALLOCATION_PARAMS = {
    "inflammation_priority": "go_inflammation_slots",
    "lipid_metabolism_priority": "go_lipid_slots",
    "confirmed_theme_priority": "go_confirmed_theme_slots",
    "top_pathway_fill": "go_top_pathway_slots",
}


def _go_priority_allocation_from_params(params, *, top_n):
    """Return the display allocation after strict, non-statistical validation."""
    allocation_mode = str(
        params.get("go_priority_allocation_mode", "balanced") or "balanced"
    ).strip().lower()
    if allocation_mode not in {"balanced", "custom"}:
        raise ValueError("go_priority_allocation_mode 必须为 balanced 或 custom")
    if allocation_mode == "balanced":
        return allocation_mode, None

    allocation = {}
    for reason, param_key in _GO_PRIORITY_ALLOCATION_PARAMS.items():
        raw_value = params.get(param_key, _GO_PRIORITY_ALLOCATION_DEFAULTS[reason])
        try:
            numeric_value = float(raw_value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{param_key} 必须为整数") from exc
        if not np.isfinite(numeric_value) or not numeric_value.is_integer():
            raise ValueError(f"{param_key} 必须为整数")
        allocation[reason] = int(numeric_value)
    if any(value < 0 or value > 30 for value in allocation.values()):
        raise ValueError("GO 三分区各类展示名额必须为 0–30 的整数")
    if sum(allocation.values()) > int(top_n):
        raise ValueError("GO 三分区各类展示名额之和不能超过每张富集图展示 term 数")
    return allocation_mode, allocation


def _completed_significant_go_rows(frame):
    """Keep only full-library FDR-significant ORA rows for GO visualisation."""
    required_columns = {"status", "Significant", "gene_set", "Term", "method", "direction"}
    if frame is None or not required_columns.issubset(frame.columns):
        return pd.DataFrame()
    completed = frame.loc[frame["status"].eq("completed")].copy()
    significant = completed["Significant"]
    if significant.dtype != bool:
        significant = significant.astype(str).str.lower().isin({"true", "1", "yes"})
    selected = completed.loc[
        significant & completed["gene_set"].isin(GO_FOCUS_OVERVIEW_DATABASES)
    ].copy()
    selected["Database"] = selected["gene_set"].map(GO_FOCUS_OVERVIEW_DATABASES)
    selected["Method"] = selected["method"].astype(str).str.upper()
    selected["Direction"] = selected["direction"].astype(str)
    return selected


def _focus_significant_rows(frame, focus_terms):
    """Compatibility helper for exact focus-term table exports.

    The new GO triptych deliberately uses theme-priority plus Top-pathway
    fill.  This helper retains the original exact-match behaviour for callers
    that need only the focused subset and for stable downstream integrations.
    """
    selected = _completed_significant_go_rows(frame)
    if selected.empty or not focus_terms:
        return selected.iloc[0:0].copy()
    return selected.loc[_focus_term_mask(selected, focus_terms)].copy()


def _go_display_priority(frame, focus_terms):
    """Tag rows without changing their test family or FDR values."""
    if frame.empty:
        result = frame.copy()
        result["selection_reason"] = pd.Series(dtype=str)
        result["selection_priority"] = pd.Series(dtype=int)
        return result
    display = frame["Term"].astype(str).str.replace(
        r"\s*\([^)]*\)\s*$", "", regex=True,
    ).str.lower().str.strip()
    inflammation = display.map(
        lambda value: any(keyword in value for keyword in _INFLAMMATION_PRIORITY_KEYWORDS)
    )
    lipid = display.map(
        lambda value: any(keyword in value for keyword in _LIPID_PRIORITY_KEYWORDS)
    )
    confirmed = _focus_term_mask(frame, focus_terms)
    result = frame.copy()
    result["selection_reason"] = "top_pathway_fill"
    result["selection_priority"] = 3
    result.loc[confirmed, ["selection_reason", "selection_priority"]] = (
        "confirmed_theme_priority", 2,
    )
    result.loc[lipid, ["selection_reason", "selection_priority"]] = (
        "lipid_metabolism_priority", 1,
    )
    result.loc[inflammation, ["selection_reason", "selection_priority"]] = (
        "inflammation_priority", 0,
    )
    return result


def _overlap_ratio(value):
    try:
        left, right = str(value).split("/", 1)
        denominator = float(right)
        return float(left) / denominator if denominator > 0 else np.nan
    except (TypeError, ValueError):
        return np.nan


def _go_priority_budgets(subset, *, top_n, priority_allocation=None):
    """Reserve Top-pathway capacity while retaining each available theme.

    For the default cap of eight, inflammation and lipid candidates receive a
    4:2 split and two rows remain for non-theme FDR Top pathways.  Empty
    categories automatically return their allocation, so sparse ontologies
    still show all available significant rows.
    """
    total = max(1, int(top_n))
    if priority_allocation is not None:
        requested = {
            "inflammation_priority": int(priority_allocation.get("inflammation_priority", 0)),
            "lipid_metabolism_priority": int(priority_allocation.get("lipid_metabolism_priority", 0)),
            "confirmed_theme_priority": int(priority_allocation.get("confirmed_theme_priority", 0)),
            "top_pathway_fill": int(priority_allocation.get("top_pathway_fill", 0)),
        }
        if any(value < 0 for value in requested.values()):
            raise ValueError("GO 三分区展示名额不能为负数")
        if sum(requested.values()) > total:
            raise ValueError("GO 三分区各类展示名额之和不能超过每本体展示 term 数")
        # The four values are *reserved display capacities*, not a statistical
        # filter.  Unused capacity is returned later to FDR-ranked ordinary
        # pathways first, then to remaining theme hits, and the requested
        # allocation is preserved verbatim in the audit JSON.
        return requested, {
            "theme_priority_budget": int(
                requested["inflammation_priority"]
                + requested["lipid_metabolism_priority"]
                + requested["confirmed_theme_priority"]
            ),
            "top_pathway_fill_reserved": int(requested["top_pathway_fill"]),
            "allocation_mode": "custom",
            "requested_allocation": requested,
        }
    active = [
        ("inflammation_priority", 2),
        ("lipid_metabolism_priority", 1),
        ("confirmed_theme_priority", 1),
    ]
    active = [
        (reason, weight) for reason, weight in active
        if subset["selection_reason"].eq(reason).any()
    ]
    if not active:
        return {"top_pathway_fill": total}, {
            "theme_priority_budget": 0, "top_pathway_fill_reserved": total,
            "allocation_mode": "balanced",
        }
    top_reserve = max(1, int(np.ceil(total * .25))) if total > 1 else 0
    theme_budget = max(1, total - top_reserve)
    weights = sum(weight for _reason, weight in active)
    caps = {
        reason: max(1, int(np.floor(theme_budget * weight / weights)))
        for reason, weight in active
    }
    while sum(caps.values()) > theme_budget:
        # Remove surplus from the lowest-priority category that still has more
        # than one reserved row, preserving at least one visible hit per theme.
        for reason, _weight in reversed(active):
            if caps[reason] > 1:
                caps[reason] -= 1
                break
        else:
            break
    while sum(caps.values()) < theme_budget:
        # Distribute remainder in priority order, giving inflammation the
        # additional capacity before lipid/other confirmed themes.
        for reason, _weight in active:
            caps[reason] += 1
            if sum(caps.values()) >= theme_budget:
                break
    caps["top_pathway_fill"] = top_reserve
    return caps, {
        "theme_priority_budget": int(theme_budget),
        "top_pathway_fill_reserved": int(top_reserve),
        "allocation_mode": "balanced",
    }


def _select_go_priority_rows(frame, focus_terms, *, top_n, priority_allocation=None):
    """Select display rows per ontology with a deterministic, recorded order."""
    candidates = _go_display_priority(_completed_significant_go_rows(frame), focus_terms)
    selected_parts = []
    audit = {"selection_policy": _GO_PRIORITY_POLICY, "ontologies": {}}
    for library, database in GO_FOCUS_OVERVIEW_DATABASES.items():
        subset = candidates.loc[candidates.get("gene_set", pd.Series(dtype=str)).eq(library)].copy()
        ontology_audit = {
            "database": database,
            "n_fdr_significant_candidates": int(len(subset)),
            "n_inflammation_candidates": int(subset.get("selection_reason", pd.Series(dtype=str)).eq("inflammation_priority").sum()),
            "n_lipid_metabolism_candidates": int(subset.get("selection_reason", pd.Series(dtype=str)).eq("lipid_metabolism_priority").sum()),
            "n_confirmed_theme_candidates": int(subset.get("selection_reason", pd.Series(dtype=str)).eq("confirmed_theme_priority").sum()),
            "display_cap": int(top_n),
        }
        if subset.empty:
            ontology_audit["n_selected"] = 0
            audit["ontologies"][database] = ontology_audit
            continue
        subset["_fdr_sort"] = pd.to_numeric(
            subset.get("Adjusted P-value", pd.Series(np.nan, index=subset.index)), errors="coerce",
        )
        subset["_ratio_sort"] = subset.get(
            "Overlap", pd.Series(np.nan, index=subset.index),
        ).map(_overlap_ratio)
        subset["_count_sort"] = subset.get("Overlap", pd.Series(np.nan, index=subset.index)).map(
            lambda value: _safe_int(str(value).split("/", 1)[0], 0)
        )
        ranked = subset.sort_values(
            ["selection_priority", "_fdr_sort", "_ratio_sort", "_count_sort", "Term"],
            ascending=[True, True, False, False, True], na_position="last", kind="mergesort",
        )
        caps, budget_audit = _go_priority_budgets(
            ranked, top_n=int(top_n), priority_allocation=priority_allocation,
        )
        ontology_audit.update(budget_audit)
        ontology_audit["selection_caps"] = {key: int(value) for key, value in caps.items()}
        selected_indices = []
        for reason in (
            "inflammation_priority", "lipid_metabolism_priority", "confirmed_theme_priority",
            "top_pathway_fill",
        ):
            selected_indices.extend(
                ranked.loc[ranked["selection_reason"].eq(reason)].head(caps.get(reason, 0)).index.tolist()
            )
        # Fill unused slots from ordinary Top pathways first.  Only if that
        # list is exhausted do extra theme candidates return, avoiding a large
        # inflammatory list silently crowding out non-theme biology.
        remaining = max(0, int(top_n) - len(selected_indices))
        if remaining:
            ordinary = ranked.loc[
                ranked["selection_reason"].eq("top_pathway_fill")
                & ~ranked.index.isin(selected_indices)
            ]
            selected_indices.extend(ordinary.head(remaining).index.tolist())
        remaining = max(0, int(top_n) - len(selected_indices))
        if remaining:
            fallback = ranked.loc[~ranked.index.isin(selected_indices)]
            selected_indices.extend(fallback.head(remaining).index.tolist())
        ordered = ranked.loc[selected_indices].copy()
        ordered["selection_rank"] = np.arange(1, len(ordered) + 1, dtype=int)
        ordered["selection_policy"] = _GO_PRIORITY_POLICY
        ontology_audit["n_selected"] = int(len(ordered))
        ontology_audit["selected_terms"] = [
            {"term": str(row.Term), "reason": str(row.selection_reason), "rank": int(row.selection_rank)}
            for row in ordered.itertuples(index=False)
        ]
        audit["ontologies"][database] = ontology_audit
        selected_parts.append(ordered.drop(columns=["_fdr_sort", "_ratio_sort", "_count_sort"]))
    if not selected_parts:
        empty = candidates.copy()
        for column, dtype in (("selection_rank", int), ("selection_policy", str)):
            if column not in empty:
                empty[column] = pd.Series(dtype=dtype)
        return empty, audit
    return pd.concat(selected_parts, ignore_index=True, sort=False), audit


def _select_pathway_triptych_rows(frame, *, top_n):
    """Select FDR-ranked rows separately for each Human pathway library.

    Unlike the GO theme-priority display, this view has no ontology-specific
    biological keyword policy.  It shows the top already-significant pathway
    terms from each database, and records the exact rows and order in its
    source table/audit.  The routine is deliberately display-only.
    """
    required_columns = {"status", "Significant", "gene_set", "Term", "method", "direction"}
    audit = {
        "selection_policy": "Each Human pathway library contributes FDR-ranked significant terms only; no cross-database FDR pooling or re-testing.",
        "databases": {},
    }
    if frame is None or not required_columns.issubset(frame.columns):
        return pd.DataFrame(), audit

    completed = frame.loc[frame["status"].eq("completed")].copy()
    significant = completed["Significant"]
    if significant.dtype != bool:
        significant = significant.astype(str).str.lower().isin({"true", "1", "yes"})
    candidates = completed.loc[
        significant & completed["gene_set"].isin(PATHWAY_TRIPTYCH_DATABASES)
    ].copy()
    candidates["Database"] = candidates["gene_set"].map(PATHWAY_TRIPTYCH_DATABASES)

    selected_parts = []
    for library, database in PATHWAY_TRIPTYCH_DATABASES.items():
        subset = candidates.loc[candidates["gene_set"].eq(library)].copy()
        database_audit = {
            "library": library,
            "database": database,
            "n_fdr_significant_candidates": int(len(subset)),
            "display_cap": int(top_n),
        }
        if subset.empty:
            database_audit["n_selected"] = 0
            audit["databases"][database] = database_audit
            continue
        subset["_fdr_sort"] = pd.to_numeric(
            subset.get("Adjusted P-value", pd.Series(np.nan, index=subset.index)),
            errors="coerce",
        )
        subset["_ratio_sort"] = subset.get(
            "Overlap", pd.Series(np.nan, index=subset.index),
        ).map(_overlap_ratio)
        subset["_count_sort"] = subset.get(
            "Overlap", pd.Series(np.nan, index=subset.index),
        ).map(lambda value: _safe_int(str(value).split("/", 1)[0], 0))
        selected = subset.sort_values(
            ["_fdr_sort", "_ratio_sort", "_count_sort", "Term"],
            ascending=[True, False, False, True], na_position="last", kind="mergesort",
        ).head(int(top_n)).copy()
        selected["selection_rank"] = np.arange(1, len(selected) + 1, dtype=int)
        selected["selection_reason"] = "fdr_top_pathway"
        selected["selection_priority"] = 0
        selected["selection_policy"] = audit["selection_policy"]
        database_audit["n_selected"] = int(len(selected))
        database_audit["selected_terms"] = [
            {"term": str(row.Term), "rank": int(row.selection_rank)}
            for row in selected.itertuples(index=False)
        ]
        audit["databases"][database] = database_audit
        selected_parts.append(selected.drop(columns=["_fdr_sort", "_ratio_sort", "_count_sort"]))

    if not selected_parts:
        return candidates.iloc[0:0].copy(), audit
    return pd.concat(selected_parts, ignore_index=True, sort=False), audit


def _emit_go_focus_overviews(analysis, frames, *, focus_terms, focus_label,
                             show_plots, plot_top_n, priority_allocation,
                             allocation_mode, package_dirs, output_prefix,
                             result_files, figure_audits, plot_warnings):
    """Export reference-style, theme-prioritised GO triptych dot and bar plots.

    Enrichment itself has already completed separately for all GO terms in
    each ontology.  This function does no testing.  It deterministically
    ranks only FDR-significant rows for display and writes the exact source
    table plus a per-ontology decision audit for reproduction.
    """
    if not frames:
        return 0
    combined = pd.concat(frames, ignore_index=True, sort=False)
    observed_libraries = set(combined.get("gene_set", pd.Series(dtype=str)).astype(str))
    required_libraries = set(GO_FOCUS_OVERVIEW_DATABASES)
    missing_libraries = sorted(required_libraries - observed_libraries)
    if missing_libraries:
        plot_warnings.append(
            "未生成主题优先 GO 三分区图：本次运行未同时包含 "
            + "、".join(missing_libraries)
        )
        return 0

    from figure_engine import NatureFigureDirector, export_registered_figure
    import matplotlib.pyplot as plt

    director = NatureFigureDirector()
    top_n = min(12, max(1, int(plot_top_n)))
    source_rows = combined.loc[
        combined.get("gene_set", pd.Series("", index=combined.index)).isin(
            GO_FOCUS_OVERVIEW_DATABASES
        )
    ].copy()
    group_columns = ["comparison_id", "comparison", "deg_scope", "cluster", "direction", "method"]
    selection_frames = []
    selection_audit = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "focus_label": focus_label,
        "requested_focus_terms": list(focus_terms),
        "statistics_note": "ORA 在完整 GO 本体分别检验，Adjusted P-value 为各完整本体内的 FDR；本审计仅记录展示排序，不重新校正。",
        "selection_policy": _GO_PRIORITY_POLICY,
        "allocation_mode": allocation_mode,
        "requested_allocation": priority_allocation,
        "display_cap_per_ontology": int(top_n),
        "units": [],
    }
    groups = list(source_rows.groupby(group_columns, dropna=False, sort=True, observed=True))
    for group_key, source_group in groups:
        comparison_id, comparison, deg_scope, cluster, direction, method = group_key
        if str(method).upper() != "ORA":
            plot_warnings.append(
                f"{comparison} / {deg_scope} / {direction}: 参考式 GO 气泡图和 Count 柱状图仅适用于 ORA；GSEA 保留原有图型。"
            )
            continue
        selected, unit_audit = _select_go_priority_rows(
            source_group, focus_terms, top_n=top_n,
            priority_allocation=priority_allocation,
        )
        scope_label = _analysis_scope_label(source_group, deg_scope, cluster)
        direction_label = "上调" if str(direction) == "Up" else "下调" if str(direction) == "Down" else str(direction)
        unit_audit.update({
            "comparison_id": str(comparison_id), "comparison": str(comparison),
            "deg_scope": str(deg_scope), "cluster": str(cluster), "direction": str(direction),
        })
        selection_audit["units"].append(unit_audit)
        if not selected.empty:
            selection_frames.append(selected)
        if not show_plots:
            continue
        n_terms = int(len(selected))
        height_mm = min(225, max(132, 62 + 6 * max(n_terms, 6)))
        key_base = (
            f"{output_prefix}_GO_priority_{_safe_name(deg_scope)}_"
            f"{_safe_name(comparison_id)}_{_safe_name(cluster, 'All')}_"
            f"{_safe_name(direction, 'All')}"
        )
        label_base = f"{focus_label} GO 主题优先 · {scope_label} · {direction_label}"
        views = (
            ("go_priority_dotplot", "dotplot", "GO Enrichment Plot"),
            ("go_priority_barplot", "barplot", f"GO enrichment of {scope_label}"),
        )
        for plot_type, suffix, heading in views:
            figure = None
            try:
                spec = director.create_spec(
                    plot_type, width="double", height_mm=height_mm, top_n=top_n,
                    formats=("svg", "pdf", "png"), database_scope=("GO_BP", "GO_CC", "GO_MF"),
                    title=f"{heading}\n{comparison} · {scope_label} · {direction_label} · {focus_label}",
                )
                figure = director.render(spec, selected)
                exported, report = export_registered_figure(
                    figure, os.path.join(analysis.ensure_plots_dir(), f"{key_base}_{suffix}"), spec,
                    category=plot_type,
                    label=f"{label_base} {('气泡图' if suffix == 'dotplot' else '柱状图')}（炎症→脂代谢→Top；全库 FDR）",
                    qa_path=os.path.join(
                        package_dirs["go"], f"{key_base}_{suffix}_nature_readiness.json",
                    ),
                )
                result_files.extend(exported)
                figure_audits.append({
                    "status": "pass" if report.ready else "warning",
                    "nature_readiness_score": report.score,
                    "issues": [issue.message for issue in report.issues],
                    "scope": scope_label,
                    "direction": direction_label,
                    "plot_type": plot_type,
                })
                if not report.ready:
                    plot_warnings.append(
                        f"{comparison} / {scope_label} / {direction_label}: {suffix} readiness 检查需复核。"
                    )
            except Exception as exc:
                plot_warnings.append(
                    f"{comparison} / {scope_label} / {direction_label}: GO 主题优先 {suffix} 生成失败（{exc}）"
                )
            finally:
                if figure is not None:
                    plt.close(figure)

    source_columns = [
        "comparison_id", "comparison", "deg_scope", "cluster", "direction", "gene_set", "Database", "Term",
        "Overlap", "Adjusted P-value", "Genes", "selection_rank", "selection_reason", "selection_priority", "selection_policy",
    ]
    selected_source = (
        pd.concat(selection_frames, ignore_index=True, sort=False)
        if selection_frames else pd.DataFrame(columns=source_columns)
    )
    source_path = _available_path(package_dirs["go"], f"{output_prefix}_GO_priority_source", ".csv")
    selected_source.to_csv(source_path, index=False)
    result_files.append({
        "file_path": source_path, "file_type": "csv", "category": "table",
        "label": "GO 主题优先三分区图来源表（炎症→脂代谢→Top；全库 FDR）",
    })
    audit_path = _available_path(package_dirs["go"], f"{output_prefix}_GO_priority_selection_audit", ".json")
    Path(audit_path).write_text(json.dumps(selection_audit, ensure_ascii=False, indent=2), encoding="utf-8")
    result_files.append({
        "file_path": audit_path, "file_type": "json", "category": "info",
        "label": "GO 主题优先三分区图选择审计（可复现）",
    })
    return sum(2 for _key, group in groups if str(_key[-1]).upper() == "ORA" and show_plots)


def _emit_pathway_triptych_overviews(analysis, frames, *, show_plots, plot_top_n,
                                     package_dirs, output_prefix, result_files,
                                     figure_audits, plot_warnings):
    """Export an auditable KEGG/Reactome/WikiPathways Human triptych.

    All three libraries must be present for every comparison/cluster/direction
    unit.  Their ORA tables remain separate statistical families; the source
    CSV contains only the selected FDR-ranked display terms.
    """
    if not frames:
        return 0
    combined = pd.concat(frames, ignore_index=True, sort=False)
    observed_libraries = set(combined.get("gene_set", pd.Series(dtype=str)).astype(str))
    required_libraries = set(PATHWAY_TRIPTYCH_DATABASES)
    missing_libraries = sorted(required_libraries - observed_libraries)
    if missing_libraries:
        plot_warnings.append(
            "未生成 KEGG/Reactome/WikiPathways Human 三分区图：本次运行未同时包含 "
            + "、".join(missing_libraries)
        )
        return 0

    from figure_engine import NatureFigureDirector, export_registered_figure
    import matplotlib.pyplot as plt

    director = NatureFigureDirector()
    top_n = min(12, max(1, int(plot_top_n)))
    source_rows = combined.loc[
        combined.get("gene_set", pd.Series("", index=combined.index)).isin(
            PATHWAY_TRIPTYCH_DATABASES
        )
    ].copy()
    group_columns = ["comparison_id", "comparison", "deg_scope", "cluster", "direction", "method"]
    selection_frames = []
    selection_audit = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "statistics_note": (
            "KEGG、Reactome 与 WikiPathways Human 分别在完整库内检验；"
            "Adjusted P-value 为各库内独立 FDR。本审计只记录展示排序，不重新校正或合并 FDR。"
        ),
        "selection_policy": "每个数据库仅展示已显著 term 的 FDR Top 通路。",
        "display_cap_per_database": int(top_n),
        "units": [],
    }
    groups = list(source_rows.groupby(group_columns, dropna=False, sort=True, observed=True))
    for group_key, source_group in groups:
        comparison_id, comparison, deg_scope, cluster, direction, method = group_key
        if str(method).upper() != "ORA":
            continue
        present = set(source_group["gene_set"].dropna().astype(str))
        missing = sorted(required_libraries - present)
        scope_label = _analysis_scope_label(source_group, deg_scope, cluster)
        direction_label = "上调" if str(direction) == "Up" else "下调" if str(direction) == "Down" else str(direction)
        unit_audit = {
            "comparison_id": str(comparison_id), "comparison": str(comparison),
            "deg_scope": str(deg_scope), "cluster": str(cluster), "direction": str(direction),
            "available_libraries": sorted(present),
        }
        if missing:
            message = (
                f"{comparison} / {scope_label} / {direction_label}: 缺少 "
                + "、".join(missing) + " 的完整结果，未生成 Human 通路三分区图。"
            )
            unit_audit.update({"status": "skipped", "reason": message})
            selection_audit["units"].append(unit_audit)
            plot_warnings.append(message)
            continue

        selected, unit_selection_audit = _select_pathway_triptych_rows(
            source_group, top_n=top_n,
        )
        unit_audit.update(unit_selection_audit)
        unit_audit["status"] = "completed"
        selection_audit["units"].append(unit_audit)
        if not selected.empty:
            selection_frames.append(selected)
        if not show_plots:
            continue

        height_mm = min(225, max(132, 62 + 6 * max(len(selected), 6)))
        key_base = (
            f"{output_prefix}_pathway_triptych_{_safe_name(deg_scope)}_"
            f"{_safe_name(comparison_id)}_{_safe_name(cluster, 'All')}_"
            f"{_safe_name(direction, 'All')}"
        )
        title_suffix = f"{comparison} · {scope_label} · {direction_label}"
        for plot_type, suffix, heading in (
            ("pathway_triptych_dotplot", "dotplot", "Human Pathway Enrichment Plot"),
            ("pathway_triptych_barplot", "barplot", "Human pathway enrichment"),
        ):
            figure = None
            try:
                spec = director.create_spec(
                    plot_type, width="double", height_mm=height_mm, top_n=top_n,
                    formats=("svg", "pdf", "png"),
                    database_scope=tuple(PATHWAY_TRIPTYCH_DATABASES.values()),
                    title=f"{heading}\n{title_suffix}",
                )
                figure = director.render(spec, selected)
                exported, report = export_registered_figure(
                    figure, os.path.join(analysis.ensure_plots_dir(), f"{key_base}_{suffix}"), spec,
                    category=plot_type,
                    label=(
                        f"KEGG/Reactome/WikiPathways Human 三分区"
                        f" {('气泡图' if suffix == 'dotplot' else '柱状图')} · "
                        f"{scope_label} · {direction_label}（各库独立 FDR）"
                    ),
                    qa_path=os.path.join(
                        package_dirs["go"], f"{key_base}_{suffix}_nature_readiness.json",
                    ),
                )
                result_files.extend(exported)
                figure_audits.append({
                    "status": "pass" if report.ready else "warning",
                    "nature_readiness_score": report.score,
                    "issues": [issue.message for issue in report.issues],
                    "scope": scope_label,
                    "direction": direction_label,
                    "plot_type": plot_type,
                })
                if not report.ready:
                    plot_warnings.append(
                        f"{comparison} / {scope_label} / {direction_label}: {suffix} readiness 检查需复核。"
                    )
            except Exception as exc:
                plot_warnings.append(
                    f"{comparison} / {scope_label} / {direction_label}: Human 通路三分区 {suffix} 生成失败（{exc}）"
                )
            finally:
                if figure is not None:
                    plt.close(figure)

    source_columns = [
        "comparison_id", "comparison", "deg_scope", "cluster", "direction", "gene_set", "Database", "Term",
        "Overlap", "Adjusted P-value", "Genes", "selection_rank", "selection_reason", "selection_priority", "selection_policy",
    ]
    selected_source = (
        pd.concat(selection_frames, ignore_index=True, sort=False)
        if selection_frames else pd.DataFrame(columns=source_columns)
    )
    source_path = _available_path(package_dirs["go"], f"{output_prefix}_pathway_triptych_source", ".csv")
    selected_source.to_csv(source_path, index=False)
    result_files.append({
        "file_path": source_path, "file_type": "csv", "category": "table",
        "label": "KEGG/Reactome/WikiPathways Human 三分区图来源表（各库独立 FDR Top）",
    })
    audit_path = _available_path(package_dirs["go"], f"{output_prefix}_pathway_triptych_selection_audit", ".json")
    Path(audit_path).write_text(json.dumps(selection_audit, ensure_ascii=False, indent=2), encoding="utf-8")
    result_files.append({
        "file_path": audit_path, "file_type": "json", "category": "info",
        "label": "KEGG/Reactome/WikiPathways Human 三分区图选择审计（可复现）",
    })
    return sum(
        2 for group_key, _source_group in groups
        if show_plots and str(group_key[-1]).upper() == "ORA"
    )


class SCCellGOEnrichment(BaseAnalysis):
    """Run enrichment for cell-level or sample-level DEG exports."""

    MODULE_NAME = "sc_cell_go"
    DISPLAY_NAME = "单细胞 DEG 通路富集"
    DESCRIPTION = "基于单细胞探索性或样本级 pseudobulk DEG，运行 GO/KEGG/Reactome/WikiPathways Human ORA / GSEA"

    def run(self, input_path):
        # Keep this module chained like other analyses, even though it consumes
        # the persisted DEG CSV package rather than expression from the h5ad.
        export_folder = str(self.params.get("export_folder", "sc_batch_results") or "sc_batch_results")
        source_level = str(self.params.get("source_level", "cell_level") or "cell_level").lower()
        if source_level not in {"cell_level", "pseudobulk"}:
            raise ValueError("source_level 必须为 cell_level 或 pseudobulk")
        default_prefix = "sc_pseudobulk" if source_level == "pseudobulk" else "sc_cell_level"
        deg_prefix = _safe_name(
            self.params.get("deg_prefix", default_prefix), default_prefix,
        )
        output_prefix = _safe_name(
            self.params.get("export_prefix", "sc_cell_go"), "sc_cell_go",
        )
        # Register complete pathway tables by default so each enrichment unit
        # is directly downloadable from the task result page.  Callers can
        # still opt out for internal-only, compact pipeline runs.
        export_full_tables = str(self.params.get("export_full_tables", "true")).strip().lower() not in {
            "", "0", "false", "no", "off",
        }
        focus_terms = _parse_focus_terms(self.params.get("focus_terms"))
        focus_label = str(
            self.params.get("focus_label", "炎症与脂代谢优先") or "炎症与脂代谢优先"
        ).strip()
        if len(focus_label) > 80:
            raise ValueError("focus_label 不能超过 80 个字符")
        focus_plot_mode = str(
            self.params.get("focus_plot_mode", "overview_only")
            or "overview_only"
        ).strip().lower()
        if focus_plot_mode not in {"individual_and_overview", "overview_only"}:
            raise ValueError("focus_plot_mode 必须为 individual_and_overview 或 overview_only")
        method = str(self.params.get("method", "ORA") or "ORA").upper()
        if method not in {"ORA", "GSEA"}:
            raise ValueError("method 必须为 ORA 或 GSEA")
        libraries = _requested_gene_sets(self.params)
        go_priority_overview_requested = (
            method == "ORA"
            and focus_plot_mode == "overview_only"
            and set(GO_FOCUS_OVERVIEW_DATABASES).issubset(set(libraries))
        )
        pathway_triptych_overview_requested = (
            method == "ORA"
            and focus_plot_mode == "overview_only"
            and set(PATHWAY_TRIPTYCH_DATABASES).issubset(set(libraries))
        )
        organism = str(self.params.get("organism", "Human") or "Human")
        execution_mode = str(self.params.get("execution_mode", "local") or "local").lower()
        if execution_mode not in {"local", "enrichr"}:
            raise ValueError("execution_mode 必须为 local 或 enrichr")
        if source_level == "pseudobulk" and execution_mode != "local":
            raise ValueError(
                "pseudobulk 富集必须使用 local 基因集，以保留真实 tested-gene background；"
                "Enrichr 在线接口不能满足该统计契约。"
            )
        if method == "GSEA" and source_level != "pseudobulk":
            raise ValueError("正式 GSEA 需要样本级 pseudobulk DEG；细胞级结果仅支持探索性 ORA")
        if method == "GSEA" and execution_mode != "local":
            raise ValueError("GSEA 仅支持 local 基因集；不能把完整排序表发送到 Enrichr")
        if organism.lower() != "human":
            raise ValueError("当前富集统计仅支持 Human；小鼠数据需先完成物种映射")
        local_gene_set_dir = _resolve_local_gene_set_dir(
            self.params.get("local_gene_set_dir", ""), self.project_dir,
        )
        padj_cutoff = float(self.params.get("padj_cutoff", 0.05))
        log2fc_cutoff = float(self.params.get("log2fc_cutoff", 0.25))
        if not np.isfinite(padj_cutoff) or not 0 < padj_cutoff <= 1:
            raise ValueError("padj_cutoff 必须位于 (0, 1]")
        if not np.isfinite(log2fc_cutoff) or log2fc_cutoff < 0:
            raise ValueError("log2fc_cutoff 必须是非负有限数值")
        min_genes = max(3, int(self.params.get("min_genes", 5)))
        top_n = max(1, int(self.params.get("top_n", 30)))
        direction_mode = str(
            self.params.get(
                "direction_mode",
                "up_down" if source_level == "pseudobulk" else "up",
            ) or "up"
        ).lower()
        if direction_mode not in {"up", "down", "up_down", "all"}:
            raise ValueError("direction_mode 必须为 up、down、up_down 或 all")
        if method == "GSEA":
            directions = ["All"]
        elif direction_mode == "up_down":
            directions = ["Up", "Down"]
        elif direction_mode == "down":
            directions = ["Down"]
        elif direction_mode == "all":
            directions = ["All"]
        else:
            directions = ["Up"]
        ranking_metric = str(self.params.get("ranking_metric", "auto") or "auto")
        ora_min_size = max(1, int(self.params.get("ora_min_size", 5)))
        ora_max_size = max(ora_min_size, int(self.params.get("ora_max_size", 500)))
        gsea_min_size = max(2, int(self.params.get("gsea_min_size", 15)))
        gsea_max_size = max(gsea_min_size, int(self.params.get("gsea_max_size", 500)))
        permutation_num = min(5000, max(100, int(self.params.get("permutation_num", 1000))))
        gsea_seed = int(self.params.get("gsea_seed", 112))
        gsea_weight = float(self.params.get("gsea_weight", 1.0))
        if gsea_seed < 0:
            raise ValueError("gsea_seed 必须为非负整数")
        if not np.isfinite(gsea_weight) or gsea_weight < 0:
            raise ValueError("gsea_weight 必须是非负有限数值")
        show_enrichment_plots = self.params.get("show_enrichment_plots", True)
        if isinstance(show_enrichment_plots, str):
            show_enrichment_plots = show_enrichment_plots.strip().lower() not in {
                "", "0", "false", "no", "off",
            }
        plot_top_n = max(1, min(30, int(self.params.get("plot_top_n", 8))))
        allocation_mode, priority_allocation = _go_priority_allocation_from_params(
            self.params, top_n=plot_top_n,
        )
        plot_cluster_limit = int(self.params.get("plot_max_clusters", 0))
        if plot_cluster_limit < 0:
            raise ValueError("plot_max_clusters 必须为 0 或正整数")
        plot_suite = str(
            self.params.get("enrichment_plot_suite", "完整 Nature 套图")
            or "完整 Nature 套图"
        )
        if plot_suite.lower() not in {
            "完整 nature 套图", "full", "complete", "核心图", "核心", "core", "minimal",
        }:
            raise ValueError("enrichment_plot_suite 必须为完整 Nature 套图或核心图")
        enrichment_max_genes = min(
            80, max(5, int(self.params.get("enrichment_max_genes", 30)))
        )
        similarity_threshold = float(self.params.get("similarity_threshold", 0.15))
        if not np.isfinite(similarity_threshold) or not 0 <= similarity_threshold <= 1:
            raise ValueError("similarity_threshold 必须位于 [0, 1]")
        running_term_n = min(4, max(1, int(self.params.get("running_term_n", 2))))

        package_root, package_dirs = batch_csv_package_dirs(
            self.project_dir, export_folder, included_keys=("deg", "go"),
        )
        selected_sources = self.params.get("deg_source_files", [])
        source_files = _selected_deg_files(selected_sources, self.project_dir)
        source_selection = "explicit_task_binding" if source_files else "legacy_package_discovery"
        if not source_files:
            source_files = _latest_deg_files(package_dirs["deg"], deg_prefix)
        if not source_files:
            raise ValueError(
                "未找到可绑定的 DEG 结果；请先完成细胞类型/cluster marker、细胞级或 pseudobulk 差异分析，"
                "并在富集页面选择该任务。"
            )
        requested_clusters = _requested_clusters(self.params.get("target_clusters", ""))
        available_clusters = set()
        for source_path in source_files.values():
            try:
                header = pd.read_csv(source_path, nrows=0)
                if "cluster" not in header.columns:
                    available_clusters.add("All")
                    continue
                cluster_values = pd.read_csv(source_path, usecols=["cluster"])["cluster"]
                available_clusters.update(
                    str(value).strip() for value in cluster_values.dropna()
                    if str(value).strip()
                )
            except (OSError, ValueError, pd.errors.ParserError) as exc:
                raise ValueError(f"无法读取 DEG cluster 范围: {source_path.name}") from exc
        unknown_clusters = sorted(set(requested_clusters) - available_clusters)
        if unknown_clusters:
            raise ValueError(
                "所选 DEG 任务中不存在以下 cluster：" + "、".join(unknown_clusters[:5])
            )
        internal_go_dir = os.path.join(package_dirs["go"], ".internal")
        os.makedirs(internal_go_dir, exist_ok=True)
        gp = None
        if execution_mode == "enrichr" or source_level == "cell_level":
            try:
                import gseapy as gp
            except ImportError as exc:
                raise ImportError("缺少 gseapy，无法运行细胞级/Enrichr ORA") from exc

        result_files = []
        statuses = []
        plot_warnings = []
        figure_audits = []
        go_focus_overview_frames = []
        pathway_triptych_frames = []
        provenance_warnings = []
        source_groupings = []
        if source_selection == "legacy_package_discovery":
            provenance_warnings.append(
                "使用旧版结果包发现模式；后续请从富集页面明确选择 DEG 任务以固定溯源。"
            )
        plots_dir = self.ensure_plots_dir() if show_enrichment_plots else None
        total_outputs = max(1, len(source_files) * len(libraries))
        output_index = 0
        for (comparison_id, deg_scope), source_path in sorted(source_files.items()):
            deg = pd.read_csv(source_path)
            if deg.empty:
                raise ValueError(f"所选 DEG 结果为空: {source_path.name}")
            source_comparison_ids = set(deg["comparison_id"].dropna().astype(str))
            if deg["comparison_id"].isna().any() or source_comparison_ids != {str(comparison_id)}:
                raise ValueError(
                    f"DEG 文件混有多个或不匹配的 comparison_id: {source_path.name}"
                )
            if "deg_scope" in deg.columns:
                source_scopes = set(deg["deg_scope"].dropna().astype(str))
                if deg["deg_scope"].isna().any() or source_scopes != {str(deg_scope)}:
                    raise ValueError(
                        f"DEG 文件混有多个或不匹配的 deg_scope: {source_path.name}"
                    )
            inference_units = set(deg.get("inference_unit", pd.Series(dtype=str)).astype(str))
            known_inference_units = inference_units.intersection({"cell", "biological_sample"})
            if len(known_inference_units) > 1:
                raise ValueError(f"DEG 文件混用了 cell 与 biological_sample 推断单位: {source_path.name}")
            inferred_level = (
                "pseudobulk" if "biological_sample" in inference_units else
                "cell_level" if "cell" in inference_units else ""
            )
            declared_level = str(self.params.get("deg_source_level", "") or "").strip().lower()
            if declared_level and inferred_level and declared_level != inferred_level:
                raise ValueError(
                    f"所选 DEG 任务层级为 {inferred_level}，但提交参数声明为 {declared_level}；"
                    "为避免混用，富集已停止。"
                )
            if inferred_level and inferred_level != source_level:
                raise ValueError(
                    f"所选 DEG 结果为 {inferred_level}，但当前 source_level 为 {source_level}；"
                    "请重新选择对应的 DEG 任务。"
                )
            comparison = _single_text_value(deg, ("comparison",), comparison_id)
            experimental = _single_text_value(
                deg, ("experimental_group", "group_1"), "",
            )
            control = _single_text_value(
                deg, ("control_group", "group_2"), "",
            )
            source_statistical_status = _single_text_value(
                deg,
                ("statistical_status",),
                "biological_replicate_model" if source_level == "pseudobulk"
                else "exploratory_no_biological_replicates",
            )
            source_analysis_grouping = _single_text_value(
                deg, ("analysis_grouping",), "",
            )
            source_analysis_group_key = _single_text_value(
                deg, ("analysis_group_key",), "",
            )
            source_analysis_group_label = _single_text_value(
                deg, ("analysis_group_label",), "cluster",
            )
            source_groupings.append({
                "mode": source_analysis_grouping,
                "key": source_analysis_group_key,
                "label": source_analysis_group_label,
            })
            gene_column = "gene" if "gene" in deg.columns else "names"
            padj_column = next(
                (column for column in ("padj", "p.adjust", "p_val_adj") if column in deg.columns),
                None,
            )
            fc_column = next(
                (column for column in ("log2FC", "avg_log2FC", "logFC") if column in deg.columns),
                None,
            )
            if gene_column not in deg.columns or padj_column is None or fc_column is None:
                raise ValueError(
                    f"DEG 文件缺少 gene/{fc_column or 'log2FC'}/{padj_column or 'padj'} 列: "
                    f"{source_path.name}"
                )
            if deg_scope == "per_cluster" and "cluster" in deg.columns:
                units = [(str(cluster), table.copy()) for cluster, table in deg.groupby("cluster", sort=True)]
            elif "cluster" in deg.columns and deg["cluster"].dropna().astype(str).nunique() == 1:
                units = [(_single_text_value(deg, ("cluster",), "All"), deg)]
            else:
                units = [("All", deg)]
            if requested_clusters:
                requested_set = set(requested_clusters)
                units = [
                    (cluster, table) for cluster, table in units
                    if str(cluster) in requested_set
                ]
            # A per-cluster selection may deliberately omit a companion
            # all-cells DEG source from the same task.  It should produce no
            # empty table for that excluded scope.
            if not units:
                continue
            for library in libraries:
                output_index += 1
                gene_set_path = None
                local_pathways = None
                gene_set_error = ""
                if execution_mode == "local":
                    try:
                        gene_set_path = _local_gene_set_path(library, local_gene_set_dir)
                        local_pathways = load_human_genesets(gene_set_path)
                    except Exception as exc:
                        gene_set_error = str(exc)
                self.progress(
                    5 + int((output_index - 1) * 80 / total_outputs),
                    f"正在运行 {comparison_id} / {deg_scope} / {library} 的 {method} 富集...",
                )
                scope_results = []
                mapping_qc_by_cluster = {}
                gsea_ranking_qc_by_cluster = {}
                gsea_results_by_cluster = {}
                ora_testing_family_by_unit = {}
                for cluster, unit_deg in units:
                    unit_deg = unit_deg.copy()
                    fc = pd.to_numeric(unit_deg[fc_column], errors="coerce")
                    padj = pd.to_numeric(unit_deg[padj_column], errors="coerce")
                    # ORA's universe is the set of genes with a finite effect
                    # and finite adjusted test result in this exact DEG unit.
                    # Rows that never produced a usable test must not inflate N.
                    tested_mask = pd.Series(
                        np.isfinite(fc.to_numpy(dtype=float))
                        & np.isfinite(padj.to_numpy(dtype=float)),
                        index=unit_deg.index,
                    )
                    raw_background = unit_deg.loc[
                        tested_mask, gene_column
                    ].dropna().astype(str).tolist()
                    background = list(dict.fromkeys(item for item in raw_background if item.strip()))
                    prepared = None
                    identifier_map = {}
                    if source_level == "pseudobulk":
                        prepared, identifier_map, mapping_qc = prepare_deg_gene_table(unit_deg)
                        mapping_qc["n_tested_rows_for_ora_background"] = int(tested_mask.sum())
                        mapping_qc_by_cluster[str(cluster)] = mapping_qc
                        background = prepared.loc[
                            tested_mask, "_analysis_gene"
                        ].dropna().astype(str).tolist()
                        background = list(dict.fromkeys(item for item in background if item.strip()))

                    ranking = None
                    ranking_error = ""
                    if method == "GSEA":
                        try:
                            ranking, ranking_qc = prepare_gsea_ranking(
                                prepared, ranking_metric=ranking_metric,
                                identifier_map=identifier_map,
                            )
                            gsea_ranking_qc_by_cluster[str(cluster)] = ranking_qc
                        except Exception as exc:
                            ranking_error = str(exc)

                    for direction in directions:
                        if method == "GSEA":
                            selected = (
                                ranking["gene_name"].astype(str).tolist()
                                if ranking is not None else []
                            )
                            n_input = len(selected)
                        else:
                            if direction == "Up":
                                selected_mask = (padj <= padj_cutoff) & (fc >= log2fc_cutoff)
                            elif direction == "Down":
                                selected_mask = (padj <= padj_cutoff) & (fc <= -log2fc_cutoff)
                            else:
                                selected_mask = (padj <= padj_cutoff) & (fc.abs() >= log2fc_cutoff)
                            selected_source = (
                                prepared["_analysis_gene"]
                                if source_level == "pseudobulk" else unit_deg[gene_column]
                            )
                            selected = selected_source.loc[
                                selected_mask & tested_mask
                            ].dropna().astype(str).tolist()
                            selected = list(dict.fromkeys(item for item in selected if item.strip()))
                            n_input = len(selected)
                        analysis_level = (
                            "sample_level_pseudobulk" if source_level == "pseudobulk"
                            else "cell_level_exploratory"
                        )
                        statistical_status = source_statistical_status
                        base = {
                            "comparison_id": comparison_id, "comparison": comparison,
                            "experimental_group": experimental, "control_group": control,
                            "analysis_level": analysis_level,
                            "statistical_status": statistical_status,
                            "deg_scope": deg_scope, "cluster": cluster,
                            "analysis_grouping": source_analysis_grouping,
                            "analysis_group_key": source_analysis_group_key,
                            "analysis_group_label": source_analysis_group_label,
                            "direction": direction, "gene_set": library,
                            "n_input_genes": n_input,
                            "n_background_genes": (
                                n_input if method == "GSEA" else len(set(background))
                            ),
                            "gene_universe": (
                                "complete_valid_ranking" if method == "GSEA"
                                else "finite_effect_and_adjusted_test_rows"
                            ),
                            "method": method,
                        }
                        if method == "GSEA":
                            runnable = len(selected) >= min_genes
                            not_runnable_reason = (
                                f"完整排序基因少于 {min_genes} 个，无法进行 GSEA"
                            )
                        else:
                            runnable = len(selected) >= min_genes
                            not_runnable_reason = (
                                f"{direction} 显著基因少于 {min_genes} 个，无法进行稳定 ORA"
                            )
                        if gene_set_error:
                            result = pd.DataFrame([{
                                **base, "status": "failed",
                                "reason": f"本地基因集不可用: {gene_set_error}",
                            }])
                        elif ranking_error:
                            result = pd.DataFrame([{
                                **base, "status": "failed",
                                "reason": f"GSEA 排名表无效: {ranking_error}",
                            }])
                        elif not runnable:
                            result = pd.DataFrame([{
                                **base, "status": "not_runnable", "reason": not_runnable_reason,
                            }])
                        else:
                            try:
                                if source_level == "pseudobulk" and execution_mode == "local":
                                    if method == "ORA":
                                        enrichment = run_ora_full(
                                            selected, background, local_pathways,
                                            pvalue_cutoff=padj_cutoff,
                                            min_size=ora_min_size, max_size=ora_max_size,
                                        )
                                        ora_testing_family_by_unit[
                                            f"{cluster}|{direction}"
                                        ] = dict(enrichment.attrs.get("testing_family", {}))
                                        enrichment["gene_set"] = library
                                    else:
                                        pre_res = run_gsea_prerank(
                                            ranking, local_pathways, permutation_num=permutation_num,
                                            seed=gsea_seed, weight=gsea_weight,
                                            min_size=gsea_min_size, max_size=gsea_max_size,
                                        )
                                        gsea_results_by_cluster[str(cluster)] = pre_res
                                        enrichment = pre_res.res2d.copy()
                                        if "Term" not in enrichment.columns:
                                            index_name = enrichment.index.name or "index"
                                            enrichment = enrichment.reset_index().rename(
                                                columns={index_name: "Term"}
                                            )
                                        fdr_column = next(
                                            (column for column in ("fdr", "FDR", "FDR q-val")
                                             if column in enrichment.columns), None,
                                        )
                                        if fdr_column is None:
                                            raise ValueError("GSEA 结果缺少 FDR 列")
                                        enrichment["Adjusted P-value"] = pd.to_numeric(
                                            enrichment[fdr_column], errors="coerce"
                                        )
                                        enrichment["Significant"] = (
                                            enrichment["Adjusted P-value"] < padj_cutoff
                                        )
                                        lead_column = next(
                                            (column for column in ("lead_genes", "Lead_genes",
                                                                   "leading_edge", "core_enrichment")
                                             if column in enrichment.columns), None,
                                        )
                                        if lead_column and "Genes" not in enrichment.columns:
                                            enrichment["Genes"] = enrichment[lead_column]
                                        enrichment = add_gsea_leading_edge_metrics(enrichment)
                                        enrichment["gene_set"] = library
                                    if enrichment is None or enrichment.empty:
                                        result = pd.DataFrame([{
                                            **base, "status": "empty", "reason": "没有匹配的通路",
                                        }])
                                    else:
                                        result = enrichment.copy()
                                        if method == "ORA":
                                            result = result.sort_values(
                                                ["Adjusted P-value", "P-value"], kind="stable"
                                            )
                                        # Pseudobulk enrichment is the formal
                                        # inferential branch: retain every
                                        # tested pathway so the adjusted-pvalue
                                        # family is auditable. ``top_n`` only
                                        # limits legacy cell-level exports and
                                        # the plot renderer.
                                        if source_level == "cell_level":
                                            result = result.head(top_n).copy()
                                        for key, value in base.items():
                                            result[key] = value
                                        result["status"] = "completed"
                                        result["reason"] = ""
                                else:
                                    enrichment = _run_ora(
                                        gp, selected, library, organism, execution_mode,
                                        local_gene_set_dir, background=background,
                                    )
                                    if enrichment is None or enrichment.empty:
                                        result = pd.DataFrame([{
                                            **base, "status": "empty", "reason": "未返回显著/匹配 GO term",
                                        }])
                                    else:
                                        result = enrichment.head(top_n).copy()
                                        for key, value in base.items():
                                            result[key] = value
                                        result["status"] = "completed"
                                        result["reason"] = ""
                            except Exception as exc:
                                source = "本地基因集" if execution_mode == "local" else "GO 服务"
                                result = pd.DataFrame([{
                                    **base, "status": "failed", "reason": f"{source}不可用: {exc}",
                                }])
                        scope_results.append(result)
                result = pd.concat(scope_results, ignore_index=True)
                adjusted = pd.to_numeric(
                    result.get("Adjusted P-value", pd.Series(np.nan, index=result.index)),
                    errors="coerce",
                )
                if "Significant" in result.columns:
                    existing_significant = result["Significant"]
                    if existing_significant.dtype == bool:
                        significant = existing_significant.fillna(False)
                    else:
                        significant = existing_significant.astype(str).str.lower().isin(
                            {"true", "1", "yes"}
                        )
                    result["Significant"] = significant | adjusted.lt(padj_cutoff)
                else:
                    result["Significant"] = adjusted.lt(padj_cutoff)
                for column in GO_RESULT_COLUMNS:
                    if column not in result.columns:
                        result[column] = ""
                result = result[GO_RESULT_COLUMNS + [
                    column for column in result.columns if column not in GO_RESULT_COLUMNS
                ]]
                output_path = _available_path(
                    package_dirs["go"] if export_full_tables else internal_go_dir,
                    f"{output_prefix}_go_{_safe_name(library)}_{_safe_name(deg_scope)}_"
                    f"{_safe_name(comparison_id, 'comparison')}", ".csv",
                )
                result.to_csv(output_path, index=False)
                label = (
                    f"Cell-level GO enrichment ({library}, {deg_scope}): {comparison}"
                    if source_level == "cell_level" and method == "ORA" else
                    f"{'Cell-level' if source_level == 'cell_level' else 'Pseudobulk'} "
                    f"{method} enrichment ({library}, {deg_scope}): {comparison}"
                )
                if export_full_tables:
                    result_files.append({
                        "file_path": output_path, "file_type": "csv", "category": "table",
                        "label": label,
                    })
                audit = {
                    "source_file": source_path.name,
                    "source_task_id": str(self.params.get("source_task_id", "") or ""),
                    "source_selection": source_selection,
                    "source_level": source_level,
                    "method": method,
                    "comparison_id": comparison_id,
                    "comparison": comparison,
                    "experimental_group": experimental,
                    "control_group": control,
                    "deg_scope": deg_scope,
                    "analysis_grouping": source_analysis_grouping,
                    "analysis_group_key": source_analysis_group_key,
                    "analysis_group_label": source_analysis_group_label,
                    "selected_clusters": requested_clusters,
                    "direction_mode": direction_mode,
                    "padj_cutoff": padj_cutoff,
                    "log2fc_cutoff": log2fc_cutoff,
                    "ora_min_size": ora_min_size,
                    "ora_max_size": ora_max_size,
                    "gsea_min_size": gsea_min_size,
                    "gsea_max_size": gsea_max_size,
                    "permutation_num": permutation_num,
                    "gsea_seed": gsea_seed,
                    "ranking_metric": ranking_metric,
                    "execution_mode": execution_mode,
                    "gene_set": library,
                    "identifier_mapping": mapping_qc_by_cluster,
                    "ora_testing_family": ora_testing_family_by_unit,
                    "gsea_ranking_qc": gsea_ranking_qc_by_cluster,
                    "units": [],
                }
                if execution_mode == "local" and gene_set_path is not None and local_pathways:
                    audit["gene_set_manifest"] = geneset_manifest(
                        gene_set_path, local_pathways,
                        source_type="project_or_platform_local",
                    )
                elif execution_mode == "local":
                    audit["gene_set_error"] = gene_set_error or "本地基因集不可用"
                for (audit_cluster, audit_direction), audit_frame in result.groupby(
                    ["cluster", "direction"], dropna=False, observed=True
                ):
                    first = audit_frame.iloc[0]
                    audit["units"].append({
                        "cluster": str(audit_cluster),
                        "direction": str(audit_direction),
                        "status": str(first.get("status", "")),
                        "n_input_genes": _safe_int(first.get("n_input_genes", 0)),
                        "n_background_genes": _safe_int(first.get("n_background_genes", 0)),
                        "n_terms_exported": (
                            int(len(audit_frame))
                            if str(first.get("status", "")) == "completed" else 0
                        ),
                        "n_result_rows": int(len(audit_frame)),
                        "n_significant_terms": int(
                            audit_frame.get(
                                "Significant", pd.Series(False, index=audit_frame.index)
                            ).astype(bool).sum()
                        ),
                    })
                _emit_focus_outputs(
                    self, result, focus_terms=focus_terms,
                    show_plots=(show_enrichment_plots and focus_plot_mode != "overview_only"),
                    plot_top_n=plot_top_n,
                    enrichment_max_genes=enrichment_max_genes,
                    similarity_threshold=similarity_threshold,
                    running_term_n=running_term_n,
                    package_dirs=package_dirs, output_prefix=output_prefix,
                    library=library, deg_scope=deg_scope,
                    comparison_id=comparison_id, comparison=comparison,
                    method=method, source_level=source_level,
                    padj_cutoff=padj_cutoff,
                    gsea_results_by_cluster=gsea_results_by_cluster,
                    audit=audit, result_files=result_files,
                    figure_audits=figure_audits, plot_warnings=plot_warnings,
                )
                if library in GO_FOCUS_OVERVIEW_DATABASES:
                    go_focus_overview_frames.append(result)
                if library in PATHWAY_TRIPTYCH_DATABASES:
                    pathway_triptych_frames.append(result)
                audit_path = _available_path(
                    package_dirs["go"],
                    f"{output_prefix}_enrichment_{_safe_name(library)}_"
                    f"{_safe_name(deg_scope)}_{_safe_name(comparison_id, 'comparison')}_audit", ".json",
                )
                Path(audit_path).write_text(
                    json.dumps(audit, ensure_ascii=False, indent=2, default=str),
                    encoding="utf-8",
                )
                result_files.append({
                    "file_path": audit_path, "file_type": "json", "category": "qc",
                    "label": f"Enrichment audit: {comparison} / {library} / {deg_scope}",
                })
                suppress_individual_view = (
                    (library in GO_FOCUS_OVERVIEW_DATABASES and go_priority_overview_requested)
                    or (library in PATHWAY_TRIPTYCH_DATABASES and pathway_triptych_overview_requested)
                )
                if show_enrichment_plots and not suppress_individual_view:
                    completed = result.loc[result["status"].eq("completed")].copy()
                    if not completed.empty:
                        from modules.bulk_enrichment import _render_enrichment_variants

                        significant = completed["Significant"]
                        if significant.dtype != bool:
                            significant = significant.astype(str).str.lower().isin(
                                {"true", "1", "yes"}
                            )
                        plot_completed = completed.loc[significant].copy()
                        if plot_completed.empty:
                            plot_warnings.append(
                                f"{comparison} / {library}: 没有 FDR<{padj_cutoff:g} 的通路，"
                                "不生成可能误导的富集图；完整检验结果保留在内部表。"
                            )
                            plot_completed = completed.iloc[0:0].copy()
                        cluster_priority = (
                            plot_completed.groupby("cluster", observed=True)["Adjusted P-value"]
                            .min().sort_values().index.astype(str).tolist()
                        )
                        clusters_to_plot = (
                            cluster_priority if plot_cluster_limit == 0
                            else cluster_priority[:plot_cluster_limit]
                        )
                        if len(cluster_priority) > len(clusters_to_plot):
                            plot_warnings.append(
                                f"{comparison} / {library}: GO 图仅展示最小 FDR 最优的前 "
                                f"{plot_cluster_limit} 个细胞类型 / cluster；完整结果保留在受控内部表。"
                            )
                        for cluster in clusters_to_plot:
                            level_title = (
                                "Cell-level exploratory GO" if source_level == "cell_level"
                                else f"Pseudobulk {method}"
                            )
                            cluster_frame = plot_completed.loc[
                                plot_completed["cluster"].astype(str) == cluster
                            ].copy()
                            scope_label = _analysis_scope_label(cluster_frame, deg_scope, cluster)
                            cluster_suffix = "" if deg_scope == "all_cells" else f" · {scope_label}"
                            for plot_direction, plot_frame in cluster_frame.groupby(
                                "direction", sort=False, observed=True
                            ):
                                if plot_frame.empty:
                                    continue
                                direction_suffix = (
                                    "" if method == "GSEA" else f" · {plot_direction}"
                                )
                                title = (
                                    f"{level_title} · {comparison}{cluster_suffix}{direction_suffix}"
                                )
                                stem = (
                                    f"sc_cell_go_{_safe_name(library)}_{_safe_name(deg_scope)}_"
                                    f"{_safe_name(comparison_id)}_{_safe_name(cluster, 'All')}_"
                                    f"{_safe_name(plot_direction, 'All')}"
                                )
                                plot_params = dict(self.params)
                                plot_params.update({
                                    "top_n": plot_top_n,
                                    "enrichment_plot_suite": plot_suite,
                                    "enrichment_max_genes": enrichment_max_genes,
                                    "similarity_threshold": similarity_threshold,
                                    "running_term_n": running_term_n,
                                })
                                semantic_warnings = ()
                                if source_level == "cell_level":
                                    semantic_warnings = (
                                        "该富集基于细胞级探索性 DEG；细胞不是独立生物学重复，"
                                        "正式条件比较应以样本级 pseudobulk 结果验证。",
                                    )
                                try:
                                    exported, audits = _render_enrichment_variants(
                                        self, plot_frame, method=method, title=title,
                                        base_output_key=stem,
                                        label=(
                                            f"{level_title}: {comparison}{cluster_suffix}"
                                            f"{direction_suffix}"
                                        ),
                                        params=plot_params,
                                        pre_res=gsea_results_by_cluster.get(str(cluster)),
                                        semantic_warnings=semantic_warnings,
                                    )
                                    result_files.extend(exported)
                                    figure_audits.extend(audits)
                                    for figure_audit in audits:
                                        if figure_audit.get("status") == "warning":
                                            issues = "; ".join(
                                                str(item) for item in figure_audit.get("issues", [])[:2]
                                            ) or "图形 readiness 检查未完全通过"
                                            plot_warnings.append(
                                                f"{comparison} / {library} / cluster {cluster} / "
                                                f"{plot_direction}: {issues}"
                                            )
                                except Exception as exc:
                                    plot_warnings.append(
                                        f"{comparison} / {library} / {scope_label} / "
                                        f"{plot_direction}: 富集绘图失败（{exc}）"
                                    )
                statuses.append({
                    "comparison_id": comparison_id, "comparison": comparison,
                    "gene_set": library, "deg_scope": deg_scope, "method": method,
                    "source_level": source_level,
                    "execution_mode": execution_mode,
                    "status": "completed" if result["status"].eq("completed").any() else str(result["status"].iloc[0]),
                    "output_file": os.path.basename(output_path),
                    "n_clusters": int(result["cluster"].nunique()),
                })

        n_go_focus_overviews = _emit_go_focus_overviews(
            self, go_focus_overview_frames,
            focus_terms=focus_terms, focus_label=focus_label,
            show_plots=show_enrichment_plots, plot_top_n=plot_top_n,
            priority_allocation=priority_allocation, allocation_mode=allocation_mode,
            package_dirs=package_dirs, output_prefix=output_prefix,
            result_files=result_files, figure_audits=figure_audits,
            plot_warnings=plot_warnings,
        )
        n_pathway_triptych_overviews = 0
        if method == "ORA":
            n_pathway_triptych_overviews = _emit_pathway_triptych_overviews(
                self, pathway_triptych_frames,
                show_plots=show_enrichment_plots, plot_top_n=plot_top_n,
                package_dirs=package_dirs, output_prefix=output_prefix,
                result_files=result_files, figure_audits=figure_audits,
                plot_warnings=plot_warnings,
            )
        if export_full_tables:
            manifest_path = _available_path(
                package_dirs["go"], f"{output_prefix}_go_comparison_manifest", ".csv"
            )
            pd.DataFrame(statuses).to_csv(manifest_path, index=False)
            result_files.append({
                "file_path": manifest_path, "file_type": "csv", "category": "table",
                "label": "Cell/pseudobulk enrichment comparison registry",
            })
        self.progress(100, "单细胞差异基因富集分析完成")
        all_failed = bool(statuses) and all(item.get("status") == "failed" for item in statuses)
        summary = {
            "n_comparisons": int(len(statuses)), "csv_package_dir": package_root,
            "analysis_level": (
                "cell_level_exploratory" if source_level == "cell_level"
                else "sample_level_pseudobulk"
            ),
            "source_level": source_level, "method": method, "direction_mode": direction_mode,
            "gene_sets": libraries, "organism": organism, "execution_mode": execution_mode,
            "enrichment_plot_suite": plot_suite,
            "n_figure_audits": int(len(figure_audits)),
            "plot_warnings": plot_warnings[:20],
            "provenance_warnings": provenance_warnings,
            "source_selection": source_selection,
            "source_task_id": str(self.params.get("source_task_id", "") or ""),
            "source_analysis_groupings": list({
                (item["mode"], item["key"], item["label"]): item
                for item in source_groupings
            }.values()),
            "selected_clusters": requested_clusters,
            "cluster_execution_scope": "selected" if requested_clusters else "all",
            "full_table_exported": export_full_tables,
        }
        summary["go_priority_display"] = {
            "allocation_mode": allocation_mode,
            "requested_allocation": priority_allocation,
            "display_cap_per_ontology": int(plot_top_n),
            "overview_requested": go_priority_overview_requested,
        }
        summary["pathway_triptych_display"] = {
            "libraries": list(PATHWAY_TRIPTYCH_DATABASES.values()),
            "display_cap_per_database": int(min(12, plot_top_n)),
            "overview_requested": pathway_triptych_overview_requested,
            "n_overviews": int(n_pathway_triptych_overviews),
            "statistics_note": "各 Human 通路库保留独立完整库内 FDR；三分区图不合并或重新校正 FDR。",
        }
        summary["n_go_focus_overviews"] = int(n_go_focus_overviews)
        if focus_terms:
            summary["focus_terms"] = focus_terms
            summary["focus_note"] = "focus_terms 仅增加已确认主题的展示优先级；全量统计与全量表保持不变。"
            summary["focus_label"] = focus_label
        summary["focus_plot_mode"] = focus_plot_mode
        if all_failed:
            summary["error"] = "所有指定富集单元均执行失败；请查看通路审计中的基因集与错误信息。"
        return {
            "output_adata": input_path,
            "result_files": result_files,
            "summary": summary,
        }
