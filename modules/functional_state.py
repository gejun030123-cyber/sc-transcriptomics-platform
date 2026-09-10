"""Sample-aware regulatory and functional-state analysis for scRNA-seq.

This module deliberately keeps three biological layers separate:

* gene expression is a direct value from a log-normalised expression matrix;
* TF activity is an optional, *network-backed* signed target score;
* pathway activity is a control-gene-adjusted expression score.

It never turns expression of a transcription-factor gene into TF activity, and
it aggregates cells to sample × cell-type units before reporting a group p
value.  Existing ``sc_pseudobulk_deg`` and ``sc_cell_go`` remain the sole
sources of count-level condition DEG and GSEA evidence.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import textwrap
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from modules.base import BaseAnalysis
from modules.figure_style import NATURE_PALETTE, nature_continuous_cmap
from modules.gene_set_registry import load_selected_managed_gene_sets
from modules.native_figures import heatmap_figure, marker_dotplot_figure
from modules.sc_de_utils import _matrix_is_raw_counts, log1p_expression_from_counts


# These are compact, transparent panels, not a substitute for a versioned
# Hallmark/Reactome library.  Researchers can add larger project-local sets
# through ``custom_gene_sets`` and their exact contents are written to the
# manifest on every run.
METABOLIC_TF_GENES = (
    "PPARA", "RXRA", "HNF4A", "ESRRA", "NRF1", "FOXA2", "KLF15", "PPARD",
    "PPARG", "NFE2L2", "TFEB", "CREB1", "CREB3L3", "HNF1A", "FOXO1",
)
INFLAMMATORY_TF_GENES = ("RELA", "NFKB1", "STAT1", "STAT3", "IRF1", "IRF7", "JUN", "FOS")
METABOLIC_EFFECTOR_GENES = (
    "CPT1A", "CPT2", "SLC25A20", "ACADM", "ACADVL", "ACADS", "HADHA",
    "HADHB", "ECHS1", "ETFDH", "ACAA2", "ACOX1", "EHHADH", "ACAA1",
    "HSD17B4", "ABCD3", "SCP2", "FABP1", "FABP2", "SLC27A2", "HMGCS2",
    "HMGCL", "BDH1",
)
# The generic inflammation programme and the TNF--NF-kB programme must remain
# distinct.  They deliberately have some shared downstream chemokines, but are
# not substitutes: the former emphasises cytokine/chemokine output and the
# latter contains proximal canonical NF-kB feedback and signalling targets.
# Keeping these definitions separate avoids producing duplicate scores that
# would otherwise create spurious duplicated correlation/concordance columns.
INFLAMMATORY_GENES = (
    "IL1B", "IL6", "CCL2", "CXCL8", "CXCL2", "CXCL3", "CCL20", "CXCL1",
    "CXCL5", "S100A8", "S100A9",
)
TNF_NFKB_RESPONSE_GENES = (
    "TNFAIP3", "NFKBIA", "NFKBIZ", "TNIP1", "NFKB2", "RELB", "BIRC3",
    "TRAF1", "ICAM1", "CCL20", "CXCL8", "CXCL2", "CXCL3",
)
GENE_EXPRESSION_PANEL_SECTIONS = {
    "Regulatory TF genes": METABOLIC_TF_GENES,
    "Carnitine shuttle": ("CPT1A", "CPT2", "SLC25A20"),
    "Mitochondrial beta-oxidation": ("ACADM", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "ETFDH", "ACAA2"),
    "Peroxisomal lipid metabolism (FAO)": ("ACOX1", "EHHADH", "ACAA1", "HSD17B4", "ABCD3", "SCP2"),
    "Intracellular epithelial FA transport": ("FABP1", "FABP2", "SLC27A2"),
    "Ketogenesis": ("HMGCS2", "HMGCL", "BDH1"),
    "Acute myeloid chemokine inflammation": INFLAMMATORY_GENES,
}
DEFAULT_EXPRESSION_GENES = tuple(dict.fromkeys(
    gene for genes in GENE_EXPRESSION_PANEL_SECTIONS.values() for gene in genes
))

PPARA_TARGET_MODULE_NAME = "PPARA-associated lipid-oxidation programme"
CARNITINE_SHUTTLE_MODULE_NAME = "Carnitine shuttle (mitochondrial FA entry)"
PEROXISOMAL_LIPID_METABOLISM_MODULE_NAME = "Peroxisomal lipid metabolism (FAO)"
INTRACELLULAR_FA_TRANSPORT_MODULE_NAME = "Intracellular epithelial FA transport"
ACUTE_MYELOID_INFLAMMATION_MODULE_NAME = "Acute myeloid chemokine inflammation"
TYPE_I_II_IFN_MODULE_NAME = "Type I and II IFN response"

# A renamed built-in module must not invalidate saved heatmap, correlation or
# concordance selections.  These aliases only resolve user-provided feature
# names; new outputs and manifests always contain the more precise names.
LEGACY_BUILTIN_PATHWAY_NAMES = {
    "Curated PPARA-target module": PPARA_TARGET_MODULE_NAME,
    "FA import": CARNITINE_SHUTTLE_MODULE_NAME,
    "Peroxisomal FAO": PEROXISOMAL_LIPID_METABOLISM_MODULE_NAME,
    "FA transport": INTRACELLULAR_FA_TRANSPORT_MODULE_NAME,
    "Inflammatory response": ACUTE_MYELOID_INFLAMMATION_MODULE_NAME,
    "IFN response": TYPE_I_II_IFN_MODULE_NAME,
}
DEFAULT_CONCORDANCE_TF_FEATURES = (
    "PPARA activity", "HNF4A activity", "ESRRA activity",
    "RELA activity", "STAT1 activity", "IRF1 activity",
)
DEFAULT_CONCORDANCE_PATHWAY_FEATURES = (
    f"{CARNITINE_SHUTTLE_MODULE_NAME} score", "Mitochondrial beta-oxidation score",
    f"{PEROXISOMAL_LIPID_METABOLISM_MODULE_NAME} score", f"{ACUTE_MYELOID_INFLAMMATION_MODULE_NAME} score",
    "TNF-NFkB response score", f"{TYPE_I_II_IFN_MODULE_NAME} score",
)


BUILTIN_PATHWAYS = {
    # This is a score_genes expression module.  It is deliberately named
    # differently from ``PPARA activity`` (network-backed ULM) so figures do
    # not imply that the two metrics are interchangeable.
    PPARA_TARGET_MODULE_NAME: (
        "CPT1A", "CPT2", "SLC25A20", "ACADM", "ACADVL", "ACADS", "HADHA",
        "HADHB", "ECHS1", "ETFDH", "ACAA2", "ACOX1", "EHHADH", "ACAA1",
        "HSD17B4", "ABCD3", "SCP2", "FABP1", "FABP2", "SLC27A2", "HMGCS2",
        "HMGCL", "BDH1",
    ),
    CARNITINE_SHUTTLE_MODULE_NAME: ("CPT1A", "CPT2", "SLC25A20"),
    "Mitochondrial beta-oxidation": (
        "ACADM", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "ETFDH", "ACAA2",
    ),
    PEROXISOMAL_LIPID_METABOLISM_MODULE_NAME: ("ACOX1", "EHHADH", "ACAA1", "HSD17B4", "ABCD3", "SCP2"),
    INTRACELLULAR_FA_TRANSPORT_MODULE_NAME: ("FABP1", "FABP2", "SLC27A2"),
    "Ketogenesis": ("HMGCS2", "HMGCL", "BDH1"),
    ACUTE_MYELOID_INFLAMMATION_MODULE_NAME: INFLAMMATORY_GENES,
    "TNF-NFkB response": TNF_NFKB_RESPONSE_GENES,
    TYPE_I_II_IFN_MODULE_NAME: ("STAT1", "IRF1", "IRF7", "CXCL10", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1"),
}

# This focused panel is opt-in.  It is intentionally not an automatic cell
# annotation and the three short epithelial signatures remain hypothesis
# scores; the standard Hallmark terms below are loaded from the frozen local
# registry when this focus is selected.
IBD_EPITHELIAL_QUICK_PATHWAYS = {
    "Mature absorptive/enterocyte differentiation": (
        "ALPI", "KRT20", "FABP1", "FABP2", "APOA4", "SI", "DPP4", "VIL1", "CA1",
    ),
    "Stem/TA regenerative state": (
        "LGR5", "OLFM4", "SMOC2", "ASCL2", "SOX9", "EPCAM", "MKI67", "TOP2A",
    ),
    "Cell-cycle proliferation": (
        "MKI67", "TOP2A", "CDK1", "CCNB1", "UBE2C", "TYMS", "HMGB2", "TUBA1B",
    ),
}
IBD_EPITHELIAL_MANAGED_TERMS = (
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_HYPOXIA",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    "HALLMARK_BILE_ACID_METABOLISM",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
)


def _safe_name(value, fallback="feature"):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return text or fallback


def _as_bool(value, default=False):
    if value is None:
        return bool(default)
    if isinstance(value, str):
        return value.strip().lower() not in {"", "0", "false", "no", "off"}
    return bool(value)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload):
    """Stable content identifier for built-in or inline non-file resources."""
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _bh(pvalues):
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(len(values), np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted
    order = np.where(valid)[0][np.argsort(values[valid])]
    running = 1.0
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        running = min(running, float(values[index]) * len(order) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _spearman_effect_size(rho):
    """Use explicit, conservative labels so tiny cell-level effects are not overread."""
    if not np.isfinite(rho):
        return "not estimable"
    magnitude = abs(float(rho))
    if magnitude < 0.10:
        return "negligible"
    if magnitude < 0.30:
        return "weak"
    if magnitude < 0.50:
        return "moderate"
    if magnitude < 0.70:
        return "strong"
    return "very strong"


def _wrapped_label(value, width=18):
    """Wrap long biological labels at semantic separators before word wrapping."""
    text = str(value or "")
    text = text.replace("/", "/\n")
    lines = []
    for part in text.splitlines() or [text]:
        lines.extend(textwrap.wrap(part, width=width, break_long_words=False) or [part])
    return "\n".join(lines)


def _pathway_display_label(feature, width=18):
    """Return a compact, wrapped label for pathway-score axes.

    Managed Hallmark terms are stored with machine-readable underscores.  Those
    identifiers must remain intact in tables, but plotting them verbatim makes
    wide concordance matrices unreadable.  This helper is display-only.
    """
    name = str(feature or "").removesuffix(" score")
    if name.startswith("HALLMARK_"):
        name = "Hallmark: " + name.removeprefix("HALLMARK_").replace("_", " ")
    else:
        name = name.replace("_", " ")
    replacements = {
        "TNFA": "TNFα",
        "NFKB": "NF-κB",
        "INTERFERON GAMMA": "IFNγ",
        "INTERFERON ALPHA": "IFNα",
        "REACTIVE OXYGEN SPECIES PATHWAY": "ROS pathway",
        "UNFOLDED PROTEIN RESPONSE": "unfolded-protein response",
        "OXIDATIVE PHOSPHORYLATION": "oxidative phosphorylation",
        "FATTY ACID METABOLISM": "fatty-acid metabolism",
        "CHOLESTEROL HOMEOSTASIS": "cholesterol homeostasis",
        "BILE ACID METABOLISM": "bile-acid metabolism",
        "INFLAMMATORY RESPONSE": "inflammatory response",
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    return _wrapped_label(name, width=width)


def _parse_gene_list(value):
    if isinstance(value, (list, tuple, set)):
        raw = value
    else:
        raw = re.split(r"[,;\n]+", str(value or ""))
    return tuple(dict.fromkeys(
        item.strip().upper() for item in raw if str(item).strip()
    ))


def _parse_custom_gene_sets(value):
    """Parse ``Name: GENE1,GENE2`` lines without interpreting expressions."""
    gene_sets = {}
    for line in str(value or "").splitlines():
        line = line.strip()
        if not line:
            continue
        if ":" not in line:
            raise ValueError("自定义基因集每行必须为“名称: GENE1,GENE2”。")
        name, genes = line.split(":", 1)
        name = re.sub(r"\s+", " ", name).strip()
        if not name or len(name) > 100:
            raise ValueError("自定义基因集名称不能为空且最长 100 个字符。")
        parsed = _parse_gene_list(genes)
        if not parsed:
            raise ValueError(f"自定义基因集 '{name}' 没有有效基因。")
        if len(parsed) > 2000:
            raise ValueError(f"自定义基因集 '{name}' 超过 2000 个基因。")
        if name in gene_sets:
            raise ValueError(f"自定义基因集名称重复: '{name}'。")
        gene_sets[name] = parsed
    return gene_sets


def _parse_feature_selection(value):
    """Parse an ordered, bounded list of score names for display only."""
    selected = []
    for item in re.split(r"[,;\n]+", str(value or "")):
        name = re.sub(r"\s+", " ", item).strip()
        if name and name not in selected:
            selected.append(name)
    if len(selected) > 40:
        raise ValueError("热图展示特征一次最多选择 40 个。")
    return tuple(selected)


def _resolve_project_file(project_dir, value, label):
    """Resolve an optional project-local resource and reject symlinks/traversal."""
    supplied = str(value or "").strip()
    if not supplied:
        return None
    root = Path(project_dir).resolve()
    candidate = Path(supplied)
    if not candidate.is_absolute():
        candidate = root / candidate
    if candidate.is_symlink():
        raise ValueError(f"{label}不能是符号链接。")
    resolved = candidate.resolve()
    if root not in resolved.parents or not resolved.is_file():
        raise ValueError(f"{label}必须是当前项目目录内的现有普通文件。")
    return resolved


def _read_tf_network(path):
    """Read a local CollecTRI/DoRothEA-like network with explicit columns."""
    if path is None:
        return pd.DataFrame(columns=["tf", "target", "weight"])
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    network = pd.read_csv(path, sep=separator)
    normalized = {str(column).strip().lower(): column for column in network.columns}
    source_column = next((normalized[key] for key in ("source", "tf", "regulator") if key in normalized), None)
    target_column = next((normalized[key] for key in ("target", "gene", "target_gene") if key in normalized), None)
    weight_column = next((normalized[key] for key in ("weight", "mor", "confidence") if key in normalized), None)
    if source_column is None or target_column is None:
        raise ValueError("TF 网络需要 source/tf 和 target/gene 列。")
    weights = pd.to_numeric(network[weight_column], errors="coerce") if weight_column else pd.Series(1.0, index=network.index)
    result = pd.DataFrame({
        "tf": network[source_column].astype(str).str.strip().str.upper(),
        "target": network[target_column].astype(str).str.strip().str.upper(),
        "weight": weights.fillna(1.0),
    })
    result = result[(result["tf"] != "") & (result["target"] != "")].copy()
    result["weight"] = result["weight"].clip(lower=-10, upper=10)
    result = result[result["weight"] != 0]
    if result.empty:
        raise ValueError("TF 网络没有可用的 source-target 记录。")
    return result.drop_duplicates(["tf", "target"], keep="last")


def _read_resource_metadata(path):
    """Read optional, non-authoritative metadata alongside a managed resource."""
    metadata_path = Path(path).with_suffix(".metadata.json")
    if metadata_path.is_symlink() or not metadata_path.is_file():
        return {}
    try:
        payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _resolve_tf_network(project_dir, params):
    """Select a controlled CollecTRI snapshot or a project-local regulon.

    A supplied project file remains backwards-compatible for callers that
    predate ``tf_network_source``.  The managed option itself never accepts a
    browser-provided path, so a user cannot redirect analysis to an arbitrary
    server file.
    """
    supplied_path = str(params.get("tf_network_path", "") or "").strip()
    requested = str(params.get("tf_network_source", "") or "").strip().lower()
    if not requested:
        requested = "project_local" if supplied_path else "managed_collectri"
    aliases = {
        "managed_collectri": "managed_collectri",
        "project_local": "project_local",
        "project_file": "project_local",
        "not_configured": "not_configured",
        "disabled": "not_configured",
    }
    if requested not in aliases:
        raise ValueError("TF 网络来源只能为 managed_collectri、project_local 或 not_configured。")
    source = aliases[requested]
    if source == "not_configured":
        return None, source, {}
    if source == "project_local":
        if not supplied_path:
            raise ValueError("选择项目内 TF 网络时，必须填写当前项目内网络文件路径。")
        return _resolve_project_file(project_dir, supplied_path, "TF 网络文件"), source, {}

    # Importing Config initialises the application-controlled Numba and
    # Matplotlib cache directories before decoupler is imported below.
    from config import Config

    candidate = Path(Config.functional_state_resource_dir()) / "collectri_human.tsv"
    if candidate.is_symlink():
        raise ValueError("平台 CollecTRI 资源不能是符号链接。")
    if not candidate.is_file():
        return None, source, {}
    metadata = _read_resource_metadata(candidate)
    expected_sha = str(metadata.get("sha256", "") or "").strip().lower()
    if expected_sha and expected_sha != _sha256(candidate):
        raise ValueError("平台 CollecTRI 快照校验失败；请由管理员重新下载并生成元数据。")
    return candidate.resolve(), source, metadata


def _read_selected_gmt(path, selected_terms):
    """Read an explicit, bounded set of project-local GMT pathways."""
    if path is None:
        return {}
    selected = [item.strip() for item in re.split(r"[,;\n]+", str(selected_terms or "")) if item.strip()]
    if not selected:
        raise ValueError("填写本地 GMT 后，必须明确选择要评分的 pathway term。")
    if len(selected) > 60:
        raise ValueError("一次最多选择 60 个 GMT pathway term。")
    targets = set(selected)
    found = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in targets:
                continue
            genes = _parse_gene_list(fields[2:])
            if genes:
                found[fields[0]] = genes
    missing = [term for term in selected if term not in found]
    if missing:
        raise ValueError("本地 GMT 未找到所选 term: " + "、".join(missing[:8]))
    return found


def _pathway_gene_set_overlap_qc(pathway_sets, gene_lookup):
    """Audit pathway-set overlap in the actual expression universe.

    Two differently named signatures can collapse to the same detected genes
    after an input-specific gene intersection.  Scoring both would generate
    numerically identical results and falsely look like independent evidence,
    so this is a hard validation failure rather than a post-hoc warning.
    """
    resolved = {}
    for name, genes in pathway_sets.items():
        requested = frozenset(str(gene).strip().upper() for gene in genes if str(gene).strip())
        detected = frozenset(gene_lookup[gene] for gene in requested if gene in gene_lookup)
        resolved[name] = {"requested": requested, "detected": detected}

    rows, duplicates = [], []
    for left, right in combinations(pathway_sets, 2):
        left_sets, right_sets = resolved[left], resolved[right]
        requested_shared = left_sets["requested"] & right_sets["requested"]
        detected_shared = left_sets["detected"] & right_sets["detected"]
        requested_union = left_sets["requested"] | right_sets["requested"]
        detected_union = left_sets["detected"] | right_sets["detected"]
        identical_requested = bool(left_sets["requested"]) and left_sets["requested"] == right_sets["requested"]
        identical_detected = bool(left_sets["detected"]) and left_sets["detected"] == right_sets["detected"]
        if identical_requested:
            status = "blocked_identical_requested_genes"
        elif identical_detected:
            status = "blocked_identical_detected_genes"
        elif detected_union and len(detected_shared) / len(detected_union) >= 0.80:
            status = "high_detected_overlap_review"
        elif detected_union:
            status = "distinct_detected_genes"
        else:
            status = "no_shared_detected_genes"
        rows.append({
            "gene_set_a": left,
            "gene_set_b": right,
            "n_requested_a": len(left_sets["requested"]),
            "n_requested_b": len(right_sets["requested"]),
            "n_shared_requested": len(requested_shared),
            "requested_jaccard": len(requested_shared) / len(requested_union) if requested_union else np.nan,
            "n_detected_a": len(left_sets["detected"]),
            "n_detected_b": len(right_sets["detected"]),
            "n_shared_detected": len(detected_shared),
            "detected_jaccard": len(detected_shared) / len(detected_union) if detected_union else np.nan,
            "status": status,
        })
        if identical_requested or identical_detected:
            duplicates.append(f"{left} ↔ {right}")
    if duplicates:
        raise ValueError(
            "基因集重复 QC 阻止本次评分：以下不同名称在请求或实际检测基因层面完全相同："
            + "；".join(duplicates)
            + "。请保留一个基因集或改用彼此独立的定义。"
        )
    return pd.DataFrame(rows)


def _validate_obs_column(adata, key, label, require_multiple=False):
    key = str(key or "").strip()
    if not key or key not in adata.obs.columns:
        raise ValueError(f"{label}列 '{key}' 不在 adata.obs 中。")
    values = adata.obs[key]
    missing = values.isna() | values.astype(str).str.strip().eq("")
    if bool(missing.any()):
        raise ValueError(f"{label}列 '{key}' 有 {int(missing.sum())} 个缺失值。")
    if require_multiple and values.astype(str).nunique() < 2:
        raise ValueError(f"{label}列 '{key}' 至少需要两个取值。")
    return key


def _resolve_expression(adata, requested_layer, source_mode="auto"):
    """Choose a safe scoring matrix and record the decision without mutating X.

    Raw counts and residual/scaled matrices are never scored directly.  They
    are transformed into ``layers['_functional_log1p']`` on the output AnnData,
    preserving the input ``X`` and making the scoring scale auditable.
    """
    source_mode = str(source_mode or "auto").strip().lower()
    if source_mode not in {"auto", "x", "layers[counts]", "custom_layer"}:
        raise ValueError("表达来源只能为 auto、X、layers[counts] 或 custom_layer。")
    if source_mode == "x":
        layer = "X"
    elif source_mode == "layers[counts]":
        layer = "counts"
    else:
        layer = str(requested_layer or "").strip()
        if source_mode == "custom_layer" and layer.lower() in {"", "x"}:
            raise ValueError("选择 custom_layer 时必须填写 adata.layers 中的 log-normalized 表达层名称。")
        if source_mode == "auto":
            layer = "X"

    if layer.lower() in {"", "x"}:
        matrix = adata.X
        source = "X"
    else:
        if layer not in adata.layers:
            raise ValueError(f"表达层 '{layer}' 不在 adata.layers 中。")
        matrix = adata.layers[layer]
        source = f"layers[{layer}]"

    normalization = adata.uns.get("normalization", {}) or {}
    values = matrix.data if hasattr(matrix, "data") and hasattr(matrix, "tocsr") else np.asarray(matrix).reshape(-1)
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    is_negative = bool(finite.size and float(finite.min()) < -1e-8)
    is_raw_counts = bool(_matrix_is_raw_counts(matrix))
    qc = {
        "requested_source": source_mode,
        "requested_layer": layer,
        "selected_input": source,
        "contains_negative_values": is_negative,
        "is_integer_count_matrix": is_raw_counts,
        "has_log1p_metadata": bool(adata.uns.get("log1p")),
        "temporary_scoring_layer": None,
        "decision": "",
        "warning": "",
    }

    def rebuild_from_counts(counts, reason):
        if counts is None or not _matrix_is_raw_counts(counts):
            raise ValueError(
                "功能状态评分需要 log-normalized expression；当前来源含 residual/负值，且没有有效原始 counts 层可重建。"
            )
        transformed = log1p_expression_from_counts(counts)
        adata.layers["_functional_log1p"] = transformed
        qc.update({
            "temporary_scoring_layer": "layers[_functional_log1p]",
            "decision": reason,
        })
        return transformed, "layers[_functional_log1p]"

    # Pearson residuals are useful for embeddings, not for gene/programme/ULM
    # scores.  Use protected counts rather than silently accepting them.
    if normalization.get("x_contains") == "pearson_residuals" or is_negative:
        expression, selected = rebuild_from_counts(
            adata.layers.get("counts") if hasattr(adata, "layers") else None,
            "reconstructed_from_counts_due_to_negative_or_pearson_residuals",
        )
        return expression, selected, qc
    if is_raw_counts:
        expression, selected = rebuild_from_counts(matrix, "reconstructed_from_raw_counts")
        return expression, selected, qc

    qc["decision"] = "accepted_nonnegative_continuous_matrix"
    if not qc["has_log1p_metadata"] and "log" not in source.lower():
        qc["warning"] = (
            "当前矩阵为非负非整数连续值，按 log-normalized expression 使用；"
            "请在上游流程或输入说明中确认其不是 scaled matrix。"
        )
    return matrix, source, qc


def _dense_columns(matrix, indices):
    selected = matrix[:, indices]
    if hasattr(selected, "toarray"):
        selected = selected.toarray()
    return np.asarray(selected, dtype=float)


def _score_gene_set(expression, var_names, detected_genes, *, score_name, requested_method):
    """Use the explicitly selected pathway-score method with an honest fallback."""
    if str(requested_method or "").strip() != "scanpy_score_genes":
        raise ValueError("当前仅实现 scanpy_score_genes；不能把未实现的方法写入结果。")
    from anndata import AnnData
    import scanpy as sc

    work = AnnData(
        X=expression,
        obs=pd.DataFrame(index=np.arange(expression.shape[0]).astype(str)),
        var=pd.DataFrame(index=pd.Index(var_names)),
    )
    try:
        sc.tl.score_genes(
            work, gene_list=list(detected_genes), score_name="score",
            ctrl_size=min(50, max(1, len(detected_genes))), n_bins=25,
            random_state=0, use_raw=False,
        )
        values = pd.to_numeric(work.obs["score"], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(values).any():
            return values, "scanpy_score_genes"
    except Exception:
        pass
    positions = [int(pd.Index(var_names).get_loc(gene)) for gene in detected_genes]
    return np.nanmean(_dense_columns(expression, positions), axis=1), "mean_expression_fallback"


def _ulm_regulon_scores(expression, var_names, network, requested_tfs, *, network_name):
    """Infer signed TF activity with decoupler's ULM on log-normalised RNA.

    The returned values are ULM slope t-statistics, not TF-gene expression and
    not an ad-hoc target mean.  The five-target eligibility floor is fixed so
    the four coverage tiers have the same meaning across analysis runs.
    """
    from anndata import AnnData
    from config import Config

    Config.configure_runtime_tmpdir()
    try:
        import decoupler as dc
    except ImportError as exc:
        raise ImportError("TF 活性需要 decoupler；请安装 requirements.txt 中的固定版本。") from exc

    gene_lookup = {str(gene).upper(): str(gene) for gene in var_names}
    min_targets = 5
    records, scores, definitions = [], {}, {}
    requested = tuple(dict.fromkeys(tf.upper() for tf in requested_tfs))
    coverage = {}
    for tf in requested:
        full = network.loc[network["tf"].eq(tf), ["target", "weight"]].copy()
        full_targets = tuple(full["target"].astype(str).tolist())
        detected = [gene_lookup[gene] for gene in full_targets if gene in gene_lookup]
        missing = [gene for gene in full_targets if gene not in gene_lookup]
        coverage[tf] = (full_targets, detected, missing)

    eligible = [tf for tf, (_, detected, _) in coverage.items() if len(detected) >= int(min_targets)]
    ulm_scores = pd.DataFrame(index=np.arange(expression.shape[0]).astype(str))
    if eligible:
        net = network.loc[network["tf"].isin(eligible), ["tf", "target", "weight"]].copy()
        net["target"] = net["target"].map(gene_lookup)
        net = net.dropna(subset=["target"]).rename(columns={"tf": "source"})
        work = AnnData(
            X=expression,
            obs=pd.DataFrame(index=np.arange(expression.shape[0]).astype(str)),
            var=pd.DataFrame(index=pd.Index(var_names)),
        )
        try:
            dc.mt.ulm(data=work, net=net, tmin=int(min_targets), raw=False, verbose=False)
            raw_scores = work.obsm.get("score_ulm")
            if raw_scores is not None:
                ulm_scores = pd.DataFrame(raw_scores, index=work.obs_names)
        except Exception as exc:
            raise RuntimeError(f"decoupler ULM 运行失败：{exc}") from exc

    for tf in requested:
        full_targets, detected, missing = coverage[tf]
        status, column = "not_in_network", ""
        if not full_targets:
            status = "not_in_network"
        elif len(detected) < min_targets:
            status = "below_minimum_coverage"
            coverage_level = "below_5_not_scored"
        elif tf not in ulm_scores.columns:
            status = "ulm_no_score"
            coverage_level = "eligible_but_no_score"
        else:
            values = pd.to_numeric(ulm_scores[tf], errors="coerce").to_numpy(dtype=float)
            if not np.isfinite(values).any():
                status = "ulm_no_score"
                coverage_level = "eligible_but_no_score"
            else:
                column = "fs_tf_" + _safe_name(tf) + "_activity"
                scores[tf] = values
                definitions[f"{tf} activity"] = tuple(detected)
                if len(detected) < 10:
                    status = "low_coverage"
                    coverage_level = "low_5_9"
                elif len(detected) < 20:
                    status = "acceptable"
                    coverage_level = "acceptable_10_19"
                else:
                    status = "good"
                    coverage_level = "good_20_plus"
        if not full_targets:
            coverage_level = "not_in_network"
        records.append({
            "tf": tf, "network": network_name, "n_network_targets": len(full_targets),
            "n_detected_targets": len(detected),
            "coverage": len(detected) / max(len(full_targets), 1),
            "network_targets": len(full_targets), "detected_targets": len(detected),
            "coverage_fraction": len(detected) / max(len(full_targets), 1),
            "minimum_required_targets": min_targets,
            "coverage_level": coverage_level,
            "missing_targets": ";".join(missing), "status": status,
            "activity_status": status,
            "display_label": f"{tf}*" if status == "low_coverage" else tf,
            "score_column": column, "method": "decoupler_ulm",
            "score_interpretation": "ULM slope t-statistic; positive=activated, negative=repressed",
        })
    return pd.DataFrame(records), scores, definitions


def _choose_comparison(conditions, reference, comparison):
    labels = list(dict.fromkeys(str(value) for value in conditions))
    reference, comparison = str(reference or "").strip(), str(comparison or "").strip()
    if reference or comparison:
        if reference not in labels or comparison not in labels or reference == comparison:
            raise ValueError("参考条件和比较条件必须是 condition 列中两个不同的实际取值。")
        return reference, comparison
    if len(labels) != 2:
        raise ValueError("condition 有两个以上取值时，必须明确填写参考条件和比较条件。")
    return tuple(sorted(labels))


def _sample_feature_table(obs, feature_values, feature_types, *, sample_key, condition_key, celltype_key, min_cells):
    base = pd.DataFrame({
        "sample_id": obs[sample_key].astype(str).to_numpy(),
        "condition": obs[condition_key].astype(str).to_numpy(),
        "celltype": obs[celltype_key].astype(str).to_numpy(),
    })
    rows = []
    for feature, values in feature_values.items():
        frame = base.copy()
        frame["score"] = np.asarray(values, dtype=float)
        grouped = frame.groupby(["sample_id", "condition", "celltype"], observed=True)["score"]
        table = grouped.agg(["mean", "count"]).reset_index()
        table = table.rename(columns={"mean": "score", "count": "n_cells"})
        table["feature"] = feature
        table["feature_type"] = feature_types[feature]
        table["eligible_for_statistics"] = table["n_cells"].ge(int(min_cells))
        rows.append(table)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["sample_id", "condition", "celltype", "score", "n_cells", "feature", "feature_type", "eligible_for_statistics"]
    )


def _sample_size_evidence_tier(n_reference, n_comparison):
    """Classify biological replication without mistaking cells for samples."""
    n_min = min(int(n_reference), int(n_comparison))
    if n_min < 2:
        return "not_runnable"
    if n_min == 2:
        return "exploratory_n2"
    if n_min < 5:
        return "standard_n3_4"
    return "preferred_n5_plus"


def _feature_statistics(sample_table, reference, comparison, min_samples, method):
    from scipy.stats import mannwhitneyu, ttest_ind

    rows = []
    groups = sample_table.groupby(["feature_type", "feature", "celltype"], observed=True)
    for (feature_type, feature, celltype), frame in groups:
        eligible = frame[frame["eligible_for_statistics"]].copy()
        ref_values = eligible.loc[eligible["condition"].eq(reference), "score"].to_numpy(dtype=float)
        cmp_values = eligible.loc[eligible["condition"].eq(comparison), "score"].to_numpy(dtype=float)
        ref_values = ref_values[np.isfinite(ref_values)]
        cmp_values = cmp_values[np.isfinite(cmp_values)]
        evidence_tier = _sample_size_evidence_tier(len(ref_values), len(cmp_values))
        row = {
            "feature_type": feature_type, "feature": feature, "celltype": celltype,
            "reference_condition": reference, "comparison_condition": comparison,
            "n_samples_reference": len(ref_values), "n_samples_comparison": len(cmp_values),
            "mean_reference": float(np.mean(ref_values)) if len(ref_values) else np.nan,
            "mean_comparison": float(np.mean(cmp_values)) if len(cmp_values) else np.nan,
            "mean_difference": (float(np.mean(cmp_values) - np.mean(ref_values))
                                if len(ref_values) and len(cmp_values) else np.nan),
            "test": method, "statistic": np.nan, "pvalue": np.nan,
            "evidence_tier": evidence_tier,
            "status": "insufficient_sample_units",
        }
        if len(ref_values) >= int(min_samples) and len(cmp_values) >= int(min_samples):
            try:
                if method == "mann_whitney":
                    result = mannwhitneyu(cmp_values, ref_values, alternative="two-sided")
                    row.update({"statistic": float(result.statistic), "pvalue": float(result.pvalue)})
                else:
                    result = ttest_ind(cmp_values, ref_values, equal_var=False, nan_policy="omit")
                    row.update({"statistic": float(result.statistic), "pvalue": float(result.pvalue)})
                row["status"] = (
                    "sample_level_exploratory" if evidence_tier == "exploratory_n2"
                    else "sample_level_test"
                )
            except Exception:
                row["status"] = "test_unavailable"
        rows.append(row)
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr_bh"] = _bh(result["pvalue"])
    return result


def _resolve_feature_name(value, feature_values):
    text = str(value or "").strip()
    suffix = " score" if text.endswith(" score") else ""
    bare_name = text[:-len(suffix)] if suffix else text
    text = LEGACY_BUILTIN_PATHWAY_NAMES.get(bare_name, bare_name) + suffix
    if text in feature_values:
        return text
    normalized = re.sub(r"[^a-z0-9]+", "", text.lower())
    for feature in feature_values:
        if re.sub(r"[^a-z0-9]+", "", feature.lower()) == normalized:
            return feature
    return None


def _select_heatmap_features(value, feature_values, feature_types):
    """Return display features in user order without altering computed scores."""
    available = [
        feature for feature in feature_values
        if feature_types.get(feature) in {"pathway_score", "tf_activity"}
    ]
    requested = _parse_feature_selection(value)
    if not requested:
        return available, [], "all_available_scores"

    selected, missing = [], []
    for requested_name in requested:
        resolved = _resolve_feature_name(requested_name, feature_values)
        if resolved in available:
            if resolved not in selected:
                selected.append(resolved)
        else:
            missing.append(requested_name)
    return selected, missing, "user_ordered_selection"


def _tf_expression_activity_table(gene_statistics, tf_statistics, tf_coverage):
    """Join TF expression and ULM activity contrasts without conflating them."""
    columns = [
        "tf", "celltype", "expression_mean_difference", "expression_fdr_bh",
        "expression_evidence_tier", "activity_mean_difference", "activity_fdr_bh",
        "activity_evidence_tier", "activity_status", "coverage_level",
        "detected_targets", "network_targets",
    ]
    if tf_coverage.empty:
        return pd.DataFrame(columns=columns)

    requested_tfs = set(tf_coverage["tf"].astype(str).str.upper())
    expression = gene_statistics.copy()
    activity = tf_statistics.copy()
    expression["tf"] = expression["feature"].str.replace(r" expression$", "", regex=True)
    activity["tf"] = activity["feature"].str.replace(r" activity$", "", regex=True)
    expression["tf"] = expression["tf"].str.upper()
    activity["tf"] = activity["tf"].str.upper()
    expression = expression[expression["tf"].isin(requested_tfs)]
    activity = activity[activity["tf"].isin(requested_tfs)]
    expression = expression[["tf", "celltype", "mean_difference", "fdr_bh", "evidence_tier"]].rename(columns={
        "mean_difference": "expression_mean_difference",
        "fdr_bh": "expression_fdr_bh",
        "evidence_tier": "expression_evidence_tier",
    })
    activity = activity[["tf", "celltype", "mean_difference", "fdr_bh", "evidence_tier"]].rename(columns={
        "mean_difference": "activity_mean_difference",
        "fdr_bh": "activity_fdr_bh",
        "evidence_tier": "activity_evidence_tier",
    })
    merged = expression.merge(activity, on=["tf", "celltype"], how="outer")
    coverage = tf_coverage[[
        "tf", "activity_status", "coverage_level", "n_detected_targets", "n_network_targets",
    ]].copy().rename(columns={
        "n_detected_targets": "detected_targets", "n_network_targets": "network_targets",
    })
    coverage["tf"] = coverage["tf"].astype(str).str.upper()
    return merged.merge(coverage, on="tf", how="left").reindex(columns=columns).sort_values(
        ["tf", "celltype"], kind="stable", na_position="last"
    ).reset_index(drop=True)


def _parse_correlation_levels(value):
    allowed = {
        "cell_level_descriptive", "sample_x_celltype", "within_celltype_sample",
        "within_celltype_cell_level_descriptive",
    }
    raw = [item.strip() for item in re.split(r"[,;\n]+", str(value or "")) if item.strip()]
    selected = raw or [
        "cell_level_descriptive", "sample_x_celltype", "within_celltype_sample",
        "within_celltype_cell_level_descriptive",
    ]
    unknown = sorted(set(selected).difference(allowed))
    if unknown:
        raise ValueError("未知相关性层级: " + "、".join(unknown))
    return tuple(dict.fromkeys(selected))


def _requested_correlation_feature(params, axis):
    selected = str(params.get(f"correlation_{axis}", "") or "").strip()
    if selected == "__custom__":
        selected = str(params.get(f"correlation_{axis}_custom", "") or "").strip()
    if not selected:
        return "PPARA activity" if axis == "x" else f"{ACUTE_MYELOID_INFLAMMATION_MODULE_NAME} score"
    return selected


def _correlations(
    obs, feature_values, feature_definitions, sample_table, *, x_name, y_name,
    celltype_key, condition_key, conditions, min_samples, levels,
):
    from scipy.stats import spearmanr

    rows = []
    x_feature = _resolve_feature_name(x_name, feature_values)
    y_feature = _resolve_feature_name(y_name, feature_values)
    if x_feature is None or y_feature is None:
        missing = x_name if x_feature is None else y_name
        return pd.DataFrame([{
            "analysis_level": "not_available", "condition": "all_conditions",
            "x_feature": x_name, "y_feature": y_name,
            "n_units": 0, "n_unique_samples": 0, "repeated_donor_measures": False,
            "rho": np.nan, "effect_size_interpretation": "not estimable",
            "pvalue": np.nan, "fdr_bh": np.nan, "fdr_interpretation": "not_reported",
            "status": f"feature_not_available: {missing}",
        }]), pd.DataFrame()

    def add_row(level, x, y, status, *, condition="all_conditions", celltype="All", n_unique_samples=np.nan,
                repeated_donor_measures=False, fdr_interpretation="not_reported"):
        finite = np.isfinite(x) & np.isfinite(y)
        rho, pvalue = np.nan, np.nan
        if int(finite.sum()) >= 3:
            result = spearmanr(x[finite], y[finite])
            rho, pvalue = float(result.statistic), float(result.pvalue)
        # Cell-level rows are visual/descriptive effect-size checks.  Cells
        # are not independent biological replicates, so retaining Spearman's
        # mathematical P value in an exported CSV invites pseudoreplication.
        if "cell_level_descriptive" in level:
            pvalue = np.nan
        rows.append({
            "analysis_level": level, "condition": condition, "celltype": celltype, "x_feature": x_feature,
            "y_feature": y_feature, "n_units": int(finite.sum()), "rho": rho,
            "n_unique_samples": n_unique_samples,
            "repeated_donor_measures": bool(repeated_donor_measures),
            "effect_size_interpretation": _spearman_effect_size(rho),
            "pvalue": pvalue, "fdr_bh": np.nan,
            "fdr_interpretation": fdr_interpretation, "status": status,
        })

    # Cells are retained as a visual/descriptive check only; their p value is
    # not an independent-replicate inference.
    if "cell_level_descriptive" in levels:
        add_row(
            "cell_level_descriptive", np.asarray(feature_values[x_feature]),
            np.asarray(feature_values[y_feature]), "descriptive_only",
            n_unique_samples=int(sample_table["sample_id"].astype(str).nunique()),
        )

    # Use only sample × cell-type units meeting the cell-count threshold; the
    # remaining rows stay in the source table for transparent QC but must not
    # alter an aggregate association.
    eligible_sample_table = sample_table[sample_table["eligible_for_statistics"]].copy()
    pivot = eligible_sample_table.pivot_table(
        index=["sample_id", "condition", "celltype"], columns="feature", values="score", aggfunc="first",
    )
    if x_feature in pivot and y_feature in pivot:
        required_units = max(3, 2 * int(min_samples))
        if "sample_x_celltype" in levels:
            status = (
                "exploratory_repeated_measures"
                if len(pivot) >= required_units else "insufficient_sample_celltype_units"
            )
            add_row(
                "sample_x_celltype", pivot[x_feature].to_numpy(float), pivot[y_feature].to_numpy(float), status,
                n_unique_samples=int(pivot.index.get_level_values("sample_id").nunique()),
                repeated_donor_measures=True,
            )
            selected_conditions = tuple(dict.fromkeys(str(value) for value in (conditions or ())))
            if not selected_conditions:
                selected_conditions = tuple(sorted(pivot.index.get_level_values("condition").astype(str).unique()))
            for condition in selected_conditions:
                table = pivot[pivot.index.get_level_values("condition").astype(str) == condition]
                n_unique_samples = int(table.index.get_level_values("sample_id").nunique()) if len(table) else 0
                status = (
                    "exploratory_repeated_measures_condition_stratified"
                    if len(table) >= required_units and n_unique_samples >= 2
                    else "insufficient_condition_stratified_sample_celltype_units"
                )
                add_row(
                    "sample_x_celltype_by_condition", table[x_feature].to_numpy(float),
                    table[y_feature].to_numpy(float), status, condition=condition,
                    n_unique_samples=n_unique_samples, repeated_donor_measures=True,
                )
        if "within_celltype_sample" in levels:
            for celltype, table in pivot.groupby(level="celltype", observed=True):
                status = (
                    "exploratory_unadjusted_within_celltype"
                    if len(table) >= required_units else "insufficient_sample_units"
                )
                add_row(
                    "within_celltype_sample", table[x_feature].to_numpy(float), table[y_feature].to_numpy(float), status,
                    celltype=str(celltype),
                    n_unique_samples=int(table.index.get_level_values("sample_id").nunique()),
                    # Within one cell type there is one point per donor, but
                    # this bivariate test is still unadjusted for condition.
                    repeated_donor_measures=False,
                    fdr_interpretation="exploratory_unadjusted",
                )
    if "within_celltype_cell_level_descriptive" in levels:
        celltype_values = obs[celltype_key].astype(str)
        condition_values = obs[condition_key].astype(str)
        selected_conditions = tuple(dict.fromkeys(str(value) for value in (conditions or ())))
        if not selected_conditions:
            selected_conditions = tuple(sorted(condition_values.unique()))
        for celltype in sorted(celltype_values.unique()):
            mask = celltype_values.eq(celltype).to_numpy()
            n_samples = int(
                sample_table.loc[sample_table["celltype"].astype(str).eq(celltype), "sample_id"]
                .astype(str).nunique()
            )
            add_row(
                "within_celltype_cell_level_descriptive",
                np.asarray(feature_values[x_feature])[mask], np.asarray(feature_values[y_feature])[mask],
                "descriptive_only", celltype=celltype, n_unique_samples=n_samples,
            )
            for condition in selected_conditions:
                condition_mask = (celltype_values.eq(celltype) & condition_values.eq(condition)).to_numpy()
                n_condition_samples = int(
                    sample_table.loc[
                        sample_table["celltype"].astype(str).eq(celltype)
                        & sample_table["condition"].astype(str).eq(condition), "sample_id"
                    ].astype(str).nunique()
                )
                add_row(
                    "within_celltype_cell_level_descriptive_by_condition",
                    np.asarray(feature_values[x_feature])[condition_mask],
                    np.asarray(feature_values[y_feature])[condition_mask],
                    "descriptive_only_condition_stratified", condition=condition, celltype=celltype,
                    n_unique_samples=n_condition_samples,
                )

    correlations = pd.DataFrame(rows)
    if not correlations.empty:
        # Only within-cell-type sample correlations have one point per donor.
        # They remain exploratory because condition is not adjusted, but an
        # exploratory BH value can help rank the panel without overstating it.
        exploratory = correlations["analysis_level"].eq("within_celltype_sample")
        correlations.loc[exploratory, "fdr_bh"] = _bh(correlations.loc[exploratory, "pvalue"])

    overlap = pd.DataFrame()
    genes_x = set(feature_definitions.get(x_feature, ()))
    genes_y = set(feature_definitions.get(y_feature, ()))
    if genes_x or genes_y:
        shared = sorted(genes_x.intersection(genes_y))
        union = genes_x.union(genes_y)
        overlap = pd.DataFrame([{
            "feature_a": x_feature, "feature_b": y_feature, "n_a": len(genes_x), "n_b": len(genes_y),
            "n_overlap": len(shared), "jaccard": len(shared) / len(union) if union else np.nan,
            "overlap_genes": ";".join(shared),
            "warning": "Shared genes can partly drive this correlation." if shared else "",
        }])
    return correlations, overlap


def _pair_umap_figure(adata, feature_values, feature_x, feature_y):
    import matplotlib.pyplot as plt

    if "X_umap" not in adata.obsm or feature_x not in feature_values or feature_y not in feature_values:
        return None
    coords = np.asarray(adata.obsm["X_umap"], dtype=float)[:, :2]
    values = np.concatenate([np.asarray(feature_values[feature_x]), np.asarray(feature_values[feature_y])])
    finite = values[np.isfinite(values)]
    if not len(finite):
        return None
    vmin, vmax = np.nanpercentile(finite, [1, 99])
    if np.isclose(vmin, vmax):
        vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite) + 1e-6)
    # A dedicated colourbar axis is essential here.  Attaching a shared
    # colourbar to the two plot axes lets Matplotlib draw it over the right
    # UMAP panel for wide UMAP coordinates or long vertical labels.
    fig = plt.figure(figsize=(11.8, 4.65), dpi=150)
    grid = fig.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 0.040), wspace=0.055)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    colorbar_axis = fig.add_subplot(grid[0, 2])
    for axis, name in zip(axes, (feature_x, feature_y)):
        artist = axis.scatter(coords[:, 0], coords[:, 1], c=np.asarray(feature_values[name]), s=4,
                              cmap=nature_continuous_cmap(), vmin=vmin, vmax=vmax,
                              alpha=0.72, linewidths=0, rasterized=True)
        axis.set_title(name, loc="left", fontsize=10, fontweight="semibold")
        axis.set_xticks([]); axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
    colorbar = fig.colorbar(artist, cax=colorbar_axis)
    colorbar.outline.set_visible(False)
    colorbar.set_label("Score\n(shared 1–99% range)", fontsize=8, labelpad=7)
    colorbar.ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=0.025, right=0.965, bottom=0.055, top=0.91)
    return fig


def _condition_umap_figure(adata, values, condition_key, reference, comparison, title):
    import matplotlib.pyplot as plt

    if "X_umap" not in adata.obsm:
        return None
    coords = np.asarray(adata.obsm["X_umap"], dtype=float)[:, :2]
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if not len(finite):
        return None
    vmin, vmax = np.nanpercentile(finite, [1, 99])
    if np.isclose(vmin, vmax):
        vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite) + 1e-6)
    labels = adata.obs[condition_key].astype(str).to_numpy()
    fig = plt.figure(figsize=(11.8, 4.65), dpi=150)
    grid = fig.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 0.040), wspace=0.055)
    axes = [fig.add_subplot(grid[0, 0])]
    axes.append(fig.add_subplot(grid[0, 1], sharex=axes[0], sharey=axes[0]))
    colorbar_axis = fig.add_subplot(grid[0, 2])
    for axis, condition in zip(axes, (reference, comparison)):
        mask = labels == condition
        artist = axis.scatter(coords[mask, 0], coords[mask, 1], c=values[mask], s=4,
                              cmap=nature_continuous_cmap(), vmin=vmin, vmax=vmax,
                              alpha=0.72, linewidths=0, rasterized=True)
        axis.set_title(f"{condition}: {title}", loc="left", fontsize=10, fontweight="semibold")
        axis.set_xticks([]); axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
    colorbar = fig.colorbar(artist, cax=colorbar_axis)
    colorbar.outline.set_visible(False)
    colorbar.set_label("Score\n(shared 1–99% range)", fontsize=8, labelpad=7)
    colorbar.ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=0.025, right=0.965, bottom=0.055, top=0.91)
    return fig


def _sample_fdr_label(statistics, feature, celltype):
    """Translate the already-computed sample-level BH-FDR into a plot label."""
    if statistics is None or statistics.empty:
        return "n/a"
    rows = statistics[
        statistics["feature"].eq(feature)
        & statistics["celltype"].astype(str).eq(str(celltype))
        & statistics["status"].isin(["sample_level_test", "sample_level_exploratory"])
    ]
    if rows.empty:
        return "n/a"
    if str(rows.iloc[0].get("evidence_tier", "")) == "exploratory_n2":
        return "n=2"
    fdr = float(rows.iloc[0].get("fdr_bh", np.nan))
    if not np.isfinite(fdr):
        return "n/a"
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    if fdr < 0.05:
        return "*"
    return "ns"


def _score_violin_figure(
    obs, values, celltype_key, condition_key, reference, comparison, title, *,
    statistics=None, feature=None, show_sample_fdr=True,
):
    import matplotlib.pyplot as plt

    groups = sorted(obs[celltype_key].astype(str).unique())
    if not groups:
        return None
    conditions = (reference, comparison)
    fig, axis = plt.subplots(figsize=(max(9.0, 1.55 * len(groups) + 4.0), 5.3), dpi=150)
    positions = np.arange(len(groups), dtype=float)
    width = 0.30
    values = np.asarray(values, dtype=float)
    group_maxima = {}
    all_values = []
    for condition_index, condition in enumerate(conditions):
        datasets, used_positions = [], []
        for group_index, group in enumerate(groups):
            mask = obs[celltype_key].astype(str).eq(group) & obs[condition_key].astype(str).eq(condition)
            current = values[np.asarray(mask)]
            current = current[np.isfinite(current)]
            if len(current):
                datasets.append(current)
                used_positions.append(positions[group_index] + (condition_index - 0.5) * width)
                group_maxima[group] = max(group_maxima.get(group, -np.inf), float(np.max(current)))
                all_values.extend(current.tolist())
        if not datasets:
            continue
        violin = axis.violinplot(datasets, positions=used_positions, widths=width * 0.88,
                                 showmedians=True, showextrema=False)
        color = NATURE_PALETTE[condition_index % len(NATURE_PALETTE)]
        for body in violin["bodies"]:
            body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(0.55)
        violin["cmedians"].set_color(color)
        axis.plot([], [], color=color, linewidth=7, alpha=0.55, label=condition)

    if show_sample_fdr and feature and all_values:
        y_min, y_max = float(np.min(all_values)), float(np.max(all_values))
        span = max(y_max - y_min, 1e-6)
        bracket_height = 0.030 * span
        top = y_max
        for group_index, group in enumerate(groups):
            if group not in group_maxima:
                continue
            label = _sample_fdr_label(statistics, feature, group)
            baseline = group_maxima[group] + 0.075 * span
            x_left = positions[group_index] - 0.5 * width
            x_right = positions[group_index] + 0.5 * width
            axis.plot(
                [x_left, x_left, x_right, x_right],
                [baseline, baseline + bracket_height, baseline + bracket_height, baseline],
                color="#344054", linewidth=0.8, clip_on=False,
            )
            axis.text(
                positions[group_index], baseline + bracket_height + 0.012 * span, label,
                ha="center", va="bottom", fontsize=8, color="#344054",
            )
            top = max(top, baseline + bracket_height + 0.060 * span)
        axis.set_ylim(top=top)
    axis.set_xticks(positions, [_wrapped_label(group, width=18) for group in groups])
    axis.tick_params(axis="x", labelsize=8, pad=6)
    axis.set_ylabel("Cell-level score (descriptive)")
    axis.set_xlabel("Cell type")
    axis.set_title(title, loc="left", fontsize=10, fontweight="semibold")
    axis.legend(frameon=False, fontsize=8)
    axis.grid(axis="y", color="#E4E7EC", linewidth=0.55)
    fig.tight_layout(pad=1.1)
    return fig


def _cell_level_correlation_figure(obs, feature_values, x_feature, y_feature, condition_key, title):
    """Plot a cell-level association explicitly as descriptive evidence only."""
    from scipy.stats import spearmanr
    import matplotlib.pyplot as plt

    frame = pd.DataFrame({
        "condition": obs[condition_key].astype(str).to_numpy(),
        "x": np.asarray(feature_values[x_feature], dtype=float),
        "y": np.asarray(feature_values[y_feature], dtype=float),
    }).replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 3:
        return None
    result = spearmanr(frame["x"].to_numpy(float), frame["y"].to_numpy(float))
    rho, pvalue = float(result.statistic), float(result.pvalue)
    fig, axis = plt.subplots(figsize=(7.5, 5.15), dpi=150)
    for index, condition in enumerate(sorted(frame["condition"].unique())):
        local = frame[frame["condition"].eq(condition)]
        axis.scatter(local["x"], local["y"], s=4.2, alpha=0.30,
                     color=NATURE_PALETTE[index % len(NATURE_PALETTE)], label=condition,
                     edgecolor="none", rasterized=True)
    axis.text(
        0.03, 0.97,
        "Cell-level descriptive only\n"
        f"Spearman ρ = {rho:.2f} ({_spearman_effect_size(rho)})\n"
        f"n = {len(frame):,}; P omitted (cells are non-independent)",
        transform=axis.transAxes, va="top", fontsize=8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D0D5DD", "alpha": 0.92},
    )
    axis.set_xlabel(x_feature); axis.set_ylabel(y_feature)
    axis.set_title(title, loc="left", fontsize=10, fontweight="semibold")
    axis.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    axis.grid(color="#E4E7EC", linewidth=0.55)
    fig.tight_layout(pad=1.1, rect=(0.0, 0.0, 0.80, 1.0))
    return fig


def _correlation_figure(sample_table, x_feature, y_feature, title):
    """Show a sample × cell-type scatter without calling repeated units independent."""
    from scipy.stats import spearmanr
    import matplotlib.pyplot as plt

    eligible = sample_table[sample_table["eligible_for_statistics"]].copy()
    pivot = eligible.pivot_table(
        index=["sample_id", "condition", "celltype"], columns="feature", values="score", aggfunc="first",
    ).reset_index()
    if x_feature not in pivot or y_feature not in pivot:
        return None
    frame = pivot[["sample_id", "condition", "celltype", x_feature, y_feature]].dropna()
    if len(frame) < 3:
        return None
    fig, axis = plt.subplots(figsize=(7.5, 5.15), dpi=150)
    for index, condition in enumerate(sorted(frame["condition"].astype(str).unique())):
        local = frame[frame["condition"].astype(str).eq(condition)]
        axis.scatter(local[x_feature], local[y_feature], s=26, alpha=0.78,
                     color=NATURE_PALETTE[index % len(NATURE_PALETTE)], label=condition,
                     edgecolor="white", linewidth=0.35)
    x, y = frame[x_feature].to_numpy(float), frame[y_feature].to_numpy(float)
    result = spearmanr(x, y)
    rho = float(result.statistic)
    axis.set_xlabel(x_feature); axis.set_ylabel(y_feature)
    axis.set_title(title, loc="left", fontsize=10, fontweight="semibold")
    axis.text(
        0.03, 0.97,
        "Exploratory aggregate only\n"
        f"Spearman ρ = {rho:.2f} ({_spearman_effect_size(rho)})\n"
        f"{frame['sample_id'].nunique()} donors; repeated across cell types\n"
        "P/FDR omitted; adjusted mixed model required",
        transform=axis.transAxes, va="top", fontsize=8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D0D5DD", "alpha": 0.92},
    )
    axis.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    axis.grid(color="#E4E7EC", linewidth=0.55)
    fig.tight_layout(pad=1.1, rect=(0.0, 0.0, 0.80, 1.0))
    return fig


def _within_celltype_cell_level_correlation_figure(
    obs, feature_values, x_feature, y_feature, *, celltype_key, condition_key,
):
    """Faceted cell-level views distinguish within-state patterns from mixing states."""
    from scipy.stats import spearmanr
    import matplotlib.pyplot as plt

    frame = pd.DataFrame({
        "celltype": obs[celltype_key].astype(str).to_numpy(),
        "condition": obs[condition_key].astype(str).to_numpy(),
        "x": np.asarray(feature_values[x_feature], dtype=float),
        "y": np.asarray(feature_values[y_feature], dtype=float),
    }).replace([np.inf, -np.inf], np.nan).dropna()
    celltypes = sorted(frame["celltype"].unique())
    if not celltypes:
        return None
    ncols = min(2, len(celltypes))
    nrows = int(np.ceil(len(celltypes) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.7 * ncols, 3.75 * nrows + 0.55), dpi=150,
        sharex=True, sharey=True, squeeze=False,
    )
    for index, celltype in enumerate(celltypes):
        axis = axes.ravel()[index]
        local = frame[frame["celltype"].eq(celltype)]
        for condition_index, condition in enumerate(sorted(local["condition"].unique())):
            subset = local[local["condition"].eq(condition)]
            axis.scatter(
                subset["x"], subset["y"], s=3.6, alpha=0.30,
                color=NATURE_PALETTE[condition_index % len(NATURE_PALETTE)],
                label=condition if index == 0 else None, edgecolor="none", rasterized=True,
            )
        rho = np.nan
        if len(local) >= 3:
            rho = float(spearmanr(local["x"], local["y"]).statistic)
        axis.text(
            0.03, 0.97, f"n = {len(local):,}\nρ = {rho:.2f} ({_spearman_effect_size(rho)})",
            transform=axis.transAxes, va="top", fontsize=7.2,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D5DD", "alpha": 0.90},
        )
        axis.set_title(_wrapped_label(celltype, width=19), loc="left", fontsize=9, fontweight="semibold")
        axis.grid(color="#E4E7EC", linewidth=0.50)
    for axis in axes.ravel()[len(celltypes):]:
        axis.set_visible(False)
    fig.suptitle("Within-cell-type descriptive functional-score correlations", x=0.01, ha="left",
                 fontsize=10, fontweight="semibold")
    # Reserve a true footer below the lowest tick labels.  Without that space,
    # the shared x label can collide with the lower row in a 2×2 facet grid.
    fig.text(0.5, 0.004, x_feature, ha="center", va="bottom", fontsize=8.5)
    fig.text(0.015, 0.5, y_feature, va="center", rotation="vertical", fontsize=9)
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.995, 0.98),
                   frameon=False, fontsize=8)
    fig.tight_layout(pad=1.05, rect=(0.05, 0.09, 0.94, 0.93))
    return fig


def _within_celltype_condition_cell_level_correlation_figure(
    obs, feature_values, x_feature, y_feature, *, celltype_key, condition_key, sample_key,
    reference, comparison,
):
    """Show each within-cell-type cell-level association separately by condition."""
    from scipy.stats import spearmanr
    import matplotlib.pyplot as plt

    frame = pd.DataFrame({
        "sample_id": obs[sample_key].astype(str).to_numpy(),
        "celltype": obs[celltype_key].astype(str).to_numpy(),
        "condition": obs[condition_key].astype(str).to_numpy(),
        "x": np.asarray(feature_values[x_feature], dtype=float),
        "y": np.asarray(feature_values[y_feature], dtype=float),
    }).replace([np.inf, -np.inf], np.nan).dropna()
    conditions = [condition for condition in (reference, comparison) if condition in set(frame["condition"])]
    celltypes = sorted(frame["celltype"].unique())
    if not conditions or not celltypes:
        return None

    fig = plt.figure(figsize=(4.55 * len(conditions), 3.05 * len(celltypes) + 0.95), dpi=150)
    grid = fig.add_gridspec(len(celltypes), len(conditions), wspace=0.12, hspace=0.24)
    axes = np.empty((len(celltypes), len(conditions)), dtype=object)
    for row_index, celltype in enumerate(celltypes):
        for column_index, condition in enumerate(conditions):
            shared = axes[0, 0] if row_index or column_index else None
            axis = fig.add_subplot(
                grid[row_index, column_index],
                sharex=shared if shared is not None else None,
                sharey=shared if shared is not None else None,
            )
            axes[row_index, column_index] = axis
            local = frame[frame["celltype"].eq(celltype) & frame["condition"].eq(condition)]
            color = NATURE_PALETTE[column_index % len(NATURE_PALETTE)]
            axis.scatter(local["x"], local["y"], s=3.6, alpha=0.30, color=color,
                         edgecolor="none", rasterized=True)
            rho = np.nan
            if len(local) >= 3:
                rho = float(spearmanr(local["x"], local["y"]).statistic)
            axis.text(
                0.03, 0.97,
                f"n cells = {len(local):,}\nn donors = {local['sample_id'].nunique()}\n"
                f"ρ = {rho:.2f} ({_spearman_effect_size(rho)})",
                transform=axis.transAxes, va="top", fontsize=7.0,
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#D0D5DD", "alpha": 0.90},
            )
            axis.set_title(
                f"{_wrapped_label(celltype, width=17)}\n{condition}",
                loc="left", fontsize=8.4, fontweight="semibold",
            )
            axis.grid(color="#E4E7EC", linewidth=0.50)
            if row_index < len(celltypes) - 1:
                axis.tick_params(labelbottom=False)
            if column_index > 0:
                axis.tick_params(labelleft=False)
    fig.suptitle(
        "Within-cell-type functional-score correlations, stratified by condition",
        x=0.01, ha="left", fontsize=10, fontweight="semibold",
    )
    fig.text(
        0.01, 0.955,
        "Cell-level descriptive views; P omitted because cells are non-independent.",
        ha="left", va="top", fontsize=7.5, color="#475467",
    )
    fig.text(0.5, 0.010, x_feature, ha="center", va="bottom", fontsize=8.5)
    fig.text(0.012, 0.5, y_feature, va="center", rotation="vertical", fontsize=8.5)
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.075, top=0.91)
    return fig


def _condition_sample_celltype_correlation_figure(
    sample_table, x_feature, y_feature, *, reference, comparison,
):
    """Compare repeated-measures aggregates separately within each condition."""
    from scipy.stats import spearmanr
    import matplotlib.pyplot as plt

    eligible = sample_table[sample_table["eligible_for_statistics"]].copy()
    pivot = eligible.pivot_table(
        index=["sample_id", "condition", "celltype"], columns="feature", values="score", aggfunc="first",
    ).reset_index()
    if x_feature not in pivot or y_feature not in pivot:
        return None
    frame = pivot[["sample_id", "condition", "celltype", x_feature, y_feature]].dropna()
    conditions = [condition for condition in (reference, comparison) if condition in set(frame["condition"].astype(str))]
    celltypes = sorted(frame["celltype"].astype(str).unique())
    if not conditions or not celltypes:
        return None

    markers = ("o", "s", "^", "D", "P", "X", "v", "<")
    fig, axes = plt.subplots(1, len(conditions), figsize=(5.10 * len(conditions), 5.10), dpi=150,
                             sharex=True, sharey=True, squeeze=False)
    for condition_index, condition in enumerate(conditions):
        axis = axes[0, condition_index]
        local = frame[frame["condition"].astype(str).eq(condition)]
        for celltype_index, celltype in enumerate(celltypes):
            subset = local[local["celltype"].astype(str).eq(celltype)]
            axis.scatter(
                subset[x_feature], subset[y_feature], s=27, alpha=0.82,
                color=NATURE_PALETTE[celltype_index % len(NATURE_PALETTE)],
                marker=markers[celltype_index % len(markers)],
                label=celltype if condition_index == 0 else None,
                edgecolor="white", linewidth=0.35,
            )
        rho = np.nan
        if len(local) >= 3:
            rho = float(spearmanr(local[x_feature], local[y_feature]).statistic)
        axis.text(
            0.03, 0.97,
            "Condition-stratified exploratory aggregate\n"
            f"Spearman ρ = {rho:.2f} ({_spearman_effect_size(rho)})\n"
            f"{local['sample_id'].nunique()} donors; {len(local)} sample × cell-type units\n"
            "P/FDR omitted; donor-aware model required",
            transform=axis.transAxes, va="top", fontsize=7.5,
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#D0D5DD", "alpha": 0.92},
        )
        axis.set_title(str(condition), loc="left", fontsize=9.5, fontweight="semibold")
        axis.grid(color="#E4E7EC", linewidth=0.55)
    fig.suptitle("Sample × cell-type functional-score correlation, stratified by condition", x=0.01,
                 ha="left", fontsize=10, fontweight="semibold")
    fig.text(0.5, 0.010, x_feature, ha="center", va="bottom", fontsize=8.5)
    fig.text(0.012, 0.5, y_feature, va="center", rotation="vertical", fontsize=8.5)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        legend_axis = fig.add_axes((0.10, 0.045, 0.86, 0.095))
        legend_axis.set_axis_off()
        legend_axis.set_in_layout(False)
        legend_axis.legend(handles, labels, loc="center", frameon=False, fontsize=7.2,
                           ncol=min(2, len(handles)), title="Cell type", title_fontsize=7.2)
    # The common exporter applies tight_layout.  This layout contract reserves
    # a dedicated bottom strip for the cell-type key and shared x-axis title.
    fig._native_layout_rect = (0.0, 0.17, 1.0, 0.93)
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.22, top=0.91, wspace=0.16)
    return fig


def _select_typed_features(requested, feature_values, feature_types, expected_type, fallback):
    """Resolve an ordered user panel while silently excluding unavailable scores."""
    requested = _parse_feature_selection(requested) or tuple(fallback)
    selected, missing = [], []
    for item in requested:
        feature = _resolve_feature_name(item, feature_values)
        if feature and feature_types.get(feature) == expected_type:
            if feature not in selected:
                selected.append(feature)
        else:
            missing.append(item)
    return selected, missing


def _regulator_pathway_concordance(sample_table, tf_features, pathway_features, min_samples, *, condition=None):
    """Compute the transparent, repeated-measures-aware TF × pathway matrix."""
    from scipy.stats import spearmanr

    columns = [
        "condition", "tf_feature", "pathway_feature", "n_sample_celltype_units", "n_unique_samples",
        "rho", "effect_size_interpretation", "pvalue", "fdr_bh",
        "status", "analysis_note",
    ]
    if not tf_features or not pathway_features:
        return pd.DataFrame(columns=columns)
    needed = set(tf_features).union(pathway_features)
    eligible = sample_table[
        sample_table["eligible_for_statistics"] & sample_table["feature"].isin(needed)
    ].copy()
    condition_label = "all_conditions" if condition is None else str(condition)
    if condition is not None:
        eligible = eligible[eligible["condition"].astype(str).eq(condition_label)].copy()
    pivot = eligible.pivot_table(
        index=["sample_id", "condition", "celltype"], columns="feature", values="score", aggfunc="first",
    )
    required_units = max(3, 2 * int(min_samples))
    rows = []
    for tf_feature in tf_features:
        for pathway_feature in pathway_features:
            if tf_feature not in pivot or pathway_feature not in pivot:
                continue
            frame = pivot[[tf_feature, pathway_feature]].replace([np.inf, -np.inf], np.nan).dropna()
            rho, pvalue = np.nan, np.nan
            if len(frame) >= 3:
                result = spearmanr(frame[tf_feature], frame[pathway_feature])
                rho, pvalue = float(result.statistic), float(result.pvalue)
            rows.append({
                "condition": condition_label, "tf_feature": tf_feature, "pathway_feature": pathway_feature,
                "n_sample_celltype_units": int(len(frame)),
                "n_unique_samples": int(frame.index.get_level_values("sample_id").nunique()),
                "rho": rho, "effect_size_interpretation": _spearman_effect_size(rho),
                "pvalue": pvalue, "fdr_bh": np.nan,
                "status": (
                    ("exploratory_repeated_measures_condition_stratified" if condition is not None
                     else "exploratory_repeated_measures")
                    if len(frame) >= required_units else "insufficient_sample_celltype_units"
                ),
                "analysis_note": (
                    "Condition-stratified sample × cell-type units reuse donors across cell types; "
                    "effect sizes are exploratory and require a donor-aware mixed model for adjusted inference."
                    if condition is not None else
                    "Sample × cell-type units reuse donors across cell types; effect sizes are "
                    "exploratory and require a donor-aware mixed model for adjusted inference."
                ),
            })
    result = pd.DataFrame(rows, columns=columns)
    # Retain raw P only as an audit field.  A BH-adjusted P would visually
    # imply independent units, so it is deliberately not calculated here.
    return result


def _regulator_pathway_concordance_figure(concordance, tf_features, pathway_features, tf_labels=None):
    """Render a readable effect-size matrix, without significance stars."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if concordance.empty or not tf_features or not pathway_features:
        return None
    matrix = np.full((len(tf_features), len(pathway_features)), np.nan, dtype=float)
    for row_index, tf_feature in enumerate(tf_features):
        for column_index, pathway_feature in enumerate(pathway_features):
            rows = concordance[
                concordance["tf_feature"].eq(tf_feature)
                & concordance["pathway_feature"].eq(pathway_feature)
            ]
            if not rows.empty:
                matrix[row_index, column_index] = float(rows.iloc[0]["rho"])
    cmap = LinearSegmentedColormap.from_list(
        "functional_concordance", ["#2166AC", "#F7F7F7", "#B2182B"], N=256,
    )
    # Pathway labels can be long (especially managed Hallmark identifiers).
    # Allocate enough horizontal room for one readable wrapped label per score.
    width = max(9.5, min(20.0, 4.0 + 1.35 * len(pathway_features)))
    height = max(4.8, min(11.0, 2.8 + 0.43 * len(tf_features)))
    fig, axis = plt.subplots(figsize=(width, height), dpi=150)
    image = axis.imshow(matrix, aspect="auto", cmap=cmap, vmin=-1, vmax=1)
    labels = tf_labels or {}
    axis.set_yticks(np.arange(len(tf_features)), [labels.get(item, item) for item in tf_features], fontsize=8)
    axis.set_xticks(
        np.arange(len(pathway_features)),
        [_pathway_display_label(item, width=18) for item in pathway_features],
        fontsize=8,
    )
    axis.tick_params(axis="x", pad=5)
    for x in range(len(pathway_features) + 1):
        axis.axvline(x - 0.5, color="white", linewidth=0.70, alpha=0.85)
    for y in range(len(tf_features) + 1):
        axis.axhline(y - 0.5, color="white", linewidth=0.70, alpha=0.85)
    if matrix.size <= 56:
        for row_index in range(matrix.shape[0]):
            for column_index in range(matrix.shape[1]):
                value = matrix[row_index, column_index]
                if np.isfinite(value):
                    axis.text(
                        column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(value) >= 0.55 else "#344054",
                    )
    axis.set_xlabel("Pathway score", fontsize=9)
    axis.set_ylabel("Network-backed TF activity", fontsize=9)
    fig.suptitle("Regulator–pathway concordance", x=0.01, y=0.99, ha="left",
                 fontsize=10, fontweight="semibold")
    fig.text(
        0.01, 0.945,
        "Exploratory sample × cell-type means; same donor contributes repeated cell-type units. "
        "Colours and values are Spearman ρ, not significance claims.",
        ha="left", va="top", fontsize=7.5, color="#475467",
    )
    colorbar = fig.colorbar(image, ax=axis, fraction=0.035, pad=0.025, aspect=30)
    colorbar.outline.set_visible(False)
    colorbar.set_label("Spearman ρ", fontsize=8, labelpad=5)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    fig.tight_layout(pad=1.1, rect=(0.0, 0.0, 1.0, 0.91))
    return fig


def _regulator_pathway_condition_concordance_figure(
    concordance, tf_features, pathway_features, *, reference, comparison, tf_labels=None,
):
    """Render condition-stratified TF × pathway matrices on one colour scale."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if concordance.empty or "condition" not in concordance or not tf_features or not pathway_features:
        return None
    conditions = [condition for condition in (reference, comparison) if bool(
        concordance["condition"].astype(str).eq(str(condition)).any()
    )]
    if not conditions:
        return None
    cmap = LinearSegmentedColormap.from_list(
        "functional_concordance_condition", ["#2166AC", "#F7F7F7", "#B2182B"], N=256,
    )
    # Stack conditions vertically: with many pathway columns, a side-by-side
    # layout leaves no room for readable x-axis labels.
    width = max(9.5, min(20.0, 4.0 + 1.35 * len(pathway_features)))
    height = max(5.9, 2.45 * len(conditions) + 2.5)
    fig = plt.figure(figsize=(width, height), dpi=150)
    grid = fig.add_gridspec(len(conditions), 2, width_ratios=[1.0, 0.035], hspace=0.52, wspace=0.05)
    axes = [fig.add_subplot(grid[index, 0]) for index in range(len(conditions))]
    colorbar_axis = fig.add_subplot(grid[:, 1])
    labels = tf_labels or {}
    image = None
    for index, (axis, condition) in enumerate(zip(axes, conditions)):
        local = concordance[concordance["condition"].astype(str).eq(str(condition))]
        matrix = np.full((len(tf_features), len(pathway_features)), np.nan, dtype=float)
        for row_index, tf_feature in enumerate(tf_features):
            for column_index, pathway_feature in enumerate(pathway_features):
                rows = local[
                    local["tf_feature"].eq(tf_feature)
                    & local["pathway_feature"].eq(pathway_feature)
                ]
                if not rows.empty:
                    matrix[row_index, column_index] = float(rows.iloc[0]["rho"])
        image = axis.imshow(matrix, aspect="auto", cmap=cmap, vmin=-1, vmax=1)
        axis.set_xticks(
            np.arange(len(pathway_features)),
            [_pathway_display_label(item, width=18) for item in pathway_features], fontsize=7.6,
        )
        axis.tick_params(axis="x", pad=4)
        if index == 0:
            axis.set_yticks(np.arange(len(tf_features)), [labels.get(item, item) for item in tf_features], fontsize=7.8)
            axis.set_ylabel("Network-backed TF activity", fontsize=8.5)
        else:
            axis.set_yticks(np.arange(len(tf_features)), [])
        n_donors = int(local["n_unique_samples"].max()) if not local.empty else 0
        axis.set_title(f"{condition}\n{n_donors} donors; repeated cell-type units", loc="left",
                       fontsize=8.8, fontweight="semibold")
        for x in range(len(pathway_features) + 1):
            axis.axvline(x - 0.5, color="white", linewidth=0.70, alpha=0.85)
        for y in range(len(tf_features) + 1):
            axis.axhline(y - 0.5, color="white", linewidth=0.70, alpha=0.85)
        if matrix.size <= 56:
            for row_index in range(matrix.shape[0]):
                for column_index in range(matrix.shape[1]):
                    value = matrix[row_index, column_index]
                    if np.isfinite(value):
                        axis.text(
                            column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=6.7,
                            color="white" if abs(value) >= 0.55 else "#344054",
                        )
        axis.set_xlabel("Pathway score", fontsize=8.5)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    if image is None:
        return None
    colorbar = fig.colorbar(image, cax=colorbar_axis)
    colorbar.outline.set_visible(False)
    colorbar.set_label("Spearman ρ", fontsize=8, labelpad=5)
    colorbar.ax.tick_params(labelsize=7.5)
    fig.suptitle("Regulator–pathway concordance, stratified by condition", x=0.01, y=0.995,
                 ha="left", fontsize=10, fontweight="semibold")
    fig.text(
        0.01, 0.946,
        "Exploratory condition-stratified sample × cell-type means; P/FDR omitted because donors repeat across cell types.",
        ha="left", va="top", fontsize=7.2, color="#475467",
    )
    fig.subplots_adjust(left=0.10, right=0.965, bottom=0.13, top=0.84)
    return fig


class FunctionalStateAnalysis(BaseAnalysis):
    """Integrate gene, TF-network and pathway-state evidence without pseudoreplication."""

    MODULE_NAME = "functional_state"
    DISPLAY_NAME = "转录调控与功能状态"
    DESCRIPTION = "基因表达、网络支持的 TF 活性和通路评分；样本级条件比较"
    # _validate_obs_column(..., "细胞类型") 硬性要求 celltype 列，因此输入
    # 选择器必须只提供已注释的上游输出，而不是让用户在运行期才遇到报错。
    INPUT_REQUIRES = ['celltype']

    def validate_input(self, adata):
        for key, label, multiple in (
            (self.params.get("sample_key", "sample_id"), "生物学样本", False),
            (self.params.get("condition_key", "condition"), "条件", True),
            (self.params.get("celltype_key", "celltype"), "细胞类型", False),
        ):
            try:
                _validate_obs_column(adata, key, label, multiple)
            except ValueError as exc:
                return str(exc)
        return None

    def run(self, input_path):
        self.progress(5, "加载功能状态分析输入...")
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata, self.MODULE_NAME)
        if adata.n_obs == 0:
            raise ValueError("分析范围筛选后没有可用细胞。")

        sample_key = _validate_obs_column(adata, self.params.get("sample_key", "sample_id"), "生物学样本")
        condition_key = _validate_obs_column(adata, self.params.get("condition_key", "condition"), "条件", True)
        celltype_key = _validate_obs_column(adata, self.params.get("celltype_key", "celltype"), "细胞类型")
        organism = str(self.params.get("organism", "Human") or "Human").strip()
        if organism != "Human":
            raise ValueError("当前内置功能面板、TF 网络和基因集命名仅支持 Human HGNC gene symbols。")
        technical_names = {"batch", "technical_batch", "sequencing_batch", "library_batch"}
        if (sample_key.lower() in technical_names
                and not _as_bool(self.params.get("confirm_batch_is_biological_sample"), False)):
            raise ValueError(
                f"样本列 '{sample_key}' 看起来是技术 batch；未确认其代表独立生物学样本，不能报告样本级统计。"
            )
        if (condition_key.lower() in technical_names
                and not _as_bool(self.params.get("confirm_batch_is_biological_condition"), False)):
            raise ValueError(
                f"条件列 '{condition_key}' 看起来是技术 batch；未确认其代表生物学条件。"
            )
        design = adata.obs[[sample_key, condition_key]].astype(str).drop_duplicates()
        if bool(design.groupby(sample_key, observed=True)[condition_key].nunique().gt(1).any()):
            raise ValueError("同一个 sample_id 不能同时属于多个 condition。")
        per_sample = adata.obs[sample_key].astype(str).value_counts()
        if (per_sample.size < 2 or int(per_sample.median()) < 2
                or per_sample.size > max(20, int(0.5 * adata.n_obs))):
            raise ValueError("sample_key 看起来不是可用的生物学重复列。")

        reference, comparison = _choose_comparison(
            adata.obs[condition_key].astype(str).unique(),
            self.params.get("reference_condition", ""), self.params.get("comparison_condition", ""),
        )
        min_cells = max(1, int(self.params.get("min_cells_per_sample_celltype", 20) or 20))
        min_samples = max(2, int(self.params.get("min_samples_per_group", 3) or 3))
        statistic_method = str(self.params.get("statistic_method", "welch") or "welch").strip().lower()
        statistic_method = "mann_whitney" if statistic_method == "mann_whitney" else "welch"
        expression, expression_source, expression_qc = _resolve_expression(
            adata,
            self.params.get("expression_layer", "X"),
            self.params.get("expression_source", "auto"),
        )
        var_names = [str(gene) for gene in adata.var_names]
        gene_lookup = {gene.upper(): gene for gene in var_names}

        self.progress(15, "整理功能基因集与覆盖度...")
        pathway_scoring_method = str(
            self.params.get("pathway_scoring_method", "scanpy_score_genes") or "scanpy_score_genes"
        ).strip()
        if pathway_scoring_method != "scanpy_score_genes":
            raise ValueError("当前功能状态模块只支持 scanpy_score_genes 通路评分。")
        include_builtin_panels = _as_bool(self.params.get("include_builtin_pathway_panels"), True)
        analysis_focus = str(self.params.get("analysis_focus", "custom") or "custom").strip()
        if analysis_focus not in {"custom", "ibd_organoid_epithelial"}:
            raise ValueError("未知的功能状态分析重点。")
        pathway_sets = (
            {name: tuple(genes) for name, genes in BUILTIN_PATHWAYS.items()}
            if include_builtin_panels else {}
        )
        pathway_sources = {name: "builtin_lipid_inflammation" for name in pathway_sets}
        if analysis_focus == "ibd_organoid_epithelial":
            pathway_sets.update({
                name: tuple(genes) for name, genes in IBD_EPITHELIAL_QUICK_PATHWAYS.items()
            })
            pathway_sources.update({
                name: "builtin_ibd_epithelial_hypothesis"
                for name in IBD_EPITHELIAL_QUICK_PATHWAYS
            })
        managed_term_requests = [
            self.params.get("managed_gene_set_terms", ""),
            self.params.get("managed_gene_set_terms_advanced", ""),
        ]
        if analysis_focus == "ibd_organoid_epithelial":
            managed_term_requests.append(";".join(IBD_EPITHELIAL_MANAGED_TERMS))
        managed_sets, managed_provenance = load_selected_managed_gene_sets([
            *managed_term_requests,
        ])
        duplicated_managed = sorted(set(pathway_sets).intersection(managed_sets))
        if duplicated_managed:
            raise ValueError("托管基因集不能覆盖内置名称: " + "、".join(duplicated_managed))
        pathway_sets.update(managed_sets)
        pathway_sources.update({
            name: managed_provenance[name]["source"] for name in managed_sets
        })
        custom_sets = _parse_custom_gene_sets(self.params.get("custom_gene_sets", ""))
        duplicated_custom = sorted(set(pathway_sets).intersection(custom_sets))
        if duplicated_custom:
            raise ValueError("自定义基因集不能覆盖内置名称: " + "、".join(duplicated_custom))
        pathway_sets.update(custom_sets)
        pathway_sources.update({name: "custom_inline" for name in custom_sets})
        gmt_path = _resolve_project_file(self.project_dir, self.params.get("pathway_gmt_path", ""), "通路 GMT 文件")
        gmt_sets = _read_selected_gmt(gmt_path, self.params.get("pathway_gmt_terms", ""))
        duplicated_gmt = sorted(set(pathway_sets).intersection(gmt_sets))
        if duplicated_gmt:
            raise ValueError("GMT term 不能覆盖已有通路名称: " + "、".join(duplicated_gmt))
        pathway_sets.update(gmt_sets)
        pathway_sources.update({name: "project_local_gmt" for name in gmt_sets})
        pathway_overlap_qc = _pathway_gene_set_overlap_qc(pathway_sets, gene_lookup)
        min_general_genes = max(2, int(self.params.get("min_gene_set_genes", 5) or 5))
        coverage_rows, score_values, feature_types, feature_definitions = [], {}, {}, {}
        score_provenance = {}
        for name, genes in pathway_sets.items():
            detected = tuple(gene_lookup[gene] for gene in genes if gene in gene_lookup)
            missing = tuple(gene for gene in genes if gene not in gene_lookup)
            is_builtin_small = name in {
                CARNITINE_SHUTTLE_MODULE_NAME,
                INTRACELLULAR_FA_TRANSPORT_MODULE_NAME,
                "Ketogenesis",
            }
            source = pathway_sources[name]
            is_managed_standard = source.startswith("managed_")
            signature_kind = "standard_pathway"
            if source == "builtin_ibd_epithelial_hypothesis":
                signature_kind = "focused_epithelial_hypothesis_signature"
            elif len(genes) < min_general_genes:
                signature_kind = (
                    "exploratory_custom_signature" if source == "custom_inline"
                    else "small_signature"
                )
            # The three curated short modules and explicit inline signatures
            # may be useful hypotheses at 2+ detected genes, but they must
            # remain visibly distinct from a standard pathway score.
            threshold = 2 if (
                is_builtin_small
                or signature_kind in {
                    "exploratory_custom_signature", "focused_epithelial_hypothesis_signature",
                }
            ) else min_general_genes
            if is_managed_standard:
                threshold = max(threshold, 10)
            coverage_fraction = len(detected) / max(len(genes), 1)
            coverage_level = (
                "good_coverage" if coverage_fraction >= 0.70
                else "coverage_warning" if coverage_fraction >= 0.40
                else "unreliable_coverage"
            )
            status = "good"
            if len(detected) < threshold:
                status = "insufficient_detected_genes"
            elif len(genes) >= 10 and coverage_fraction < 0.40:
                # Do not turn a tiny fraction of a standard pathway into a
                # deceptively precise score.  Small curated hypotheses retain
                # their explicit exploratory/small-signature status below.
                status = "unreliable_coverage_not_scored"
            elif coverage_fraction < 0.70 and len(genes) >= 10:
                status = "coverage_warning"
            elif len(genes) < min_general_genes:
                status = signature_kind
            coverage_rows.append({
                "gene_set": name, "source": source, "signature_kind": signature_kind,
                "n_original": len(genes), "n_detected": len(detected),
                "coverage": coverage_fraction, "coverage_level": coverage_level,
                "minimum_detected_genes": threshold, "missing_genes": ";".join(missing),
                "status": status, "coverage_status": coverage_level,
            })
            if status in {"insufficient_detected_genes", "unreliable_coverage_not_scored"}:
                continue
            values, method = _score_gene_set(
                expression, var_names, detected, score_name=name,
                requested_method=pathway_scoring_method,
            )
            feature = f"{name} score"
            column = "fs_pathway_" + _safe_name(name) + "_score"
            adata.obs[column] = values
            score_values[feature] = values
            feature_types[feature] = "pathway_score"
            feature_definitions[feature] = tuple(gene.upper() for gene in detected)
            score_provenance[feature] = {
                "column": column, "method": method, "source": source,
                "evidence_class": signature_kind, "coverage": coverage_fraction,
                "coverage_level": coverage_level,
            }

        self.progress(35, "计算关键基因表达与可选 TF 靶基因活性...")
        custom_expression_genes = _parse_gene_list(self.params.get("expression_genes", ""))
        expression_genes = custom_expression_genes or DEFAULT_EXPRESSION_GENES
        expression_panel_sections = (
            {"Custom gene expression": tuple(expression_genes)}
            if custom_expression_genes else GENE_EXPRESSION_PANEL_SECTIONS
        )
        gene_to_panel_section = {
            gene: section
            for section, genes in expression_panel_sections.items()
            for gene in genes
        }
        expression_by_gene = {}
        for gene in expression_genes:
            actual = gene_lookup.get(gene)
            if actual is None:
                continue
            position = var_names.index(actual)
            values = _dense_columns(expression, [position]).reshape(-1)
            feature = f"{gene} expression"
            expression_by_gene[gene] = values
            score_values[feature] = values
            feature_types[feature] = "gene_expression"
            feature_definitions[feature] = (gene,)
            score_provenance[feature] = {"column": "", "method": "log_normalized_expression", "source": "gene_expression"}

        requested_tfs = _parse_gene_list(self.params.get("tf_panel", "")) or tuple(METABOLIC_TF_GENES + INFLAMMATORY_TF_GENES)
        network_path, network_source, network_resource_metadata = _resolve_tf_network(self.project_dir, self.params)
        tf_network = _read_tf_network(network_path)
        if tf_network.empty:
            tf_coverage = pd.DataFrame([{
                "tf": tf, "network": network_source, "n_network_targets": 0, "n_detected_targets": 0,
                "coverage": np.nan, "network_targets": 0, "detected_targets": 0,
                "coverage_fraction": np.nan, "missing_targets": "", "status": "network_not_configured",
                "activity_status": "network_not_configured", "display_label": tf,
                "minimum_required_targets": 5, "coverage_level": "network_not_configured",
                "score_column": "", "method": "not_run",
            } for tf in requested_tfs])
            tf_scores, tf_definitions = {}, {}
        else:
            tf_coverage, tf_scores, tf_definitions = _ulm_regulon_scores(
                expression, var_names, tf_network, requested_tfs, network_name=network_source,
            )
        for tf, values in tf_scores.items():
            feature = f"{tf} activity"
            column = "fs_tf_" + _safe_name(tf) + "_activity"
            adata.obs[column] = values
            score_values[feature] = values
            feature_types[feature] = "tf_activity"
            feature_definitions[feature] = tuple(gene.upper() for gene in tf_definitions.get(feature, ()))
            score_provenance[feature] = {"column": column, "method": "decoupler_ulm", "source": network_source,
                                         "score_interpretation": "ULM slope t-statistic"}
        heatmap_feature_labels = {}
        for row in tf_coverage.itertuples(index=False):
            feature = f"{row.tf} activity"
            if getattr(row, "activity_status", "") == "low_coverage":
                heatmap_feature_labels[feature] = f"{row.tf} activity*"

        if not score_values:
            raise ValueError("所有功能基因集均未命中表达矩阵，无法计算功能状态。")

        self.progress(55, "按 sample × celltype 汇总并进行样本级比较...")
        sample_table = _sample_feature_table(
            adata.obs, score_values, feature_types, sample_key=sample_key,
            condition_key=condition_key, celltype_key=celltype_key, min_cells=min_cells,
        )
        statistics = _feature_statistics(sample_table, reference, comparison, min_samples, statistic_method)
        gene_statistics = statistics[statistics["feature_type"].eq("gene_expression")].copy()
        tf_statistics = statistics[statistics["feature_type"].eq("tf_activity")].copy()
        pathway_statistics = statistics[statistics["feature_type"].eq("pathway_score")].copy()
        if not gene_statistics.empty:
            gene_statistics["panel_section"] = gene_statistics["feature"].str.replace(
                r" expression$", "", regex=True,
            ).map(gene_to_panel_section).fillna("Custom gene expression")
        tf_expression_activity = _tf_expression_activity_table(gene_statistics, tf_statistics, tf_coverage)

        x_requested = _requested_correlation_feature(self.params, "x")
        y_requested = _requested_correlation_feature(self.params, "y")
        correlation_levels = _parse_correlation_levels(self.params.get("correlation_levels", ""))
        correlations, overlap = _correlations(
            adata.obs, score_values, feature_definitions, sample_table,
            x_name=x_requested, y_name=y_requested, celltype_key=celltype_key,
            condition_key=condition_key, conditions=(reference, comparison), min_samples=min_samples,
            levels=correlation_levels,
        )
        x_feature = _resolve_feature_name(x_requested, score_values)
        y_feature = _resolve_feature_name(y_requested, score_values)
        concordance_tf_features, missing_concordance_tfs = _select_typed_features(
            self.params.get("concordance_tf_features", ""), score_values, feature_types,
            "tf_activity", DEFAULT_CONCORDANCE_TF_FEATURES,
        )
        concordance_pathway_features, missing_concordance_pathways = _select_typed_features(
            self.params.get("concordance_pathway_features", ""), score_values, feature_types,
            "pathway_score", DEFAULT_CONCORDANCE_PATHWAY_FEATURES,
        )
        concordance = _regulator_pathway_concordance(
            sample_table, concordance_tf_features, concordance_pathway_features, min_samples,
        )
        condition_concordance = pd.concat([
            _regulator_pathway_concordance(
                sample_table, concordance_tf_features, concordance_pathway_features, min_samples,
                condition=condition,
            )
            for condition in (reference, comparison)
        ], ignore_index=True)

        run_key = _safe_name(self.params.get("_analysis_id"), "manual_run")
        results_dir = Path(self.project_dir) / "results" / "functional_state" / run_key
        plots_dir = self.ensure_plots_dir()
        results_dir.mkdir(parents=True, exist_ok=True)
        result_files = []

        def save_table(frame, filename, label):
            path = results_dir / filename
            frame.to_csv(path, index=False)
            result_files.append({"file_path": str(path), "file_type": "csv", "category": "table", "label": label})

        save_table(sample_table, "sample_celltype_scores.csv", "Sample × cell type functional scores")
        save_table(statistics, "functional_state_statistics.csv", "Functional-state sample-level statistics")
        save_table(gene_statistics, "gene_expression_statistics.csv", "Key gene-expression sample-level statistics")
        save_table(tf_statistics, "tf_activity_statistics.csv", "TF activity sample-level statistics")
        save_table(tf_expression_activity, "tf_expression_vs_activity.csv", "TF expression versus ULM activity contrasts")
        save_table(pathway_statistics, "pathway_statistics.csv", "Pathway-score sample-level statistics")
        save_table(tf_coverage, "tf_activity_coverage.csv", "TF regulon coverage")
        save_table(pd.DataFrame(coverage_rows), "gene_set_coverage.csv", "Gene-set coverage")
        save_table(
            pathway_overlap_qc, "pathway_gene_set_overlap_qc.csv",
            "Pathway gene-set overlap QC (input-specific detected genes)",
        )
        save_table(correlations, "activity_correlations.csv", "Functional-state correlations (descriptive/exploratory)")
        save_table(
            concordance, "regulator_pathway_concordance.csv",
            "Regulator–pathway concordance (exploratory repeated-measures aggregate)",
        )
        save_table(
            condition_concordance, "regulator_pathway_concordance_by_condition.csv",
            "Regulator–pathway concordance by condition (exploratory repeated measures)",
        )
        if not overlap.empty:
            save_table(overlap, "gene_set_overlap.csv", "Correlation gene-set overlap QC")
        if _as_bool(self.params.get("export_cell_scores"), False):
            cell_scores = pd.DataFrame(index=adata.obs_names)
            cell_scores.index.name = "cell_id"
            for feature, values in score_values.items():
                if feature_types[feature] != "gene_expression":
                    cell_scores[_safe_name(feature)] = values
            cell_path = results_dir / "cell_activity_scores.csv"
            cell_scores.reset_index().to_csv(cell_path, index=False)
            result_files.append({"file_path": str(cell_path), "file_type": "csv", "category": "table", "label": "Cell-level activity scores (optional export)"})

        self.progress(72, "生成功能状态图表...")
        # Direct gene-expression evidence: mean expression and detection rate
        # remain visible together instead of being confused with activity.
        if expression_by_gene:
            categories = []
            mean_rows, detection_rows = [], []
            celltypes = sorted(adata.obs[celltype_key].astype(str).unique())
            for celltype in celltypes:
                for condition in (reference, comparison):
                    mask = adata.obs[celltype_key].astype(str).eq(celltype) & adata.obs[condition_key].astype(str).eq(condition)
                    if not bool(mask.any()):
                        continue
                    categories.append(f"{celltype}\n{condition}")
                    mean_rows.append([float(np.nanmean(values[np.asarray(mask)])) for values in expression_by_gene.values()])
                    detection_rows.append([float(np.mean(values[np.asarray(mask)] > 0)) for values in expression_by_gene.values()])
            if categories:
                figure = marker_dotplot_figure(np.asarray(mean_rows), np.asarray(detection_rows), categories,
                                               list(expression_by_gene), title="Key functional gene expression by cell type and condition",
                                               x_label="Cell type / condition", y_label="Gene",
                                               gene_sections=[gene_to_panel_section.get(gene, "Custom gene expression")
                                                              for gene in expression_by_gene])
                result_files.extend(self.save_matplotlib_figure(
                    figure, plots_dir, "functional_state_gene_expression_dotplot.png", "dotplot",
                    "Key functional gene expression", formats=("png", "svg"), dpi=300, preserve_aspect=True,
                ))

        heat_features, missing_heatmap_features, heatmap_selection_mode = _select_heatmap_features(
            self.params.get("heatmap_features", ""), score_values, feature_types,
        )
        if heat_features and not sample_table.empty:
            means = sample_table[sample_table["feature"].isin(heat_features)].groupby(
                ["feature", "celltype", "condition"], observed=True
            )["score"].mean().reset_index()
            ordered_columns = [
                (celltype, condition) for celltype in sorted(means["celltype"].astype(str).unique())
                for condition in (reference, comparison)
                if bool(((means["celltype"].astype(str) == celltype) & (means["condition"].astype(str) == condition)).any())
            ]
            matrix = np.full((len(heat_features), len(ordered_columns)), np.nan, dtype=float)
            for row_index, feature in enumerate(heat_features):
                for column_index, (celltype, condition) in enumerate(ordered_columns):
                    found = means[(means["feature"] == feature) & (means["celltype"].astype(str) == celltype) & (means["condition"].astype(str) == condition)]["score"]
                    if not found.empty:
                        matrix[row_index, column_index] = float(found.iloc[0])
            row_mean = np.nanmean(matrix, axis=1, keepdims=True)
            row_std = np.nanstd(matrix, axis=1, keepdims=True)
            z_matrix = np.divide(matrix - row_mean, row_std, out=np.zeros_like(matrix), where=row_std > 1e-10)
            heatmap_title = (
                "Selected functional-state scores by cell type and condition"
                if heatmap_selection_mode == "user_ordered_selection"
                else "Functional-state scores by cell type and condition"
            )
            figure = heatmap_figure(
                                    z_matrix, x_labels=[condition for _, condition in ordered_columns],
                                    x_group_labels=[_wrapped_label(celltype, width=18) for celltype, _ in ordered_columns],
                                    y_labels=[heatmap_feature_labels.get(feature, feature) for feature in heat_features], title=heatmap_title,
                                    x_label="Condition", y_label="Functional feature", colorbar_label="Row z-score",
                                    vmin=-3, vmax=3, x_label_rotation=0)
            result_files.extend(self.save_matplotlib_figure(
                figure, plots_dir, "functional_state_score_heatmap.png", "heatmap",
                "Functional-state sample-mean heatmap (row z-score; * low TF regulon coverage; grouped cell-type headers)", formats=("png", "svg"), dpi=300, preserve_aspect=True,
            ))

        if x_feature and y_feature:
            figure = _pair_umap_figure(adata, score_values, x_feature, y_feature)
            if figure is not None:
                result_files.extend(self.save_matplotlib_figure(
                    figure, plots_dir, "functional_state_paired_score_umap.png", "umap",
                    "Paired functional-score UMAP", formats=("png", "svg"), dpi=300, preserve_aspect=True,
                ))
            figure = _condition_umap_figure(adata, score_values[x_feature], condition_key, reference, comparison, x_feature)
            if figure is not None:
                result_files.extend(self.save_matplotlib_figure(
                    figure, plots_dir, "functional_state_condition_score_umap.png", "umap",
                    "Condition-split functional-score UMAP", formats=("png", "svg"), dpi=300, preserve_aspect=True,
                ))
            figure = _score_violin_figure(
                adata.obs, score_values[x_feature], celltype_key, condition_key, reference, comparison,
                f"{x_feature} by cell type and condition", statistics=statistics, feature=x_feature,
                show_sample_fdr=_as_bool(self.params.get("show_violin_sample_fdr"), True),
            )
            if figure is not None:
                result_files.extend(self.save_matplotlib_figure(
                    figure, plots_dir, "functional_state_score_violin.png", "violin",
                    "Functional score by cell type and condition (descriptive cells; sample-level FDR brackets)", formats=("png", "svg"), dpi=300, preserve_aspect=True,
                ))
            if (
                "cell_level_descriptive" in correlation_levels
                and _as_bool(self.params.get("show_cell_level_correlation"), True)
            ):
                figure = _cell_level_correlation_figure(
                    adata.obs, score_values, x_feature, y_feature, condition_key,
                    "Cell-level descriptive functional-score correlation",
                )
                if figure is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        figure, plots_dir, "functional_state_cell_level_score_correlation.png", "scatter",
                        "Cell-level descriptive functional-score correlation (Spearman rho effect-size interpretation; P omitted)",
                        formats=("png", "svg"), dpi=300, preserve_aspect=True,
                    ))
            if "within_celltype_cell_level_descriptive" in correlation_levels:
                figure = _within_celltype_cell_level_correlation_figure(
                    adata.obs, score_values, x_feature, y_feature,
                    celltype_key=celltype_key, condition_key=condition_key,
                )
                if figure is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        figure, plots_dir, "functional_state_within_celltype_cell_level_correlation.png", "scatter",
                        "Within-cell-type cell-level descriptive functional-score correlations",
                        formats=("png", "svg"), dpi=300, preserve_aspect=True,
                    ))
                figure = _within_celltype_condition_cell_level_correlation_figure(
                    adata.obs, score_values, x_feature, y_feature,
                    celltype_key=celltype_key, condition_key=condition_key, sample_key=sample_key,
                    reference=reference, comparison=comparison,
                )
                if figure is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        figure, plots_dir, "functional_state_within_celltype_condition_cell_level_correlation.png", "scatter",
                        "Within-cell-type cell-level functional-score correlations, stratified by condition (descriptive; P omitted)",
                        formats=("png", "svg"), dpi=300, preserve_aspect=True,
                    ))
            if "sample_x_celltype" in correlation_levels:
                figure = _correlation_figure(sample_table, x_feature, y_feature, "Sample × cell-type functional-score correlation")
                if figure is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        figure, plots_dir, "functional_state_score_correlation.png", "scatter",
                        "Sample × cell-type functional-score correlation (exploratory repeated measures)", formats=("png", "svg"), dpi=300, preserve_aspect=True,
                    ))
                figure = _condition_sample_celltype_correlation_figure(
                    sample_table, x_feature, y_feature, reference=reference, comparison=comparison,
                )
                if figure is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        figure, plots_dir, "functional_state_condition_sample_celltype_score_correlation.png", "scatter",
                        "Sample × cell-type functional-score correlation, stratified by condition (exploratory repeated measures; P/FDR omitted)",
                        formats=("png", "svg"), dpi=300, preserve_aspect=True,
                    ))

        if not concordance.empty:
            concordance_labels = {
                feature: heatmap_feature_labels.get(feature, feature)
                for feature in concordance_tf_features
            }
            figure = _regulator_pathway_concordance_figure(
                concordance, concordance_tf_features, concordance_pathway_features,
                tf_labels=concordance_labels,
            )
            if figure is not None:
                result_files.extend(self.save_matplotlib_figure(
                    figure, plots_dir, "functional_state_regulator_pathway_concordance.png", "heatmap",
                    "Regulator–pathway concordance (exploratory repeated-measures aggregate)",
                    formats=("png", "svg"), dpi=300, preserve_aspect=True,
                ))
            figure = _regulator_pathway_condition_concordance_figure(
                condition_concordance, concordance_tf_features, concordance_pathway_features,
                reference=reference, comparison=comparison, tf_labels=concordance_labels,
            )
            if figure is not None:
                result_files.extend(self.save_matplotlib_figure(
                    figure, plots_dir, "functional_state_regulator_pathway_concordance_by_condition.png", "heatmap",
                    "Regulator–pathway concordance, stratified by condition (exploratory repeated measures; P/FDR omitted)",
                    formats=("png", "svg"), dpi=300, preserve_aspect=True,
                ))

        if network_path is None:
            tf_network_manifest = {"status": "not_configured", "requested_source": network_source}
        else:
            try:
                from importlib.metadata import version
                decoupler_version = version("decoupler")
            except Exception:
                decoupler_version = "unknown"
            tf_network_manifest = {
                "source": network_source,
                "file": network_path.name,
                "sha256": _sha256(network_path),
                "method": "decoupler.mt.ulm",
                "score_interpretation": "ULM slope t-statistic; positive=activated, negative=repressed",
                "decoupler_version": decoupler_version,
            }
            if network_source == "project_local":
                tf_network_manifest["path"] = str(network_path.relative_to(Path(self.project_dir).resolve()))
            elif network_resource_metadata:
                tf_network_manifest["resource_metadata"] = network_resource_metadata

        pathway_resource_versions = {}
        if include_builtin_panels:
            pathway_resource_versions["builtin_lipid_inflammation"] = {
                "version": "functional_state_builtin_lipid_inflammation_v2",
                "sha256": _json_sha256(BUILTIN_PATHWAYS),
            }
        if analysis_focus == "ibd_organoid_epithelial":
            pathway_resource_versions["builtin_ibd_epithelial_hypothesis"] = {
                "version": "ibd_organoid_epithelial_hypothesis_v1",
                "sha256": _json_sha256(IBD_EPITHELIAL_QUICK_PATHWAYS),
                "interpretation": (
                    "Focused epithelial differentiation, regenerative/TA and proliferation "
                    "signatures; they are hypothesis scores, not a clinical or universal cell-state ontology."
                ),
            }
        if custom_sets:
            pathway_resource_versions["custom_inline"] = {
                "version": "inline_run_definition",
                "sha256": _json_sha256(custom_sets),
            }
        if gmt_path:
            pathway_resource_versions["project_local_gmt"] = {
                "file": gmt_path.name,
                "sha256": _sha256(gmt_path),
                "selected_terms": list(gmt_sets),
            }
        for details in managed_provenance.values():
            source = details["source"]
            if source in pathway_resource_versions:
                continue
            pathway_resource_versions[source] = {
                key: details[key]
                for key in ("library_key", "resource", "version", "file", "sha256", "license", "license_url", "source_url")
            }

        manifest = {
            "module": self.MODULE_NAME,
            "input_path": os.path.basename(str(input_path)),
            "input_sha256": _sha256(input_path),
            "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
            "expression_source": expression_source,
            "expression_qc": expression_qc,
            "sample_key": sample_key, "condition_key": condition_key, "celltype_key": celltype_key,
            "organism": organism,
            "analysis_focus": analysis_focus,
            "reference_condition": reference, "comparison_condition": comparison,
            "min_cells_per_sample_celltype": min_cells, "min_samples_per_group": min_samples,
            "statistic_method": statistic_method,
            "sample_size_evidence_tiers": {
                "not_runnable": "<2 samples in either condition",
                "exploratory_n2": "2 samples in either condition; displayed as exploratory only",
                "standard_n3_4": "3–4 samples in both conditions",
                "preferred_n5_plus": "≥5 samples in both conditions",
            },
            "tf_coverage_policy": {
                "minimum_detected_targets_for_scoring": 5,
                "below_5": "not scored",
                "5_9": "low coverage; shown with * in the heatmap",
                "10_19": "acceptable",
                "20_plus": "good",
            },
            "pathway_scoring_method": {
                "requested": pathway_scoring_method,
                "per_score_methods": sorted({
                    details["method"] for feature, details in score_provenance.items()
                    if feature_types.get(feature) == "pathway_score"
                }),
            },
            "pathway_gene_sets": {name: list(genes) for name, genes in pathway_sets.items()},
            "builtin_quick_panels": {
                "included": include_builtin_panels,
                "interpretation": "Transparent FAO/inflammation hypothesis signatures; retain for focused exploratory use, but prefer managed Hallmark/Reactome for standard pathway claims.",
            },
            "ibd_epithelial_focus": {
                "enabled": analysis_focus == "ibd_organoid_epithelial",
                "hypothesis_signatures": (
                    {name: list(genes) for name, genes in IBD_EPITHELIAL_QUICK_PATHWAYS.items()}
                    if analysis_focus == "ibd_organoid_epithelial" else {}
                ),
                "managed_hallmark_terms": (
                    list(IBD_EPITHELIAL_MANAGED_TERMS)
                    if analysis_focus == "ibd_organoid_epithelial" else []
                ),
            },
            "gene_set_coverage_policy": {
                "managed_standard_min_detected_genes": 10,
                "good": "detected / original ≥ 70%",
                "warning": "40% ≤ detected / original < 70%; score retained and flagged",
                "unreliable": "detected / original < 40% for a standard set; score not calculated",
                "small_signatures": "curated short or inline exploratory signatures retain their separately labelled 2-gene minimum",
            },
            "pathway_gene_set_overlap_qc": {
                "table": "pathway_gene_set_overlap_qc.csv",
                "policy": "block identical requested or input-detected gene sets; flag detected Jaccard overlap >=0.80 for review",
                "n_pairs": int(len(pathway_overlap_qc)),
                "n_high_overlap_review": int(pathway_overlap_qc["status"].eq("high_detected_overlap_review").sum()),
            },
            "pathway_resource_versions": pathway_resource_versions,
            "gene_expression_panel": {
                "source": "custom" if custom_expression_genes else "builtin_lipid_inflammation_v2",
                "sections": {name: list(genes) for name, genes in expression_panel_sections.items()},
                "detected_genes": list(expression_by_gene),
                "missing_genes": [gene for gene in expression_genes if gene not in expression_by_gene],
            },
            "score_provenance": score_provenance,
            "heatmap_display": {
                "mode": heatmap_selection_mode,
                "features": heat_features,
                "missing_requested_features": missing_heatmap_features,
                "note": "Display selection/order only; full score tables and statistics remain unchanged. Long cell-type labels are rendered as grouped headers above short condition labels.",
            },
            "tf_network": tf_network_manifest,
            "pathway_gmt": ({"path": str(gmt_path.relative_to(Path(self.project_dir))), "sha256": _sha256(gmt_path),
                              "selected_terms": list(gmt_sets)} if gmt_path else {"status": "not_configured"}),
            "managed_gene_sets": ({
                "selected_terms": list(managed_sets),
                "term_provenance": managed_provenance,
            } if managed_sets else {"status": "not_selected"}),
            "correlation": {
                "x_feature_requested": x_requested,
                "y_feature_requested": y_requested,
                "levels": list(correlation_levels),
                "cell_level_interpretation": "descriptive only; cells are not independent biological replicates. Figures show Spearman rho plus an effect-size interpretation and intentionally omit P.",
                "sample_x_celltype_interpretation": "exploratory repeated-measures aggregate; a donor can contribute more than one cell-type point, so P/FDR are not displayed as independent-sample inference.",
                "condition_stratified_views": "For the selected reference and comparison conditions, the module also outputs within-cell-type cell-level and sample × cell-type Spearman effect-size views separately by condition. These views keep pooled results for context, omit P/FDR, and do not remove donor repetition across cell types.",
                "within_celltype_sample_interpretation": "one point per donor within cell type, but bivariate correlation remains exploratory because it does not adjust condition.",
                "confirmatory_model": "InflammationScore ~ TFActivity + Celltype + Condition + (1|Sample), only when biological replication supports stable fitting.",
            },
            "regulator_pathway_concordance": {
                "tf_features_requested": list(_parse_feature_selection(self.params.get("concordance_tf_features", "")) or DEFAULT_CONCORDANCE_TF_FEATURES),
                "pathway_features_requested": list(_parse_feature_selection(self.params.get("concordance_pathway_features", "")) or DEFAULT_CONCORDANCE_PATHWAY_FEATURES),
                "tf_features_used": concordance_tf_features,
                "pathway_features_used": concordance_pathway_features,
                "missing_tf_features": missing_concordance_tfs,
                "missing_pathway_features": missing_concordance_pathways,
                "interpretation": "Sample × cell-type repeated-measures effect-size matrix; do not interpret as independent-sample significance testing.",
                "condition_stratified_table": "regulator_pathway_concordance_by_condition.csv",
                "condition_stratified_interpretation": "Separate reference/comparison matrices retain repeated cell-type units per donor and omit P/FDR; use only as condition-specific effect-size checks.",
            },
            "pseudobulk_validation": {
                "existing_modules": ["sc_pseudobulk_deg", "sc_cell_go"],
                "rule": "Use existing raw-count pseudobulk DEG and complete ranked-statistic GSEA to confirm functional-state hypotheses; this module does not duplicate them.",
            },
        }
        manifest_path = results_dir / "analysis_manifest.json"
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
        result_files.append({"file_path": str(manifest_path), "file_type": "json", "category": "provenance", "label": "Functional-state analysis manifest"})

        self.progress(94, "保存带功能状态列的 AnnData...")
        adata.uns["functional_state"] = manifest
        output_path = self.save_output(adata, self.MODULE_NAME)
        warnings = []
        if expression_qc.get("warning"):
            warnings.append(expression_qc["warning"])
        if network_path is None:
            if network_source == "managed_collectri":
                warnings.append("平台管理的 CollecTRI 快照不可用：未计算 TF activity；pathway score 与 gene expression 仍已完成。")
            else:
                warnings.append("未启用 TF 网络：未计算 TF activity；pathway score 与 gene expression 仍已完成。")
        if any(row["status"] == "small_signature" for row in coverage_rows):
            warnings.append("部分内置小型基因集只有 3 个基因；请结合完整 Hallmark/Reactome 基因集复核。")
        if any(row["status"] == "exploratory_custom_signature" for row in coverage_rows):
            warnings.append("部分自定义基因集少于一般通路的最少基因数，仅作为探索性小型 signature 评分。")
        if missing_heatmap_features:
            warnings.append(
                "以下热图展示特征不可用，已跳过：" + "、".join(missing_heatmap_features[:8])
            )
        if missing_concordance_tfs or missing_concordance_pathways:
            missing = (missing_concordance_tfs + missing_concordance_pathways)[:8]
            warnings.append("一致性矩阵中以下特征不可用，已跳过：" + "、".join(missing))
        n_low_coverage_tf = int(tf_coverage.get("activity_status", pd.Series(dtype=str)).eq("low_coverage").sum())
        if n_low_coverage_tf:
            warnings.append(f"{n_low_coverage_tf} 个 TF 仅有 5–9 个检测靶基因；热图以 * 标记，须谨慎解释。")
        self.progress(100, "功能状态分析完成。")
        return {
            "output_adata": output_path, "result_files": result_files,
            "summary": {
                "n_cells": int(adata.n_obs), "n_pathway_scores": int(sum(value == "pathway_score" for value in feature_types.values())),
                "n_tf_activities": int(sum(value == "tf_activity" for value in feature_types.values())),
                "n_expression_features": int(sum(value == "gene_expression" for value in feature_types.values())),
                "reference_condition": reference, "comparison_condition": comparison,
                "sample_level_statistics": True, "tf_network_configured": bool(network_path),
                "tf_network_source": network_source,
                "n_heatmap_display_features": len(heat_features),
                "n_regulator_pathway_concordance_pairs": int(len(concordance)),
                "n_condition_stratified_regulator_pathway_concordance_pairs": int(len(condition_concordance)),
                "expression_source": expression_source,
                "correlation_levels": list(correlation_levels),
                "warnings": warnings, "scope_key": str(self.params.get("scope_key", "") or "") or None,
            },
        }
