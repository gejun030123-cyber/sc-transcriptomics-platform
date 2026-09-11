"""Shared statistical contracts for single-cell differential expression.

Cell-level marker/condition tests need a non-negative log-expression scale to
make Scanpy's reported ``logfoldchanges`` interpretable.  Pearson residuals
are useful for PCA and neighbourhood graphs, but their signed residual scale
is not a fold-change scale.  This module therefore always reconstructs a
temporary log1p matrix from the protected raw-count layer for cell-level DE.
"""

from __future__ import annotations


def resolve_cell_grouping(params, obs_columns):
    """Resolve the biological grouping used for single-cell DE.

    ``cluster`` remains available for workflows that deliberately operate on
    unsupervised clusters.  The analysis form, however, now selects
    ``annotated_celltype`` by default so post-annotation condition contrasts
    are grouped by the reviewed cell-type labels rather than by Leiden IDs.

    The ``auto`` mode is intentionally retained for older API/pipeline calls:
    it prefers an explicitly supplied key and otherwise uses ``celltype``
    when available before falling back to ``leiden``.  All modes return the
    legacy-compatible ``cluster`` output field via the caller, while recording
    the actual grouping semantics separately in the audit/result metadata.
    """
    params = dict(params or {})
    columns = {str(column) for column in obs_columns}
    raw_mode = str(
        params.get("grouping_mode", params.get("analysis_unit", "auto")) or "auto"
    ).strip().lower()
    aliases = {
        "celltype": "annotated_celltype",
        "cell_type": "annotated_celltype",
        "annotation": "annotated_celltype",
        "annotated": "annotated_celltype",
        "annotated_celltype": "annotated_celltype",
        "cluster": "cluster",
        "clusters": "cluster",
        "leiden": "cluster",
        "auto": "auto",
    }
    if raw_mode not in aliases:
        raise ValueError("grouping_mode 必须为 annotated_celltype、cluster 或 auto")
    mode = aliases[raw_mode]
    celltype_key = str(params.get("celltype_key", "celltype") or "").strip()
    cluster_key = str(params.get("cluster_key", "leiden") or "").strip()

    if mode == "annotated_celltype":
        key = celltype_key
        label = "细胞类型"
    elif mode == "cluster":
        key = cluster_key
        label = "聚类"
    else:
        # Preserve scripts that explicitly requested a legacy cluster key,
        # while making a bare call annotation-aware when a celltype column is
        # already present.
        explicitly_requested_celltype = "celltype_key" in params and bool(celltype_key)
        explicitly_requested_cluster = "cluster_key" in params and bool(cluster_key)
        if explicitly_requested_celltype and celltype_key in columns:
            key, mode, label = celltype_key, "annotated_celltype", "细胞类型"
        elif explicitly_requested_cluster and cluster_key in columns:
            key, mode, label = cluster_key, "cluster", "聚类"
        elif celltype_key in columns:
            key, mode, label = celltype_key, "annotated_celltype", "细胞类型"
        else:
            key, mode, label = cluster_key, "cluster", "聚类"

    if not key:
        raise ValueError(f"未指定{label}分组列")
    if key not in columns:
        if mode == "annotated_celltype":
            raise ValueError(
                f"缺少注释细胞类型列 '{key}'；请先完成细胞注释，"
                "或将“分析分组”改为聚类。"
            )
        raise ValueError(f"缺少聚类列 '{key}'")
    return {
        "mode": mode,
        "key": key,
        "label": label,
        "value_label": "细胞类型" if mode == "annotated_celltype" else "cluster",
    }


def _matrix_is_raw_counts(matrix):
    """Return whether *matrix* is finite, non-negative and integer-valued."""
    import numpy as np

    if hasattr(matrix, "data") and hasattr(matrix, "tocsr"):
        values = np.asarray(matrix.data, dtype=float)
    else:
        values = np.asarray(matrix, dtype=float).reshape(-1)
    if not values.size:
        return True
    return bool(
        np.isfinite(values).all()
        and (values >= 0).all()
        and np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6)
    )


def log1p_expression_from_counts(matrix, target_sum=10_000.0):
    """Return library-size-normalized log1p expression from raw counts.

    The input matrix is never mutated.  Sparse inputs remain sparse, which is
    important for annotation and DEG runs on ordinary single-cell datasets.
    """
    import numpy as np
    from scipy import sparse

    if not _matrix_is_raw_counts(matrix):
        raise ValueError(
            "表达矩阵不是非负整数原始计数；不能重建 library-size normalized log1p。"
        )
    target_sum = float(target_sum)
    if not np.isfinite(target_sum) or target_sum <= 0:
        raise ValueError("target_sum 必须是大于 0 的有限数值。")

    if sparse.issparse(matrix):
        work = matrix.astype(float).tocsr(copy=True)
        library_size = np.asarray(work.sum(axis=1), dtype=float).reshape(-1)
        scale = np.divide(
            target_sum,
            library_size,
            out=np.zeros_like(library_size, dtype=float),
            where=library_size > 0,
        )
        work = work.multiply(scale[:, None]).tocsr()
        work.data = np.log1p(work.data)
        work.eliminate_zeros()
        return work

    work = np.asarray(matrix, dtype=float).copy()
    library_size = work.sum(axis=1)
    scale = np.divide(
        target_sum,
        library_size,
        out=np.zeros_like(library_size, dtype=float),
        where=library_size > 0,
    )
    return np.log1p(work * scale[:, None])


def log1p_adata_for_cell_level_de(adata, target_sum=10_000.0):
    """Return a temporary AnnData with log1p-normalized raw counts in ``X``.

    The input object is never mutated.  Reconstructing the scale from
    ``layers['counts']`` makes cell-level marker and exploratory comparisons
    safe after either log1p normalization or Pearson-residual normalization.
    """
    if "counts" not in adata.layers:
        normalization = adata.uns.get("normalization", {}) or {}
        if normalization.get("x_contains") == "pearson_residuals":
            raise ValueError(
                "当前 X 是 Pearson residuals 且缺少 layers['counts']；"
                "不能计算有可解释 log2FC 的细胞级 DEG。请从保留原始 counts 的 AnnData 重新开始。"
            )
        raise ValueError(
            "细胞级 DEG 需要 layers['counts'] 中的原始 UMI counts，以确保 log2FC 可解释。"
        )
    if not _matrix_is_raw_counts(adata.layers["counts"]):
        raise ValueError(
            "layers['counts'] 不是非负整数原始计数；不能据此计算可解释的细胞级 log2FC。"
        )

    work = adata.copy()
    work.X = log1p_expression_from_counts(
        work.layers["counts"], target_sum=float(target_sum),
    )
    work.uns["log1p"] = {"base": None}
    work.uns["cell_level_de_expression"] = {
        "source": "layers[counts]",
        "transform": "library_size_normalize_then_log1p",
        "target_sum": float(target_sum),
    }
    return work


def normalization_semantics(adata):
    """Return a concise, serializable source-scale description for audit logs."""
    normalization = adata.uns.get("normalization", {}) or {}
    return {
        "input_x_contains": normalization.get("x_contains", "unknown"),
        "cell_level_de_expression": "log1p_from_layers[counts]",
    }
