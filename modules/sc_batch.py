"""Manifest-driven import helpers for many server-resident 10x datasets.

The platform deliberately separates *where data are stored* from *what the
experimental design is*.  Directory names are useful for discovery, but never
become biological conditions or replicates automatically.  A user-reviewed
manifest supplies that information before data are merged.
"""

import os
import re
from itertools import combinations

import numpy as np
import pandas as pd

from config import Config
from modules.base import BaseAnalysis
from modules.io_utils import read_10x_mtx_compat, summarize_adata_import


MANIFEST_REQUIRED_COLUMNS = ("sample_id", "matrix_dir", "condition", "replicate")
TENX_MATRIX_FILENAMES = ("matrix.mtx", "matrix.mtx.gz")
TENX_BARCODE_FILENAMES = ("barcodes.tsv", "barcodes.tsv.gz")
TENX_FEATURE_FILENAMES = (
    "features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz",
)


def _safe_name(value, fallback="batch"):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return name or fallback


def _nonempty_text(value):
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def is_10x_matrix_dir(path):
    """Return whether ``path`` contains one complete 10x Matrix Market triplet."""
    if not os.path.isdir(path):
        return False
    names = set(os.listdir(path))
    return (
        bool(names.intersection(TENX_MATRIX_FILENAMES))
        and bool(names.intersection(TENX_BARCODE_FILENAMES))
        and bool(names.intersection(TENX_FEATURE_FILENAMES))
    )


def find_10x_matrix_dirs(source_root):
    """Recursively find complete 10x directories below a validated source root."""
    root = os.path.realpath(source_root)
    if not os.path.isdir(root):
        raise FileNotFoundError(f"10x 源目录不存在: {source_root}")
    found = []
    for current, _dirs, files in os.walk(root):
        names = set(files)
        if (
            bool(names.intersection(TENX_MATRIX_FILENAMES))
            and bool(names.intersection(TENX_BARCODE_FILENAMES))
            and bool(names.intersection(TENX_FEATURE_FILENAMES))
        ):
            found.append(os.path.realpath(current))
    return sorted(found)


def discovered_manifest_frame(source_root):
    """Build an editable manifest template without guessing biological labels."""
    root = os.path.realpath(source_root)
    rows = []
    used_ids = set()
    for index, matrix_dir in enumerate(find_10x_matrix_dirs(root), start=1):
        relative = os.path.relpath(matrix_dir, root)
        base_id = _safe_name(os.path.basename(matrix_dir), f"sample_{index}")
        sample_id = base_id
        suffix = 2
        while sample_id in used_ids:
            sample_id = f"{base_id}_{suffix}"
            suffix += 1
        used_ids.add(sample_id)
        rows.append({
            "sample_id": sample_id,
            "matrix_dir": relative,
            "condition": "",
            "replicate": "",
            "batch": "",
            "notes": "请填写 condition 和 replicate；目录名不自动代表生物学分组。",
        })
    return pd.DataFrame(rows, columns=[
        "sample_id", "matrix_dir", "condition", "replicate", "batch", "notes",
    ])


def _read_manifest(path):
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(f"样本 manifest 不存在: {path}")
    lower = str(path).lower()
    if lower.endswith((".tsv", ".txt")):
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    elif lower.endswith((".xlsx", ".xls")):
        frame = pd.read_excel(path, dtype=str).fillna("")
    else:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def read_sample_manifest(manifest_path, source_root):
    """Validate a sample manifest and resolve every 10x directory safely.

    ``matrix_dir`` is normally relative to ``source_root``.  An absolute path
    is permitted only when it is still below the configured server source root.
    ``condition`` and ``replicate`` must be explicit: without them, pseudobulk
    inference cannot distinguish biology from an arbitrary directory layout.
    """
    root = Config.validate_sc_batch_source_path(source_root)
    frame = _read_manifest(manifest_path)
    missing = [name for name in MANIFEST_REQUIRED_COLUMNS if name not in frame.columns]
    if missing:
        raise ValueError("样本 manifest 缺少必需列: " + ", ".join(missing))
    if frame.empty:
        raise ValueError("样本 manifest 没有任何样本")

    frame = frame.copy()
    for column in frame.columns:
        frame[column] = frame[column].map(_nonempty_text)
    for column in MANIFEST_REQUIRED_COLUMNS:
        invalid = frame[column].eq("")
        if bool(invalid.any()):
            rows = ", ".join(str(index + 2) for index in frame.index[invalid][:5])
            raise ValueError(f"样本 manifest 的 '{column}' 存在空值（第 {rows} 行）")
    if frame["sample_id"].duplicated().any():
        duplicates = frame.loc[frame["sample_id"].duplicated(keep=False), "sample_id"].tolist()
        raise ValueError("sample_id 必须唯一；重复值: " + ", ".join(map(str, duplicates[:8])))
    if frame["sample_id"].map(lambda item: _safe_name(item) != item).any():
        raise ValueError("sample_id 仅允许字母、数字、点、下划线和连字符")

    resolved = []
    for row_index, row in frame.iterrows():
        supplied = row["matrix_dir"]
        candidate = supplied if os.path.isabs(supplied) else os.path.join(root, supplied)
        try:
            matrix_dir = Config.validate_sc_batch_source_path(candidate)
        except ValueError as exc:
            raise ValueError(
                f"manifest 第 {row_index + 2} 行 sample_id={row['sample_id']} 的 matrix_dir 无效: {exc}"
            ) from exc
        if not (matrix_dir == root or matrix_dir.startswith(root + os.sep)):
            raise ValueError(f"matrix_dir 必须位于 source_root 下: {supplied}")
        if not is_10x_matrix_dir(matrix_dir):
            raise ValueError(
                f"sample_id={row['sample_id']} 的目录不是完整 10x 矩阵: {supplied}"
            )
        resolved.append(matrix_dir)
    frame["matrix_path"] = resolved
    if "batch" not in frame.columns:
        frame["batch"] = ""
    return frame


def sample_design_from_obs(adata, sample_key="sample_id", condition_key="condition"):
    """Return one validated design row per biological sample.

    A sample occupying more than one condition is a design error, rather than
    something that should be silently split into pseudo-replicates.
    """
    for column in (sample_key, condition_key):
        if not column or column not in adata.obs.columns:
            raise ValueError(f"adata.obs 中缺少 '{column}' 列")
    design = pd.DataFrame({
        "sample_id": adata.obs[sample_key].map(_nonempty_text),
        "condition": adata.obs[condition_key].map(_nonempty_text),
    }, index=adata.obs_names)
    if bool(design["sample_id"].eq("").any()) or bool(design["condition"].eq("").any()):
        raise ValueError("sample_id 或 condition 存在空值，不能进行样本级导出/比较")
    per_sample = design.groupby("sample_id", observed=True)["condition"].nunique()
    mixed = per_sample[per_sample > 1]
    if not mixed.empty:
        raise ValueError("同一 sample_id 对应多个 condition: " + ", ".join(mixed.index[:8]))
    return design.drop_duplicates("sample_id").sort_values("sample_id").reset_index(drop=True)


def build_comparison_plan(design, requested_comparisons="", mode="all_pairwise",
                          reference_group="", min_samples_per_group=2):
    """Create an explicit, reproducible condition comparison registry."""
    groups = sorted(design["condition"].dropna().astype(str).unique().tolist())
    if len(groups) < 2:
        return pd.DataFrame(columns=[
            "comparison_id", "comparison", "group_1", "group_2",
            "n_samples_group_1", "n_samples_group_2", "status", "reason",
        ])

    pairs = []
    manual = str(requested_comparisons or "").strip()
    if manual:
        for raw in re.split(r"[;\n]+", manual):
            item = raw.strip()
            if not item:
                continue
            parts = re.split(r"\s*(?:-vs-|\s+vs\s+)\s*", item, maxsplit=1,
                             flags=re.IGNORECASE)
            if len(parts) != 2 or not all(part.strip() for part in parts):
                raise ValueError(f"比较格式无效: '{item}'；请使用 GroupA-vs-GroupB")
            pairs.append((parts[0].strip(), parts[1].strip()))
    elif mode == "vs_reference":
        reference = str(reference_group or "").strip()
        if reference not in groups:
            raise ValueError("vs_reference 模式需要选择一个存在的 reference_group")
        pairs = [(group, reference) for group in groups if group != reference]
    else:
        pairs = list(combinations(groups, 2))

    rows = []
    seen = set()
    for group_1, group_2 in pairs:
        if group_1 == group_2 or (group_1, group_2) in seen:
            continue
        seen.add((group_1, group_2))
        count_1 = int((design["condition"] == group_1).sum())
        count_2 = int((design["condition"] == group_2).sum())
        valid_groups = group_1 in groups and group_2 in groups
        enough = valid_groups and count_1 >= int(min_samples_per_group) and count_2 >= int(min_samples_per_group)
        if not valid_groups:
            reason = "比较组不在 manifest condition 中"
        elif not enough:
            reason = f"每组至少需要 {int(min_samples_per_group)} 个独立样本"
        else:
            reason = ""
        label = f"{group_1} vs {group_2}"
        rows.append({
            "comparison_id": f"{_safe_name(group_1)}_vs_{_safe_name(group_2)}",
            "comparison": label,
            "group_1": group_1,
            "group_2": group_2,
            "n_samples_group_1": count_1,
            "n_samples_group_2": count_2,
            "status": "runnable" if enough else "not_runnable",
            "reason": reason,
        })
    return pd.DataFrame(rows)


def aggregate_pseudobulk_counts(adata, sample_key="sample_id", condition_key="condition",
                                celltype_key="", min_cells=1):
    """Aggregate a count layer as sample × optional cell type × gene records."""
    from scipy import sparse

    if "counts" not in adata.layers:
        raise ValueError("缺少 layers['counts']；无法生成原始 counts pseudobulk")
    design = sample_design_from_obs(adata, sample_key, condition_key)
    matrix = adata.layers["counts"]
    if sparse.issparse(matrix):
        values = np.asarray(matrix.data, dtype=float)
    else:
        values = np.asarray(matrix, dtype=float).reshape(-1)
    if values.size and (np.nanmin(values) < 0 or not np.allclose(values[:min(values.size, 100_000)], np.round(values[:min(values.size, 100_000)]), atol=1e-6)):
        raise ValueError("counts 层不是非负整数原始计数，不能用于 pseudobulk 差异分析")

    work = pd.DataFrame({
        "sample_id": adata.obs[sample_key].map(_nonempty_text),
        "condition": adata.obs[condition_key].map(_nonempty_text),
        "row_index": np.arange(adata.n_obs),
    }, index=adata.obs_names)
    requested_celltype = str(celltype_key or "").strip()
    if requested_celltype:
        if requested_celltype not in adata.obs.columns:
            raise ValueError(f"adata.obs 中没有 celltype_key '{requested_celltype}'")
        work["celltype"] = adata.obs[requested_celltype].map(_nonempty_text)
        if bool(work["celltype"].eq("").any()):
            raise ValueError(f"celltype_key '{requested_celltype}' 存在空值")
        group_columns = ["sample_id", "condition", "celltype"]
    else:
        work["celltype"] = "All"
        group_columns = ["sample_id", "condition", "celltype"]

    rows = []
    count_rows = []
    for keys, positions in work.groupby(group_columns, observed=True)["row_index"]:
        indices = np.asarray(positions, dtype=int)
        n_cells = int(len(indices))
        if n_cells < int(min_cells):
            continue
        vector = np.asarray(matrix[indices, :].sum(axis=0)).reshape(-1)
        rows.append({
            "sample_id": str(keys[0]), "condition": str(keys[1]),
            "celltype": str(keys[2]), "n_cells": n_cells,
            "library_size": float(vector.sum()),
        })
        count_rows.append(vector)
    if not rows:
        return pd.DataFrame(columns=["sample_id", "condition", "celltype", "n_cells", "library_size"]), np.empty((0, adata.n_vars))
    return pd.DataFrame(rows), np.vstack(count_rows)


def _available_path(directory, stem, suffix):
    os.makedirs(directory, exist_ok=True)
    candidate = os.path.join(directory, f"{stem}{suffix}")
    index = 2
    while os.path.exists(candidate):
        candidate = os.path.join(directory, f"{stem}_{index}{suffix}")
        index += 1
    return candidate


class SCBatchImport(BaseAnalysis):
    """Merge arbitrary many 10x samples using an explicit biological manifest."""

    MODULE_NAME = "sc_batch_import"
    DISPLAY_NAME = "批量 10x manifest 导入"
    DESCRIPTION = "读取多个服务器 10x 矩阵并按样本 manifest 合并为标准 h5ad"

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import anndata as ad

        manifest_path = str(self.params.get("manifest_path") or "").strip()
        source_root = str(self.params.get("source_root") or "").strip()
        dataset_name = _safe_name(self.params.get("dataset_name"), "sc_batch")
        species = self.params.get("species")
        genome = self.params.get("genome")

        if not manifest_path:
            raise ValueError("缺少 manifest_path；请先下载并填写样本 manifest")
        project_root = os.path.realpath(self.project_dir)
        manifest_real = os.path.realpath(manifest_path)
        if not manifest_real.startswith(project_root + os.sep):
            raise ValueError("样本 manifest 必须位于项目目录内")
        manifest = read_sample_manifest(manifest_real, source_root)

        self.progress(5, f"已验证 {len(manifest)} 个 10x 样本及其生物学设计")
        imported = []
        import_rows = []
        total = len(manifest)
        for index, (_, row) in enumerate(manifest.iterrows(), start=1):
            sample_id = row["sample_id"]
            self.progress(5 + int((index - 1) * 65 / total),
                          f"正在读取 {index}/{total}: {sample_id}")
            adata = read_10x_mtx_compat(row["matrix_path"], var_names="gene_symbols")
            adata.var_names_make_unique()
            adata.obs_names = pd.Index([f"{sample_id}:{barcode}" for barcode in adata.obs_names.astype(str)])
            if "gene_ids" in adata.var.columns:
                adata.var["gene_id"] = adata.var["gene_ids"].astype(str)
            else:
                adata.var["gene_id"] = adata.var_names.astype(str)
            adata.var["gene_name"] = adata.var_names.astype(str)
            adata.layers["counts"] = adata.X.copy()
            for column in manifest.columns:
                if column in {"matrix_path", "matrix_dir", "notes"}:
                    continue
                adata.obs[column] = row[column]
            adata.obs["source_matrix_dir"] = row["matrix_dir"]
            imported.append(adata)
            totals = np.asarray(adata.X.sum(axis=1)).reshape(-1)
            import_rows.append({
                "sample_id": sample_id,
                "condition": row["condition"],
                "replicate": row["replicate"],
                "batch": row.get("batch", ""),
                "matrix_dir": row["matrix_dir"],
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "median_umi": float(np.median(totals)) if totals.size else 0.0,
            })

        self.progress(72, "正在按共有基因/并集基因合并所有样本...")
        combined = ad.concat(
            imported, axis=0, join="outer", merge="first", index_unique=None,
            fill_value=0,
        )
        combined.obs["sample_id"] = combined.obs["sample_id"].astype(str)
        combined.obs["condition"] = combined.obs["condition"].astype(str)
        combined.obs["replicate"] = combined.obs["replicate"].astype(str)
        combined.uns["sc_batch_import"] = {
            "manifest_file": os.path.basename(manifest_real),
            "source_root": Config.validate_sc_batch_source_path(source_root),
            "n_samples": int(len(manifest)),
            "sample_key": "sample_id",
            "condition_key": "condition",
            "replicate_key": "replicate",
            "dataset_name": dataset_name,
        }
        if species:
            combined.uns["species"] = species
        if genome:
            combined.uns["genome"] = genome
        combined.uns["input_format"] = "10x_mtx_manifest_batch"

        uploads_dir = os.path.join(self.project_dir, "uploads")
        results_dir = os.path.join(self.project_dir, "results")
        output_path = _available_path(uploads_dir, f"{dataset_name}_imported", ".h5ad")
        summary_path = _available_path(results_dir, f"{dataset_name}_sample_import_summary", ".csv")
        resolved_path = _available_path(results_dir, f"{dataset_name}_resolved_manifest", ".csv")
        result_files = []
        self.progress(85, "正在保存合并 h5ad 和样本设计表...")
        combined.write_h5ad(output_path)
        result_files.append({
            "file_path": output_path, "file_type": "h5ad", "category": "data",
            "label": "Manifest-merged single-cell h5ad",
        })
        pd.DataFrame(import_rows).to_csv(summary_path, index=False)
        result_files.append({
            "file_path": summary_path, "file_type": "csv", "category": "table",
            "label": "Batch 10x sample import summary",
        })
        manifest.drop(columns=["matrix_path"], errors="ignore").to_csv(resolved_path, index=False)
        result_files.append({
            "file_path": resolved_path, "file_type": "csv", "category": "table",
            "label": "Resolved biological sample manifest",
        })

        summary = summarize_adata_import(combined, "10x_mtx_manifest_batch", output_path)
        summary.update({
            "n_samples": int(len(manifest)),
            "sample_key": "sample_id", "condition_key": "condition",
            "conditions": sorted(manifest["condition"].astype(str).unique().tolist()),
            "sample_counts_by_condition": manifest["condition"].value_counts().sort_index().to_dict(),
        })
        self.progress(100, f"批量导入完成：{summary['n_samples']} 个样本，{summary['n_cells']} 个细胞")
        return {
            "output_adata": output_path,
            "result_files": result_files,
            "summary": summary,
        }
