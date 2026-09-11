"""SCENIC regulon-activity result analysis.

This module intentionally analyses *precomputed* AUCell matrices stored inside
an AnnData object.  SCENIC network inference requires versioned cisTarget
databases and sizeable local compute resources; silently rebuilding a network
from a log-normalised ``X`` matrix would not be reproducible.  Keeping the
AUCell matrix and its regulon definitions in the project H5AD makes the
visualisation step lightweight, local, and auditable.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from modules.base import BaseAnalysis


DEFAULT_AUC_KEYS = ("X_aucell", "aucell", "X_scenic_auc", "scenic_auc")
DEFAULT_REGULON_NAME_KEYS = (
    "scenic_regulon_names", "regulon_names", "aucell_regulon_names",
)
DEFAULT_TARGET_KEYS = ("scenic_regulon_targets", "regulon_targets")


def _safe_name(value: object, default: str = "manual_run") -> str:
    # ``\w`` retains Unicode letters under Python's default regex behaviour,
    # so a user-facing Chinese analysis ID remains a distinct safe directory
    # component rather than collapsing to the shared manual_run fallback.
    text = re.sub(r"[^\w.-]+", "_", str(value or "").strip())
    return text.strip("._") or default


def _top_indices(values, limit: int) -> np.ndarray:
    values = np.nan_to_num(np.asarray(values, dtype=float), nan=-np.inf)
    limit = max(1, min(int(limit), values.size))
    # Stable mergesort keeps the input regulon order for ties, which makes the
    # exported figure selection reproducible.
    return np.argsort(-values, kind="mergesort")[:limit]


def _row_zscore(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    mean = np.nanmean(matrix, axis=1, keepdims=True)
    std = np.nanstd(matrix, axis=1, keepdims=True)
    std[~np.isfinite(std) | (std < 1e-12)] = 1.0
    return np.clip((matrix - mean) / std, -3.0, 3.0)


def _regulon_tf(regulon: str) -> str:
    """Return a readable TF label from common pySCENIC regulon names."""
    return re.split(r"\s*\(|\s*\[|_extended$", str(regulon), maxsplit=1)[0].strip() or str(regulon)


def _positive_int(value, *, default: int, minimum: int, label: str) -> int:
    """Read a required positive integer without treating an explicit 0 as unset."""
    raw = default if value is None else value
    try:
        parsed = int(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} 必须是整数。") from exc
    if parsed < minimum:
        raise ValueError(f"{label} 必须不小于 {minimum}。")
    return parsed


class ScenicAnalysis(BaseAnalysis):
    """Create auditable SCENIC/AUCell regulon activity result panels."""

    MODULE_NAME = "scenic"
    DISPLAY_NAME = "SCENIC 调控子活性"
    DESCRIPTION = "本地预计算 AUCell 的 regulon 活性、RSS 特异性、网络和共活性模块图"
    # ``groupby`` defaults to celltype but can intentionally target a reviewed
    # cluster/annotation column.  Require the selected column at runtime
    # instead of rejecting otherwise valid AUCell inputs before the user can
    # choose that grouping.
    INPUT_REQUIRES = []

    def _resolve_auc_matrix(self, adata):
        requested = str(self.params.get("auc_obsm_key", "") or "").strip()
        keys = [requested] if requested else []
        keys.extend(key for key in DEFAULT_AUC_KEYS if key not in keys)
        selected_key = next((key for key in keys if key and key in adata.obsm), None)
        if not selected_key:
            raise ValueError(
                "未找到预计算的 AUCell 矩阵。请将每细胞 × regulon 的矩阵写入 "
                "adata.obsm['X_aucell']（推荐，DataFrame 列名为 regulon 名），"
                "或在本页指定实际 obsm 键；不能用 TF 基因表达替代 regulon 活性。"
            )

        raw = adata.obsm[selected_key]
        if isinstance(raw, pd.DataFrame):
            names = [str(item) for item in raw.columns]
            matrix = raw.to_numpy(dtype=float, copy=True)
        else:
            matrix = np.asarray(raw, dtype=float)
            requested_names_key = str(self.params.get("regulon_names_uns_key", "") or "").strip()
            name_keys = [requested_names_key] if requested_names_key else []
            name_keys.extend(key for key in DEFAULT_REGULON_NAME_KEYS if key not in name_keys)
            names = None
            for key in name_keys:
                candidate = adata.uns.get(key)
                if isinstance(candidate, (list, tuple, np.ndarray, pd.Index)):
                    names = [str(item) for item in candidate]
                    break
                if isinstance(candidate, dict):
                    for nested_key in ("regulon_names", "columns", "names"):
                        nested = candidate.get(nested_key)
                        if isinstance(nested, (list, tuple, np.ndarray, pd.Index)):
                            names = [str(item) for item in nested]
                            break
                    if names is not None:
                        break
            if names is None:
                raise ValueError(
                    "AUCell obsm 矩阵没有 regulon 列名。请使用 pandas DataFrame 写入 "
                    "obsm，或在 adata.uns['scenic_regulon_names'] 写入与列数一致的名称列表。"
                )

        if matrix.ndim != 2 or matrix.shape[0] != adata.n_obs or matrix.shape[1] < 1:
            raise ValueError("AUCell 矩阵必须是与 AnnData 细胞数一致的二维 cell × regulon 矩阵。")
        if len(names) != matrix.shape[1]:
            raise ValueError("regulon 名称数与 AUCell 矩阵列数不一致。")
        if len(set(names)) != len(names):
            raise ValueError("regulon 名称存在重复；请在写入 AUCell 矩阵前使用唯一列名。")
        if not np.isfinite(matrix).all():
            raise ValueError("AUCell 矩阵含 NaN 或无穷值；请先在 SCENIC 上游结果中处理。")
        return matrix, names, selected_key

    def _resolve_group(self, adata):
        group_key = str(self.params.get("groupby", "celltype") or "").strip()
        if not group_key or group_key not in adata.obs.columns:
            raise ValueError(f"细胞分组列 '{group_key or '（空）'}' 不在 adata.obs 中。")
        values = adata.obs[group_key]
        if pd.api.types.is_numeric_dtype(values) and not isinstance(values.dtype, pd.CategoricalDtype):
            raise ValueError("RSS 必须按分类细胞类型/cluster 计算；groupby 不能是连续数值列。")
        if values.isna().any():
            raise ValueError(f"分组列 '{group_key}' 含缺失值；请先完成或清理细胞注释。")
        labels = values.astype(str).to_numpy()
        if isinstance(values.dtype, pd.CategoricalDtype):
            groups = [str(item) for item in values.cat.categories if str(item) in set(labels)]
        else:
            groups = list(dict.fromkeys(labels.tolist()))
        if len(groups) < 2:
            raise ValueError("RSS 至少需要两个细胞类型/cluster。")
        return group_key, labels, groups

    def _resolve_targets(self, adata):
        requested = str(self.params.get("regulon_targets_uns_key", "") or "").strip()
        keys = [requested] if requested else []
        keys.extend(key for key in DEFAULT_TARGET_KEYS if key not in keys)
        for key in keys:
            raw = adata.uns.get(key)
            if not isinstance(raw, dict):
                continue
            parsed = {}
            for regulon, targets in raw.items():
                if isinstance(targets, dict):
                    targets = targets.get("targets", targets.get("genes", []))
                if isinstance(targets, (list, tuple, np.ndarray, pd.Index, set)):
                    parsed[str(regulon)] = [str(gene) for gene in targets if str(gene).strip()]
            if parsed:
                return parsed, key
        return {}, ""

    @staticmethod
    def _heatmap_cells(labels, groups, max_cells):
        max_cells = max(1, int(max_cells))
        per_group = {group: np.flatnonzero(labels == group) for group in groups}
        selected = []
        selected_set = set()
        # First allocate an equal deterministic quota, then fill remaining
        # positions in group order.  This prevents a large cell type from
        # visually erasing rare groups in a cell-level heatmap.
        base = max_cells // len(groups)
        for group in groups:
            indices = per_group[group]
            take = min(len(indices), base)
            if take:
                chosen = indices[np.linspace(0, len(indices) - 1, take, dtype=int)].tolist()
                selected.extend(chosen)
                selected_set.update(chosen)
        remaining = max_cells - len(selected)
        if remaining > 0:
            for group in groups:
                if remaining <= 0:
                    break
                indices = per_group[group]
                available = np.asarray([item for item in indices if item not in selected_set], dtype=int)
                take = min(len(available), remaining)
                if take:
                    chosen = available[np.linspace(0, len(available) - 1, take, dtype=int)]
                    selected.extend(chosen.tolist())
                    selected_set.update(chosen.tolist())
                    remaining -= take
        return np.asarray(selected, dtype=int)

    @staticmethod
    def _rss(matrix, labels, groups, chunk_size=8):
        """Regulon specificity score = 1 - Jensen-Shannon distance.

        This mirrors the SCENIC RSS interpretation: a regulon concentrated in
        one group approaches one, whereas uniform activity approaches zero.
        AUCell values are non-negative by construction; clipped values make a
        defensive interpretation possible for transformed imported matrices.
        """
        values = np.clip(np.asarray(matrix, dtype=float), 0.0, None)
        column_sums = values.sum(axis=0)
        normalizer = np.where(column_sums > 0, column_sums, 1.0)
        n_regulons = values.shape[1]
        chunk_size = max(1, int(chunk_size))
        scores = np.zeros((len(groups), n_regulons), dtype=float)
        # Vectorise across a small regulon block, not every group at once.
        # The latter would allocate a cells × regulons × groups tensor that is
        # unsafe for typical 50k-cell SCENIC runs.
        for group_index, group in enumerate(groups):
            mask = labels == group
            q = mask.astype(float)
            q /= q.sum()
            q_positive = q > 0
            for start in range(0, n_regulons, chunk_size):
                end = min(start + chunk_size, n_regulons)
                p = values[:, start:end] / normalizer[None, start:end]
                midpoint = 0.5 * (p + q[:, None])
                with np.errstate(divide="ignore", invalid="ignore"):
                    kl_p = np.sum(np.where(
                        p > 0,
                        p * (np.log2(p) - np.log2(midpoint)),
                        0.0,
                    ), axis=0)
                    # q is sparse for a cell-group indicator, so calculate
                    # its term only on in-group cells instead of broadcasting
                    # an additional full cells × regulons array.
                    q_values = q[q_positive, None]
                    kl_q = np.sum(
                        q_values * (np.log2(q_values) - np.log2(midpoint[q_positive, :])),
                        axis=0,
                    )
                js_distance = np.sqrt(np.clip(0.5 * (kl_p + kl_q), 0.0, 1.0))
                scores[group_index, start:end] = 1.0 - js_distance
        scores[:, column_sums <= 0] = 0.0
        np.clip(scores, 0.0, 1.0, out=scores)
        return scores

    def _build_aucell_heatmap(self, matrix, regulons, labels, groups, selected_regulons, max_cells):
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT

        cell_indices = self._heatmap_cells(labels, groups, max_cells)
        ordered = np.concatenate([
            cell_indices[labels[cell_indices] == group] for group in groups
        ])
        shown = matrix[ordered][:, selected_regulons].T
        shown_z = _row_zscore(shown)
        shown_labels = labels[ordered]
        group_codes = np.array([groups.index(item) for item in shown_labels], dtype=int)

        fig = plt.figure(figsize=(max(9.0, 5.6 + 0.012 * len(ordered)), max(5.2, 3.6 + 0.26 * len(selected_regulons))), dpi=150)
        grid = fig.add_gridspec(2, 1, height_ratios=(0.12, 1.0), hspace=0.08)
        ax_group = fig.add_subplot(grid[0])
        ax = fig.add_subplot(grid[1])
        group_colors = [NATURE_PALETTE[i % len(NATURE_PALETTE)] for i in range(len(groups))]
        ax_group.imshow(group_codes.reshape(1, -1), aspect="auto", interpolation="nearest",
                        cmap=ListedColormap(group_colors))
        ax_group.set_yticks([])
        ax_group.set_xticks([])
        for start, group in self._group_boundaries(shown_labels, groups):
            ax_group.axvline(start - 0.5, color="white", linewidth=0.65)
            ax.axvline(start - 0.5, color="#6B7280", linewidth=0.35, alpha=0.45)
        ax_group.set_title("Cell type / cluster", loc="left", fontsize=8, color=NATURE_TEXT, pad=2)
        ax_group.legend(
            handles=[Patch(facecolor=color, edgecolor="none", label=group)
                     for group, color in zip(groups, group_colors)],
            loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
            fontsize=7, title="Group", title_fontsize=8, ncol=1,
        )
        image = ax.imshow(shown_z, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3, interpolation="nearest")
        ax.set_yticks(np.arange(len(selected_regulons)))
        ax.set_yticklabels([regulons[index] for index in selected_regulons], fontsize=8)
        ax.set_xticks([])
        ax.set_xlabel(f"Cells (deterministic display subset: {len(ordered)}/{len(labels)})", fontsize=9)
        ax.set_title("Top variable AUCell regulons (row z-score for display)", loc="left", fontsize=10, fontweight="semibold", color=NATURE_TEXT)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.028, pad=0.02)
        colorbar.set_label("AUCell z-score", fontsize=8)
        colorbar.ax.tick_params(labelsize=7)
        fig.subplots_adjust(left=0.25, right=0.93, bottom=0.10, top=0.90)
        return fig, ordered

    @staticmethod
    def _group_boundaries(labels, groups):
        start = 0
        for group in groups:
            count = int(np.sum(labels == group))
            if count:
                yield start, group
                start += count

    def _build_rss_bubble(self, rss, regulons, groups, top_n):
        import matplotlib.pyplot as plt
        from modules.figure_style import NATURE_TEXT

        selected = set()
        for group_index in range(len(groups)):
            selected.update(_top_indices(rss[group_index], top_n).tolist())
        selected = sorted(selected, key=lambda index: (-float(rss[:, index].max()), regulons[index]))
        rows = []
        for group_index, group in enumerate(groups):
            for reg_index in selected:
                rows.append({"group": group, "regulon": regulons[reg_index], "rss": float(rss[group_index, reg_index])})
        frame = pd.DataFrame(rows)
        y_labels = list(reversed([regulons[index] for index in selected]))
        y_map = {name: index for index, name in enumerate(y_labels)}
        x_map = {group: index for index, group in enumerate(groups)}
        fig, ax = plt.subplots(figsize=(max(7.2, 1.05 * len(groups) + 4.2), max(5.0, 0.23 * len(y_labels) + 2.7)), dpi=150)
        scores = frame["rss"].to_numpy(dtype=float)
        scatter = ax.scatter(
            [x_map[item] for item in frame["group"]],
            [y_map[item] for item in frame["regulon"]],
            s=28 + 300 * scores, c=scores, cmap="YlOrRd", vmin=0,
            vmax=max(0.3, float(scores.max()) if scores.size else 0.0),
            edgecolors="#4B5563", linewidths=0.25, alpha=0.9,
        )
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=7)
        ax.set_xlabel("Cell type / cluster", fontsize=9)
        ax.set_ylabel("Regulon", fontsize=9)
        ax.set_title("Regulon specificity score (RSS)", loc="left", fontsize=10, fontweight="semibold", color=NATURE_TEXT)
        colorbar = fig.colorbar(scatter, ax=ax, fraction=0.032, pad=0.02)
        colorbar.set_label("RSS (1 − Jensen–Shannon distance)", fontsize=8)
        colorbar.ax.tick_params(labelsize=7)
        ax.grid(axis="x", color="#E5E7EB", linewidth=0.7)
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.tight_layout(pad=1.1)
        return fig, frame, selected

    def _build_rss_bar(self, rss, regulons, groups, top_n):
        import matplotlib.pyplot as plt
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT

        requested = str(self.params.get("rss_target_group", "") or "").strip()
        if requested and requested not in groups:
            raise ValueError(f"RSS 目标细胞类型 '{requested}' 不存在于 '{self.params.get('groupby', 'celltype')}'。")
        group = requested or groups[int(np.argmax(rss.max(axis=1)))]
        scores = rss[groups.index(group)]
        selected = _top_indices(scores, top_n)
        selected = selected[np.argsort(scores[selected], kind="mergesort")]
        fig, ax = plt.subplots(figsize=(8.0, max(4.3, 0.32 * len(selected) + 2.2)), dpi=150)
        ax.barh(np.arange(len(selected)), scores[selected], color=NATURE_PALETTE[0], edgecolor="#1F2937", linewidth=0.45)
        ax.set_yticks(np.arange(len(selected)))
        ax.set_yticklabels([regulons[index] for index in selected], fontsize=8)
        ax.set_xlabel("RSS", fontsize=9)
        selected_max = float(scores[selected].max()) if selected.size else 0.0
        ax.set_xlim(0, max(0.1, min(1.0, selected_max * 1.12)))
        ax.set_title(f"Top RSS regulons: {group}", loc="left", fontsize=10, fontweight="semibold", color=NATURE_TEXT)
        ax.grid(axis="x", color="#E5E7EB", linewidth=0.7)
        ax.set_axisbelow(True)
        for spine_name, spine in ax.spines.items():
            spine.set_visible(spine_name in {"left", "bottom"})
        fig.tight_layout(pad=1.1)
        return fig, group, selected

    def _build_network(self, targets, selected_regulons, regulons, rss, groups, max_targets):
        import matplotlib.pyplot as plt
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT

        if not targets:
            return None, pd.DataFrame()
        rows = []
        for reg_index in selected_regulons:
            regulon = regulons[reg_index]
            candidate_targets = targets.get(regulon, targets.get(_regulon_tf(regulon), []))
            for target in list(dict.fromkeys(candidate_targets))[:max_targets]:
                rows.append({
                    "regulon": regulon, "tf": _regulon_tf(regulon), "target": target,
                    "max_rss": float(rss[:, reg_index].max()),
                    "top_specific_group": groups[int(np.argmax(rss[:, reg_index]))],
                })
        edges = pd.DataFrame(rows)
        if edges.empty:
            return None, edges

        tfs = list(dict.fromkeys(edges["tf"].tolist()))
        target_names = list(dict.fromkeys(edges["target"].tolist()))
        fig, ax = plt.subplots(figsize=(10.0, max(5.8, 0.18 * len(target_names) + 2.5)), dpi=150)
        tf_y = np.linspace(0.92, 0.08, len(tfs))
        target_y = np.linspace(0.96, 0.04, len(target_names))
        tf_pos = {name: (0.17, tf_y[index]) for index, name in enumerate(tfs)}
        target_pos = {name: (0.80, target_y[index]) for index, name in enumerate(target_names)}
        max_rss = edges.groupby("tf")["max_rss"].max().to_dict()
        for row in edges.itertuples(index=False):
            x1, y1 = tf_pos[row.tf]
            x2, y2 = target_pos[row.target]
            ax.plot([x1, x2], [y1, y2], color="#9CA3AF", linewidth=0.55 + 1.3 * row.max_rss, alpha=0.55, zorder=1)
        for index, tf in enumerate(tfs):
            x, y = tf_pos[tf]
            ax.scatter(x, y, s=150 + 330 * max_rss[tf], color=NATURE_PALETTE[index % len(NATURE_PALETTE)], edgecolors="#111827", linewidths=0.55, zorder=2)
            ax.text(x - 0.018, y, tf, ha="right", va="center", fontsize=8, color=NATURE_TEXT)
        for target in target_names:
            x, y = target_pos[target]
            ax.scatter(x, y, s=38, color="#FFFFFF", edgecolors="#475569", linewidths=0.5, zorder=2)
            ax.text(x + 0.018, y, target, ha="left", va="center", fontsize=7, color=NATURE_TEXT)
        ax.text(0.17, 1.02, "TF / regulon", ha="center", va="bottom", fontsize=9, fontweight="semibold")
        ax.text(0.80, 1.02, "Target gene", ha="center", va="bottom", fontsize=9, fontweight="semibold")
        ax.set_title("Top RSS regulon target network", loc="left", fontsize=10, fontweight="semibold", color=NATURE_TEXT)
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.04, 1.08)
        ax.axis("off")
        fig.tight_layout(pad=1.0)
        return fig, edges

    def _build_correlation_heatmap(self, matrix, regulons, selected_regulons):
        from modules.native_figures import heatmap_figure

        selected = np.asarray(selected_regulons, dtype=int)
        frame = pd.DataFrame(matrix[:, selected], columns=[regulons[index] for index in selected])
        correlation = frame.corr(method="spearman").fillna(0.0)
        figure = heatmap_figure(
            correlation.to_numpy(), x_labels=correlation.columns.tolist(), y_labels=correlation.index.tolist(),
            title="Regulon co-activity correlation (Spearman)", x_label="Regulon", y_label="Regulon",
            colorbar_label="Spearman ρ", vmin=-1, vmax=1,
        )
        return figure, correlation

    def run(self, input_path):
        self.progress(5, "加载包含本地 AUCell 矩阵的 AnnData...")
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata, self.MODULE_NAME)
        if adata.n_obs < 3:
            raise ValueError("SCENIC regulon 结果分析至少需要 3 个细胞。")
        matrix, regulons, auc_key = self._resolve_auc_matrix(adata)
        group_key, labels, groups = self._resolve_group(adata)
        targets, target_key = self._resolve_targets(adata)

        top_n = _positive_int(self.params.get("top_n_regulons"), default=20, minimum=2,
                              label="AUCell 热图 regulon 数")
        max_heatmap_cells = _positive_int(self.params.get("max_heatmap_cells"), default=500, minimum=1,
                                           label="AUCell 热图展示细胞数")
        rss_top_n = _positive_int(self.params.get("top_n_rss"), default=10, minimum=1,
                                  label="每群 RSS Top regulon 数")
        aucell_variance = np.nanvar(matrix, axis=0)
        selected_regulons = _top_indices(aucell_variance, top_n)
        results_dir = Path(self.project_dir) / "results" / self.MODULE_NAME / _safe_name(self.params.get("_analysis_id"))
        results_dir.mkdir(parents=True, exist_ok=True)
        plots_dir = self.ensure_plots_dir()
        result_files = []

        def save_table(frame, filename, label):
            path = results_dir / filename
            frame.to_csv(path, index=False)
            result_files.append({"file_path": str(path), "file_type": "csv", "category": "table", "label": label})

        self.progress(20, "生成 AUCell 调控子活性热图...")
        heatmap_figure, heatmap_cells = self._build_aucell_heatmap(
            matrix, regulons, labels, groups, selected_regulons, max_heatmap_cells,
        )
        result_files.extend(self.save_matplotlib_figure(
            heatmap_figure, plots_dir, "scenic_aucell_heatmap.png", "heatmap",
            "SCENIC AUCell regulon activity heatmap", formats=("png", "svg"), dpi=300,
            preserve_aspect=True,
        ))
        heatmap_selection = pd.DataFrame({
            "regulon": [regulons[index] for index in selected_regulons],
            "aucell_variance": aucell_variance[selected_regulons],
            "selection_rank": np.arange(1, len(selected_regulons) + 1),
        })
        save_table(heatmap_selection, "aucell_heatmap_regulons.csv", "AUCell heatmap regulon selection")
        summary_rows = []
        for group in groups:
            group_matrix = matrix[labels == group]
            for reg_index, regulon in enumerate(regulons):
                summary_rows.append({
                    "groupby": group_key, "cell_group": group, "regulon": regulon,
                    "n_cells": int(group_matrix.shape[0]),
                    "mean_aucell": float(group_matrix[:, reg_index].mean()),
                    "median_aucell": float(np.median(group_matrix[:, reg_index])),
                })
        save_table(pd.DataFrame(summary_rows), "group_regulon_aucell_summary.csv", "Cell-group regulon AUCell summary")
        save_table(pd.DataFrame({
            "cell_id": adata.obs_names[heatmap_cells].astype(str),
            group_key: labels[heatmap_cells],
            "display_order": np.arange(1, len(heatmap_cells) + 1),
        }), "aucell_heatmap_cells.csv", "AUCell heatmap deterministic cell selection")

        self.progress(45, "计算 regulon specificity score (RSS)...")
        rss = self._rss(matrix, labels, groups)
        rss_rows = []
        for group_index, group in enumerate(groups):
            order = _top_indices(rss[group_index], len(regulons))
            ranks = np.empty(len(regulons), dtype=int)
            ranks[order] = np.arange(1, len(regulons) + 1)
            for reg_index, regulon in enumerate(regulons):
                rss_rows.append({"groupby": group_key, "cell_group": group, "regulon": regulon,
                                 "rss": float(rss[group_index, reg_index]), "rank_within_group": int(ranks[reg_index])})
        rss_frame = pd.DataFrame(rss_rows)
        save_table(rss_frame, "regulon_specificity_score.csv", "SCENIC regulon specificity score (RSS)")
        bubble_figure, bubble_frame, rss_selected = self._build_rss_bubble(rss, regulons, groups, rss_top_n)
        result_files.extend(self.save_matplotlib_figure(
            bubble_figure, plots_dir, "scenic_rss_bubble.png", "bubble",
            "SCENIC regulon specificity (RSS) bubble plot", formats=("png", "svg"), dpi=300,
            preserve_aspect=True,
        ))
        save_table(bubble_frame, "rss_bubble_regulons.csv", "RSS bubble-plot regulon selection")
        bar_figure, rss_bar_group, rss_bar_selected = self._build_rss_bar(rss, regulons, groups, rss_top_n)
        result_files.extend(self.save_matplotlib_figure(
            bar_figure, plots_dir, "scenic_rss_bar.png", "bar",
            "SCENIC top regulons by RSS", formats=("png", "svg"), dpi=300,
            preserve_aspect=True,
        ))

        network_edges = pd.DataFrame()
        if bool(self.params.get("show_network", True)):
            self.progress(65, "整理 regulon 靶基因网络...")
            max_network_targets = _positive_int(
                self.params.get("network_max_targets_per_regulon"), default=10, minimum=1,
                label="每个 regulon 展示靶基因数",
            )
            network_figure, network_edges = self._build_network(
                targets, rss_selected, regulons, rss, groups, max_network_targets,
            )
            if network_figure is not None:
                result_files.extend(self.save_matplotlib_figure(
                    network_figure, plots_dir, "scenic_regulon_network.png", "network",
                    "SCENIC regulon target network", formats=("png", "svg"), dpi=300,
                    preserve_aspect=True,
                ))
                save_table(network_edges, "regulon_network_edges.csv", "SCENIC regulon target-network edges")

        correlation = pd.DataFrame()
        if bool(self.params.get("show_correlation", True)):
            self.progress(78, "计算 regulon 共活性相关性...")
            correlation_n = _positive_int(
                self.params.get("top_n_correlation_regulons"), default=20, minimum=2,
                label="相关热图 regulon 数",
            )
            correlation_indices = _top_indices(aucell_variance, correlation_n)
            correlation_figure, correlation = self._build_correlation_heatmap(matrix, regulons, correlation_indices)
            result_files.extend(self.save_matplotlib_figure(
                correlation_figure, plots_dir, "scenic_regulon_correlation.png", "heatmap",
                "SCENIC regulon co-activity correlation heatmap", formats=("png", "svg"), dpi=300,
                preserve_aspect=True,
            ))
            correlation_out = correlation.reset_index(names="regulon")
            save_table(correlation_out, "regulon_coactivity_correlation.csv", "SCENIC regulon co-activity correlation")

        self.progress(90, "保存 SCENIC 结果溯源与 AnnData...")
        source_method = ""
        scenic_metadata = adata.uns.get("scenic")
        if isinstance(scenic_metadata, dict):
            source_method = str(scenic_metadata.get("method", scenic_metadata.get("pipeline", "")) or "")
        manifest = {
            "module": self.MODULE_NAME,
            "activity_method": "precomputed_AUCell",
            "activity_matrix": {"obsm_key": auc_key, "n_cells": int(adata.n_obs), "n_regulons": int(len(regulons))},
            "upstream_scenic_method": source_method or None,
            "grouping": {"key": group_key, "groups": groups, "n_groups": len(groups)},
            "core_outputs": {
                "aucell_heatmap": {"regulons": [regulons[index] for index in selected_regulons], "n_display_cells": int(len(heatmap_cells))},
                "rss": {"formula": "1 - Jensen-Shannon distance (base 2)", "top_n_per_group": rss_top_n, "bar_group": rss_bar_group},
            },
            "network": {
                "targets_uns_key": target_key or None,
                "status": "available" if not network_edges.empty else "not_available",
                "interpretation": "Edges are imported SCENIC regulon memberships, not inferred by this display module.",
            },
            "coactivity_correlation": {
                "status": "available" if not correlation.empty else "not_run",
                "method": "Spearman correlation across cells",
                "interpretation": "Co-activity is association, not evidence of a directed regulatory edge.",
            },
            "interpretation": [
                "AUCell/RSS represent inferred regulon activity, not TF mRNA abundance.",
                "RSS is a cell-group specificity ranking; it is not a condition-level differential test.",
                "For condition conclusions, validate hypotheses with sample-level pseudobulk DEG/GSEA and independent biological replicates.",
            ],
        }
        manifest_path = results_dir / "analysis_manifest.json"
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
        result_files.append({"file_path": str(manifest_path), "file_type": "json", "category": "provenance", "label": "SCENIC analysis manifest"})
        adata.uns["scenic_result_display"] = manifest
        output_path = self.save_output(adata, self.MODULE_NAME)
        warnings = []
        if not targets:
            warnings.append("未找到 regulon 靶基因定义：已生成 AUCell/RSS/共活性结果，但跳过网络图。")
        if len(heatmap_cells) < adata.n_obs:
            warnings.append("AUCell 热图为保证可读性使用确定性分层展示子集；完整矩阵保留在输出 H5AD 的 obsm 中。")
        self.progress(100, "SCENIC regulon 活性结果分析完成。")
        return {
            "output_adata": output_path,
            "result_files": result_files,
            "summary": {
                "n_cells": int(adata.n_obs), "n_regulons": int(len(regulons)), "n_groups": int(len(groups)),
                "aucell_obsm_key": auc_key, "groupby": group_key, "rss_bar_group": rss_bar_group,
                "n_network_edges": int(len(network_edges)), "n_correlation_regulons": int(len(correlation)),
                "warnings": warnings,
            },
        }
