import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis
from modules.io_utils import obs_grouping_info, restore_scanpy_qc_percentages


class QCReassessAnalysis(BaseAnalysis):
    MODULE_NAME = "qc_reassess"
    DISPLAY_NAME = "QC 重新评估"
    DESCRIPTION = "聚类后检查 doublet 和 QC 指标，标记低质量簇"
    INPUT_REQUIRES = ['leiden']

    @staticmethod
    def _as_bool(value, default=False):
        """把表单/API 传来的布尔值统一解析，避免字符串 'false' 被当真。"""
        if value is None:
            return default
        if isinstance(value, str):
            return value.strip().lower() not in {'', '0', 'false', 'no', 'off', 'none'}
        return bool(value)

    def run(self, input_path):
        import scanpy as sc
        import matplotlib.pyplot as plt
        from modules.native_figures import umap_figure, umap_panel_figure
        from modules.figure_style import NATURE_PALETTE, NATURE_GRID, NATURE_TEXT

        self.progress(5, "加载数据...")
        adata = self.load_adata(input_path)
        # 兼容旧版 QC 输出：旧文件的 pct_counts_* 可能被写成
        # 0-1 fraction，而本模块的阈值和图均使用 0-100 percentage。
        restore_scanpy_qc_percentages(adata)

        cluster_key = self.params.get('cluster_key', 'leiden')
        doublet_threshold = float(self.params.get('doublet_threshold', 0.3))
        doublet_score_cutoff = float(self.params.get('doublet_score_cutoff', 0.5))
        mt_threshold = float(self.params.get('mt_threshold', 15.0))
        ribo_threshold = float(self.params.get('ribosomal_threshold', 0))
        min_cells = int(self.params.get('min_cells_per_cluster', 10))
        auto_remove = self._as_bool(self.params.get('auto_remove', False), False)
        if not 0 <= doublet_threshold <= 1:
            raise ValueError('doublet_threshold 必须位于 [0, 1]。')
        if not 0 <= doublet_score_cutoff <= 1:
            raise ValueError('doublet_score_cutoff 必须位于 [0, 1]。')

        # A QC metric (for example ``log1p_n_genes_by_counts``) is numeric and
        # can have nearly one value per cell.  Treating it as a cluster column
        # creates hundreds of pseudo-clusters, marks almost every cell as
        # low-quality, and makes the QC charts unusable.  Prefer a real Leiden
        # grouping when the requested column is missing or over-cardinal.
        requested_cluster_key = str(cluster_key)
        def _valid_cluster_key(key):
            return bool(obs_grouping_info(
                adata, key, max_categories=50,
                max_numeric_categories=20, require_multiple=False,
            )['valid'])

        if not _valid_cluster_key(cluster_key):
            fallback_keys = [
                'leiden', 'leiden_0.8', 'leiden_0.6', 'leiden_1.0',
            ] + [key for key in adata.obs.columns if str(key).startswith('leiden_')]
            cluster_key = next((key for key in fallback_keys if _valid_cluster_key(key)), None)
            if cluster_key is None:
                raise ValueError(
                    f"cluster_key '{requested_cluster_key}' 不是有效的分类聚类列；"
                    "请改用 leiden 或其他低基数 cluster 列。"
                )
            self.progress(
                12,
                f"'{requested_cluster_key}' 不是有效 cluster 列，已改用 '{cluster_key}'",
            )

        self.progress(20, "计算各簇 QC 指标...")
        result_files = []
        plots_dir = self.ensure_plots_dir()
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)

        # 计算每个簇的统计
        clusters = adata.obs[cluster_key].unique()
        if 'predicted_doublet' in adata.obs.columns:
            doublet_fraction_source = 'predicted_doublet'
        elif 'doublet_score' in adata.obs.columns:
            doublet_fraction_source = f'doublet_score >= {doublet_score_cutoff:g}'
        else:
            doublet_fraction_source = 'unavailable'
        cluster_stats = []
        for c in sorted(clusters):
            mask = adata.obs[cluster_key] == c
            n_cells = mask.sum()
            doublet_frac = np.nan
            mean_doublet_score = np.nan
            if 'predicted_doublet' in adata.obs.columns:
                doublet_frac = float(
                    adata.obs.loc[mask, 'predicted_doublet'].astype(bool).mean()
                )
            elif 'doublet_score' in adata.obs.columns:
                doublet_frac = float(
                    (pd.to_numeric(adata.obs.loc[mask, 'doublet_score'], errors='coerce')
                     >= doublet_score_cutoff).mean()
                )
            if 'doublet_score' in adata.obs.columns:
                mean_doublet_score = float(pd.to_numeric(
                    adata.obs.loc[mask, 'doublet_score'], errors='coerce'
                ).mean())
            mt_mean = adata.obs.loc[mask, 'pct_counts_mt'].mean() if 'pct_counts_mt' in adata.obs.columns else 0
            counts_mean = adata.obs.loc[mask, 'total_counts'].mean() if 'total_counts' in adata.obs.columns else 0
            genes_mean = adata.obs.loc[mask, 'n_genes_by_counts'].mean() if 'n_genes_by_counts' in adata.obs.columns else 0
            is_low = False
            reasons = []
            if pd.notna(doublet_frac) and doublet_frac > doublet_threshold:
                is_low = True
                reasons.append('high_doublet')
            if mt_mean > mt_threshold:
                is_low = True
                reasons.append('high_mt')
            if ribo_threshold > 0:
                ribo_mean = adata.obs.loc[mask, 'pct_counts_ribo'].mean() if 'pct_counts_ribo' in adata.obs.columns else 0
                if ribo_mean > ribo_threshold:
                    is_low = True
                    reasons.append('high_ribo')
            if n_cells < min_cells:
                is_low = True
                reasons.append('too_few_cells')
            cluster_stats.append({
                'cluster': str(c), 'n_cells': n_cells,
                'doublet_fraction': round(doublet_frac, 4),
                'mean_doublet_score': round(mean_doublet_score, 4),
                'doublet_fraction_source': doublet_fraction_source,
                'mean_pct_mt': round(mt_mean, 2),
                'mean_total_counts': round(counts_mean, 0),
                'mean_n_genes': round(genes_mean, 0),
                'low_quality': is_low,
                'low_reasons': ', '.join(reasons),
            })

        stats_df = pd.DataFrame(cluster_stats)
        csv_path = os.path.join(results_dir, 'low_quality_clusters.csv')
        stats_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '簇 QC 统计'})

        n_low = stats_df['low_quality'].sum()

        if self.params.get('show_cluster_qc_bar', True) and not stats_df.empty:
            color_by_quality = [NATURE_PALETTE[3] if value else NATURE_PALETTE[0]
                                for value in stats_df['low_quality']]
            fig_qc_bar, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), dpi=150)
            metrics = [('n_cells', 'Cell count'), ('mean_n_genes', 'Mean detected genes'),
                       ('mean_pct_mt', 'Mean MT%'), ('doublet_fraction', 'Doublet fraction')]
            for ax, (column, label) in zip(axes.ravel(), metrics):
                ax.bar(stats_df['cluster'].astype(str), stats_df[column], color=color_by_quality,
                       alpha=0.88, edgecolor='white', linewidth=0.3)
                ax.set_title(label, loc='left', fontsize=9, color=NATURE_TEXT)
                ax.set_xlabel('Cluster', fontsize=8)
                ax.tick_params(axis='x', rotation=35, labelsize=7)
                ax.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
            fig_qc_bar.suptitle('Cluster QC Summary', x=0.05, ha='left', fontsize=11,
                                fontweight='semibold', color=NATURE_TEXT)
            result_files.extend(self.save_matplotlib_figure(
                fig_qc_bar, plots_dir, 'qc_reassess_cluster_qc_bar.png',
                'bar', 'Cluster QC Summary', formats=('png', 'svg'), dpi=300,
                preserve_aspect=True,
            ))

        # Auto-remove low quality clusters
        if auto_remove and n_low > 0:
            low_clusters = set(stats_df[stats_df['low_quality']]['cluster'].tolist())
            keep_mask = ~adata.obs[cluster_key].astype(str).isin(low_clusters)
            if int(keep_mask.sum()) == 0:
                raise ValueError(
                    'auto_remove 会移除全部细胞；请放宽 QC 阈值或关闭 auto_remove。'
                )
            n_before = adata.n_obs
            adata = adata[keep_mask].copy()
            self.progress(45, f"Removed {n_before - adata.n_obs} cells from {len(low_clusters)} low-quality clusters")

        self.progress(50, "生成 UMAP 图...")
        # UMAP with doublet score
        if 'X_umap' in adata.obsm and 'doublet_score' in adata.obs.columns:
            fig = umap_figure(adata, 'doublet_score', title='Doublet Score')
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'qc_reassess_doublet_umap.png', 'umap',
                'Doublet Score UMAP', formats=('png', 'svg'), dpi=300,
            ))

        # UMAP with MT percentage
        if 'X_umap' in adata.obsm and 'pct_counts_mt' in adata.obs.columns:
            fig = umap_figure(adata, 'pct_counts_mt', title='MT Percentage')
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'qc_reassess_mt_umap.png', 'umap',
                'MT% UMAP', formats=('png', 'svg'), dpi=300,
            ))

        # QC 指标 UMAP 面板
        if self.params.get('show_qc_umap_panel', True) and 'X_umap' in adata.obsm:
            qc_metrics = [
                ('pct_counts_mt', 'MT%'),
                ('n_genes_by_counts', 'Detected genes'),
                ('total_counts', 'Total counts'),
                ('doublet_score', 'Doublet score'),
            ]
            qc_metrics = [(col, label) for col, label in qc_metrics if col in adata.obs.columns]
            if qc_metrics:
                fig_panel = umap_panel_figure(
                    adata, [col for col, _ in qc_metrics[:4]],
                    titles=[label for _, label in qc_metrics[:4]], point_size=5,
                    opacity=0.78, ncols=2,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_panel, plots_dir, 'qc_reassess_metrics_umap_panel.png',
                    'umap', 'QC 指标 UMAP 面板', formats=('png', 'svg'), dpi=300,
                    preserve_aspect=True,
                ))

        # UMAP with cluster highlighting (low-quality clusters in red)
        if 'X_umap' in adata.obsm:
            low_quality_set = set(stats_df[stats_df['low_quality']]['cluster'].tolist())
            cluster_labels = adata.obs[cluster_key].astype(str)
            fig = umap_figure(
                adata, cluster_key,
                title=f'Clusters ({cluster_key}) — Low Quality Highlighted',
                label_categories=True,
            )
            ax = fig.axes[0]
            if low_quality_set:
                low_mask = cluster_labels.isin(low_quality_set).values
                coords = adata.obsm['X_umap'][:, :2]
                ax.scatter(coords[low_mask, 0], coords[low_mask, 1], s=22,
                           facecolors='none', edgecolors=NATURE_PALETTE[3],
                           linewidths=0.9, label='Low Quality Highlight')
                ax.legend(frameon=False, fontsize=8)
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'qc_reassess_clusters_umap.png', 'umap',
                '聚类 UMAP（低质量簇高亮）', formats=('png', 'svg'), dpi=300,
            ))

        self.progress(90, "保存输出...")
        output_path = self.save_output(adata, 'qc_reassess')

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_clusters': len(clusters),
                'cluster_key': cluster_key,
                'requested_cluster_key': requested_cluster_key,
                'n_low_quality': int(n_low),
                'low_quality_clusters': stats_df[stats_df['low_quality']]['cluster'].tolist(),
                'doublet_threshold': doublet_threshold,
                'doublet_score_cutoff': doublet_score_cutoff,
                'doublet_fraction_source': doublet_fraction_source,
                'mt_threshold': mt_threshold,
            }
        }
