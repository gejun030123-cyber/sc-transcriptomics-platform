import os
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkQCAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_qc"
    DISPLAY_NAME = "Bulk RNA-seq 质控"
    DESCRIPTION = "计数矩阵质控：文库大小、基因检测、离群值过滤"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import matplotlib.pyplot as plt
        from modules.native_figures import (
            heatmap_figure, scatter_figure, bar_figure,
            correlation_heatmap_figure, correlation_pairwise_table,
            summarize_correlation_pairs,
        )
        from modules.figure_style import NATURE_PALETTE, NATURE_GRID, NATURE_TEXT

        self.progress(5, "加载计数矩阵...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)
        from modules.io_utils import infer_expression_measurement
        input_measurement = infer_expression_measurement(adata, input_path)

        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'bulk_qc')

        min_counts = int(self.params.get('min_counts', 100000))
        min_genes = int(self.params.get('min_genes', 5000))
        max_mt_pct = float(self.params.get('max_mt_pct', 20.0))
        max_ribo_pct = float(self.params.get('max_ribo_pct', 40.0))
        min_gini = float(self.params.get('min_gini', 0))
        min_sample_expr = int(self.params.get('min_sample_expr', 0))
        min_count_threshold = int(self.params.get('min_count_threshold', 1))
        detect_outliers = self.params.get('detect_outliers', True)
        filter_strategy = self.params.get('filter_strategy', 'standard')

        # 过滤策略预设覆盖阈值
        if filter_strategy == 'strict':
            min_counts = max(min_counts, 200000)
            min_genes = max(min_genes, 8000)
            max_mt_pct = min(max_mt_pct, 15.0)
            max_ribo_pct = min(max_ribo_pct, 30.0)

        if input_measurement != 'raw_counts':
            self.progress(-1, '检测到连续或已标准化表达值：跳过 count 文库大小、MT/Ribo 百分比硬过滤。')
            min_counts, min_genes, max_mt_pct, max_ribo_pct = 0, 0, 100.0, 100.0

        self.progress(20, "计算质控指标...")
        # 优先用 gene_name 检测线粒体基因（Ensembl ID 不以 MT- 开头）
        if 'gene_name' in adata.var.columns:
            # h5ad 会将重复字符串列读为 categorical；先转 StringDtype 再填空值。
            gene_names_for_mt = adata.var['gene_name'].astype('string').fillna('')
        else:
            gene_names_for_mt = adata.var_names.astype(str)
        adata.var['mt'] = gene_names_for_mt.str.startswith('MT-')
        adata.var['ribo'] = gene_names_for_mt.str.startswith(('RPL', 'RPS'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

        n_before = adata.n_obs
        lib_sizes = adata.obs['total_counts'].values
        n_genes_detected = adata.obs['n_genes_by_counts'].values
        mt_pct = adata.obs['pct_counts_mt'].values if 'pct_counts_mt' in adata.obs.columns else np.zeros(n_before)
        ribo_pct = adata.obs['pct_counts_ribo'].values if 'pct_counts_ribo' in adata.obs.columns else np.zeros(n_before)

        # 文库复杂度 (Gini) 和新颖度
        raw_counts = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
        gini_values = np.array([_gini(raw_counts[i]) for i in range(n_before)])
        novelty_values = np.log10(n_genes_detected + 1) / np.log10(lib_sizes + 1)

        self.progress(40, "过滤样本...")
        # 分组推断
        sample_names = adata.obs.index.tolist()
        group_col = self.params.get('group_column', '').strip()
        if group_col == '_auto_group_':
            auto_mapping = self.params.get('_auto_group_mapping', {})
            if isinstance(auto_mapping, str):
                try:
                    import json as _json
                    auto_mapping = _json.loads(auto_mapping)
                except (TypeError, ValueError):
                    auto_mapping = {}
            if not auto_mapping:
                from modules.io_utils import infer_sample_group_candidates
                candidates = infer_sample_group_candidates(sample_names)
                auto_mapping = candidates[0]['mapping'] if candidates else {}
            groups = [str(auto_mapping.get(str(name), 'unknown')) for name in sample_names]
            adata.obs['_auto_group'] = groups
            group_col = '_auto_group'
        elif group_col and group_col in adata.obs.columns:
            groups = adata.obs[group_col].astype(str).tolist()
        else:
            groups = _infer_groups(sample_names)
            adata.obs['_auto_group'] = groups
            group_col = '_auto_group'

        # 样本过滤
        mask = (lib_sizes >= min_counts) & (n_genes_detected >= min_genes) & (mt_pct <= max_mt_pct) & (ribo_pct <= max_ribo_pct)
        if min_gini > 0:
            mask = mask & (gini_values >= min_gini)

        # 过滤日志
        filter_log_rows = []
        for i in range(n_before):
            fail_reasons = []
            if lib_sizes[i] < min_counts:
                fail_reasons.append(f'lib_size<{min_counts}')
            if n_genes_detected[i] < min_genes:
                fail_reasons.append(f'n_genes<{min_genes}')
            if mt_pct[i] > max_mt_pct:
                fail_reasons.append(f'mt_pct>{max_mt_pct}')
            if ribo_pct[i] > max_ribo_pct:
                fail_reasons.append(f'ribo_pct>{max_ribo_pct}')
            if min_gini > 0 and gini_values[i] < min_gini:
                fail_reasons.append(f'gini<{min_gini}')
            filter_log_rows.append({
                'sample': sample_names[i],
                'group': groups[i],
                'passed': len(fail_reasons) == 0,
                'lib_size': int(lib_sizes[i]),
                'n_genes': int(n_genes_detected[i]),
                'mt_pct': round(float(mt_pct[i]), 2),
                'ribo_pct': round(float(ribo_pct[i]), 2),
                'gini': round(float(gini_values[i]), 4),
                'novelty': round(float(novelty_values[i]), 4),
                'fail_reasons': '; '.join(fail_reasons) if fail_reasons else '',
            })

        adata_filtered = adata[mask].copy()
        gini_filtered = gini_values[mask]
        n_after = adata_filtered.n_obs

        # 输出样本指标表和过滤日志
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        metrics_df = pd.DataFrame(filter_log_rows)
        metrics_csv = os.path.join(results_dir, 'bulk_qc_sample_metrics.csv')
        metrics_df.to_csv(metrics_csv, index=False)

        filter_log_df = metrics_df[~metrics_df['passed']]
        filter_log_csv = os.path.join(results_dir, 'bulk_qc_filter_log.csv')
        filter_log_df.to_csv(filter_log_csv, index=False)

        result_files.append({'file_path': metrics_csv, 'file_type': 'csv', 'category': 'table', 'label': '样本 QC 指标'})
        result_files.append({'file_path': filter_log_csv, 'file_type': 'csv', 'category': 'table', 'label': '过滤日志'})

        # 基因层面过滤（基于 raw count 阈值）
        genes_before_filter = adata_filtered.n_vars
        gene_filter_rows = []
        if min_sample_expr > 0:
            raw_filt = adata_filtered.X.toarray() if hasattr(adata_filtered.X, 'toarray') else np.asarray(adata_filtered.X)
            expr_count_per_gene = (raw_filt >= min_count_threshold).sum(axis=0)
            gene_mask = expr_count_per_gene >= min_sample_expr

            removed_indices = np.where(~gene_mask)[0]
            gene_filter_rows = [
                {
                    'gene': adata_filtered.var_names[idx],
                    'expressed_in_n_samples': int(expr_count_per_gene[idx]),
                    'action': 'removed',
                }
                for idx in removed_indices
            ]
            adata_filtered = adata_filtered[:, gene_mask].copy()

        # 管家基因稳定性检查
        housekeeping_genes = ['GAPDH', 'ACTB', 'B2M', 'HPRT1', 'TBP', 'UBC', 'YWHAZ', 'SDHA', 'HMBS', 'RPLP0']
        found_hk = [g for g in housekeeping_genes if g in adata_filtered.var_names]

        # 基因过滤日志
        if min_sample_expr > 0 and gene_filter_rows:
            gene_filter_csv = os.path.join(results_dir, 'bulk_qc_gene_filter.csv')
            pd.DataFrame(gene_filter_rows).to_csv(gene_filter_csv, index=False)
            result_files.append({'file_path': gene_filter_csv, 'file_type': 'csv', 'category': 'table', 'label': '基因过滤日志'})

        self.progress(60, "生成质控图表...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)

        fig, axes = plt.subplots(3, 2, figsize=(10.0, 11.0), dpi=150)
        sample_idx = list(range(n_before))
        colors = ['#4caf50' if m else '#e53935' for m in mask]
        overview = [
            (lib_sizes, '文库大小分布', 'Library size'),
            (n_genes_detected, '检测基因数', 'Detected genes'),
            (mt_pct, '线粒体基因比例', 'MT%'),
            (ribo_pct, '核糖体基因比例', 'Ribo%'),
            (gini_values, 'Gini 系数', 'Gini'),
        ]
        for ax, (values, subtitle, ylabel), color in zip(
                axes.ravel()[:5], overview,
                [NATURE_PALETTE[0], NATURE_PALETTE[1], NATURE_PALETTE[3],
                 NATURE_PALETTE[4], NATURE_PALETTE[5]]):
            ax.bar(sample_idx, values, color=color, alpha=0.88,
                   edgecolor='white', linewidth=0.25)
            ax.set_title(subtitle, loc='left', fontsize=9, color=NATURE_TEXT)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
        ax = axes.ravel()[5]
        ax.scatter(lib_sizes, n_genes_detected, c=colors, s=28, alpha=0.82,
                   linewidths=0)
        ax.set_title('文库大小 vs 检测基因数', loc='left', fontsize=9, color=NATURE_TEXT)
        ax.set_xlabel('Library size', fontsize=8)
        ax.set_ylabel('Detected genes', fontsize=8)
        ax.grid(False)
        fig.suptitle('Bulk RNA-seq 质控总览', x=0.05, ha='left', fontsize=12,
                     fontweight='semibold', color=NATURE_TEXT)
        result_files.extend(self.save_matplotlib_figure(
            fig, plots_dir, 'bulk_qc_overview.png', 'qc', '质控总览',
            formats=('png', 'svg'), dpi=300,
        ))

        # 保存原始 counts 副本
        adata_raw_filtered = adata_filtered.copy()

        # 统一标准化一份副本，用于相关性热图和 PCA
        self.progress(65, "标准化数据...")
        adata_normed = adata_filtered.copy()
        if input_measurement == 'raw_counts':
            sc.pp.normalize_total(adata_normed, target_sum=1e6)
            sc.pp.log1p(adata_normed)
        else:
            adata_normed.X = np.log2(np.maximum(adata_normed.X, 0) + 1)

        self.progress(70, "生成相关性热图...")
        corr_data = adata_normed.X if not hasattr(adata_normed.X, 'toarray') else adata_normed.X.toarray()
        corr_matrix = np.corrcoef(corr_data)
        sample_labels_corr = adata_normed.obs.index.tolist()

        # 分组条与样本排序均基于过滤后的样本；相关性原始矩阵不改动。
        filtered_groups = [groups[sample_names.index(s)] for s in sample_labels_corr]
        unique_groups = sorted(set(filtered_groups))
        fig_corr, corr_metadata = correlation_heatmap_figure(
            corr_matrix, sample_labels_corr, title='样本相关性热图 (Pearson)',
            method='Pearson', group_labels=filtered_groups,
            colorscale='Blues', cluster=True, mask_diagonal=True,
        )
        result_files.extend(self.save_matplotlib_figure(
            fig_corr, plots_dir, 'bulk_qc_corr.png', 'heatmap',
            '样本相关性热图', formats=('png', 'svg'), dpi=300,
        ))
        corr_pairs = correlation_pairwise_table(
            corr_matrix, sample_labels_corr, filtered_groups, method='pearson')
        corr_pairs_csv = os.path.join(results_dir, 'bulk_qc_correlation_pairs.csv')
        corr_pairs.to_csv(corr_pairs_csv, index=False)
        corr_summary = summarize_correlation_pairs(corr_pairs, method='pearson')
        corr_summary_csv = os.path.join(results_dir, 'bulk_qc_correlation_summary.csv')
        corr_summary.to_csv(corr_summary_csv, index=False)
        result_files.extend([
            {'file_path': corr_pairs_csv, 'file_type': 'csv', 'category': 'table',
             'label': '样本两两 Pearson 相关性'},
            {'file_path': corr_summary_csv, 'file_type': 'csv', 'category': 'table',
             'label': '样本相关性组内/组间摘要'},
        ])

        if n_before > n_after:
            fig_r = bar_figure(
                ['过滤前', '过滤后'], [n_before, n_after], title='样本过滤结果',
                x_label='Stage', y_label='样本数',
                colors=[NATURE_PALETTE[3], NATURE_PALETTE[0]], rotation=0,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_r, plots_dir, 'bulk_qc_filter.png', 'qc', '样本过滤结果',
                formats=('png', 'svg'), dpi=300,
            ))

        self.progress(75, "PCA 离群检测...")
        outlier_samples = []
        n_comps = min(10, n_after - 1, adata_normed.n_vars - 1)
        if n_comps >= 2:
            sc.pp.pca(adata_normed, n_comps=n_comps)
            pc = adata_normed.obsm['X_pca']
            filtered_sample_names = adata_normed.obs.index.tolist()
            filtered_groups_pca = [groups[sample_names.index(s)] for s in filtered_sample_names]
            unique_groups_pca = sorted(set(filtered_groups_pca))
            group_color_map_pca = {g: NATURE_PALETTE[i % len(NATURE_PALETTE)]
                                   for i, g in enumerate(unique_groups_pca)}

            # 离群检测
            if detect_outliers:
                outlier_samples = _detect_outliers_mahal(pc, filtered_sample_names)

            # Reuse the display PCA canvas: sample labels are intentionally
            # omitted here because naming every replicate masks compact groups.
            from modules.bulk_pca import _pca_embedding_figure
            pca_variance_qc = adata_normed.uns.get('pca', {}).get('variance_ratio', [])
            pc1_var = float(pca_variance_qc[0]) * 100 if len(pca_variance_qc) > 0 else None
            pc2_var = float(pca_variance_qc[1]) * 100 if len(pca_variance_qc) > 1 else None
            axis_suffix = (
                f' ({pc1_var:.1f}% variance)' if pc1_var is not None else ''
            )
            fig_pca = _pca_embedding_figure(
                pc[:, :2], filtered_groups_pca, filtered_sample_names,
                title='质控后样本 PCA',
                x_label=f'PC1{axis_suffix}',
                y_label=(f'PC2 ({pc2_var:.1f}% variance)' if pc2_var is not None else 'PC2'),
                group_label='Sample group', show_labels=False,
                subtitle='样本重复以点表示；名称已隐藏以避免遮挡',
            )
            ax_pca = fig_pca.axes[0]
            # 离群点高亮
            if outlier_samples:
                out_idx = [filtered_sample_names.index(s) for s in outlier_samples if s in filtered_sample_names]
                if out_idx:
                    ax_pca.scatter(pc[out_idx, 0], pc[out_idx, 1], s=80, color='#B64342',
                                   marker='x', linewidths=1.6, label='离群样本')
                    # Rebuild the externally anchored legend to include the
                    # outlier marker without moving it back into the data area.
                    old_legend = ax_pca.get_legend()
                    if old_legend is not None:
                        old_legend.remove()
                    handles, legend_labels = ax_pca.get_legend_handles_labels()
                    legend = ax_pca.legend(
                        handles, legend_labels, title='Sample group',
                        loc='center left', bbox_to_anchor=(1.01, 0.5),
                        frameon=False, fontsize=7.5, title_fontsize=8,
                        handletextpad=0.5, borderaxespad=0,
                    )
                    for handle in legend.legend_handles:
                        if hasattr(handle, 'set_sizes'):
                            handle.set_sizes([34])
            result_files.extend(self.save_matplotlib_figure(
                fig_pca, plots_dir, 'bulk_qc_pca.png', 'pca', '样本 PCA',
                formats=('png', 'svg'), dpi=300,
            ))

        # PCA 方差解释 elbow 图（使用 scanpy 存储的 variance_ratio）
        pca_variance = adata_normed.uns.get('pca', {}).get('variance_ratio', None)
        if pca_variance is None and n_comps >= 2:
            pca_var = np.var(pc, axis=0)
            pca_variance = pca_var / pca_var.sum()
        if pca_variance is not None and len(pca_variance) > 0:
            n_pcs = len(pca_variance)
            pc_labels = [f'PC{i+1}' for i in range(n_pcs)]
            cumulative = np.cumsum(pca_variance).tolist()
            fig_elbow, ax_elbow = plt.subplots(figsize=(7.5, 5.0), dpi=150)
            x_pc = np.arange(n_pcs)
            ax_elbow.bar(x_pc, pca_variance, color=NATURE_PALETTE[0], alpha=0.88,
                         label='方差比例')
            ax2 = ax_elbow.twinx()
            ax2.plot(x_pc, cumulative, color=NATURE_PALETTE[3], linewidth=1.8,
                     marker='o', markersize=3.5, label='累积比例')
            ax_elbow.set_xticks(x_pc, pc_labels, rotation=45)
            ax_elbow.set_xlabel('主成分', fontsize=9)
            ax_elbow.set_ylabel('方差比例', fontsize=9)
            ax2.set_ylabel('累积比例', fontsize=9)
            ax_elbow.set_title('PCA 方差解释比例', loc='left', fontsize=10,
                               fontweight='semibold', color=NATURE_TEXT)
            ax_elbow.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
            handles, labels = ax_elbow.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax_elbow.legend(handles + h2, labels + l2, frameon=False, fontsize=8)
            result_files.extend(self.save_matplotlib_figure(
                fig_elbow, plots_dir, 'bulk_qc_pca_elbow.png', 'pca',
                'PCA 方差解释', formats=('png', 'svg'), dpi=300,
            ))

        # 组内 vs 组间距离箱线图
        if len(set(groups)) > 1 and n_after > 3:
            norm_arr = adata_normed.X if not hasattr(adata_normed.X, 'toarray') else adata_normed.X.toarray()
            corr_mat_all = np.corrcoef(norm_arr)
            dist_mat = 1 - corr_mat_all
            intra_dists, inter_dists = [], []
            filtered_sample_list = adata_normed.obs.index.tolist()
            for i in range(n_after):
                for j in range(i + 1, n_after):
                    gi = groups[sample_names.index(filtered_sample_list[i])]
                    gj = groups[sample_names.index(filtered_sample_list[j])]
                    if gi == gj:
                        intra_dists.append(float(dist_mat[i, j]))
                    else:
                        inter_dists.append(float(dist_mat[i, j]))

            try:
                from scipy.stats import ttest_ind
                _, pval = ttest_ind(intra_dists, inter_dists, equal_var=False)
                title_suffix = f' (p={pval:.2e})'
            except Exception:
                title_suffix = ''
            fig_dist, ax_dist = plt.subplots(figsize=(6.5, 4.8), dpi=150)
            box_data = [values for values in (intra_dists, inter_dists) if values]
            box_labels = [label for label, values in zip(['组内距离', '组间距离'],
                                                          (intra_dists, inter_dists)) if values]
            boxes = ax_dist.boxplot(box_data, tick_labels=box_labels, patch_artist=True,
                                    showfliers=False)
            for patch, color in zip(boxes['boxes'], [NATURE_PALETTE[0], NATURE_PALETTE[3]]):
                patch.set_facecolor(color)
                patch.set_alpha(0.62)
                patch.set_edgecolor(color)
            ax_dist.set_title(f'组内 vs 组间距离{title_suffix}', loc='left', fontsize=10,
                              fontweight='semibold', color=NATURE_TEXT)
            ax_dist.set_ylabel('1 - Pearson r', fontsize=9)
            ax_dist.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
            result_files.extend(self.save_matplotlib_figure(
                fig_dist, plots_dir, 'bulk_qc_group_distance.png', 'qc',
                '组内/组间距离', formats=('png', 'svg'), dpi=300,
            ))

        # QC 指标散点矩阵 (Pairs Plot)
        obs_filtered = adata_filtered.obs
        pairs_groups = [groups[sample_names.index(s)] for s in obs_filtered.index.tolist()]
        unique_pg = sorted(set(pairs_groups))
        pg_color_map = {g: NATURE_PALETTE[i % len(NATURE_PALETTE)]
                        for i, g in enumerate(unique_pg)}
        pg_colors = [pg_color_map[g] for g in pairs_groups]
        pair_values = [
            obs_filtered['total_counts'].values,
            obs_filtered['n_genes_by_counts'].values,
            obs_filtered['pct_counts_mt'].values if 'pct_counts_mt' in obs_filtered.columns else np.zeros(n_after),
            obs_filtered['pct_counts_ribo'].values if 'pct_counts_ribo' in obs_filtered.columns else np.zeros(n_after),
        ]
        pair_labels = ['Library Size', 'N Genes', 'MT%', 'Ribo%']
        fig_pairs, pair_axes = plt.subplots(4, 4, figsize=(9.0, 9.0), dpi=150)
        for row in range(4):
            for col in range(4):
                ax_pair = pair_axes[row, col]
                if row == col:
                    ax_pair.hist(pair_values[row], bins=25, color=NATURE_PALETTE[0],
                                 alpha=0.78, edgecolor='white', linewidth=0.2)
                elif row > col:
                    ax_pair.scatter(pair_values[col], pair_values[row], c=pg_colors,
                                    s=12, alpha=0.72, linewidths=0)
                else:
                    ax_pair.set_visible(False)
                if row == 3:
                    ax_pair.set_xlabel(pair_labels[col], fontsize=7)
                if col == 0 and row > col:
                    ax_pair.set_ylabel(pair_labels[row], fontsize=7)
                ax_pair.tick_params(labelsize=6)
        fig_pairs.suptitle('QC 指标散点矩阵', x=0.05, ha='left', fontsize=11,
                           fontweight='semibold', color=NATURE_TEXT)
        result_files.extend(self.save_matplotlib_figure(
            fig_pairs, plots_dir, 'bulk_qc_pairs_plot.png', 'qc',
            'QC 指标散点矩阵', formats=('png', 'svg'), dpi=300,
        ))

        # 各组 QC 指标小提琴图
        if len(unique_groups) > 1:
            violin_data = []
            for s in obs_filtered.index.tolist():
                g = groups[sample_names.index(s)]
                violin_data.append({
                    'group': g,
                    'MT%': float(obs_filtered.loc[s, 'pct_counts_mt']) if 'pct_counts_mt' in obs_filtered.columns else 0,
                    'Ribo%': float(obs_filtered.loc[s, 'pct_counts_ribo']) if 'pct_counts_ribo' in obs_filtered.columns else 0,
                    'Library Size': float(obs_filtered.loc[s, 'total_counts']),
                    'N Genes': float(obs_filtered.loc[s, 'n_genes_by_counts']),
                })
            violin_df = pd.DataFrame(violin_data)
            fig_violin, violin_axes = plt.subplots(2, 2, figsize=(9.0, 7.5), dpi=150)
            metrics = [('MT%', 1, 1), ('Ribo%', 1, 2), ('Library Size', 2, 1), ('N Genes', 2, 2)]
            for metric, row, col in metrics:
                ax_v = violin_axes[row - 1, col - 1]
                vals_by_group = [violin_df[violin_df['group'] == g][metric].values
                                  for g in unique_groups]
                parts = ax_v.violinplot(vals_by_group, positions=np.arange(len(unique_groups)),
                                        showmeans=True, showextrema=False)
                for body, index in zip(parts['bodies'], range(len(unique_groups))):
                    body.set_facecolor(NATURE_PALETTE[index % len(NATURE_PALETTE)])
                    body.set_edgecolor(NATURE_PALETTE[index % len(NATURE_PALETTE)])
                    body.set_alpha(0.62)
                ax_v.set_title(metric, loc='left', fontsize=9, color=NATURE_TEXT)
                ax_v.set_xticks(np.arange(len(unique_groups)), unique_groups, rotation=35,
                                ha='right', fontsize=7)
                ax_v.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
            fig_violin.suptitle('各组 QC 指标分布', x=0.05, ha='left', fontsize=11,
                                fontweight='semibold', color=NATURE_TEXT)
            result_files.extend(self.save_matplotlib_figure(
                fig_violin, plots_dir, 'bulk_qc_violin_by_group.png', 'qc',
                '各组 QC 指标分布', formats=('png', 'svg'), dpi=300,
            ))

        # 管家基因稳定性热图
        if found_hk:
            norm_hk = adata_raw_filtered[:, found_hk].copy()
            sc.pp.normalize_total(norm_hk, target_sum=1e6)
            sc.pp.log1p(norm_hk)
            hk_data = norm_hk.X if not hasattr(norm_hk.X, 'toarray') else norm_hk.X.toarray()
            hk_genes = norm_hk.var_names.tolist()
            hk_samples = norm_hk.obs.index.tolist()
            hk_cv = {}
            for j, g in enumerate(hk_genes):
                vals = hk_data[:, j]
                cv = float(np.std(vals) / (np.mean(vals) + 1e-10))
                hk_cv[g] = cv
            cv_labels = [f'{g} (CV={hk_cv[g]:.2f})' for g in hk_genes]
            fig_hk = heatmap_figure(
                hk_data.T, x_labels=hk_samples, y_labels=cv_labels,
                title='管家基因表达稳定性', x_label='Sample', y_label='Gene',
                colorbar_label='ln(CPM+1)',
                vmin=float(np.nanmin(hk_data)), vmax=float(np.nanmax(hk_data)),
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_hk, plots_dir, 'bulk_qc_housekeeping.png', 'heatmap',
                '管家基因稳定性', formats=('png', 'svg'), dpi=300,
            ))

        self.progress(90, "保存输出...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        # 移除临时分组列
        if '_auto_group' in adata_filtered.obs.columns:
            adata_filtered.obs.drop(columns=['_auto_group'], inplace=True)

        output_path = os.path.join(intermediate_dir, 'bulk_qc_output.h5ad')
        adata_raw_filtered.write_h5ad(output_path)

        removed_samples = [r['sample'] for r in filter_log_rows if not r['passed']]
        removed_reasons = {r['sample']: r['fail_reasons'] for r in filter_log_rows if not r['passed']}

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'samples_before': n_before,
                'samples_after': n_after,
                'samples_removed': n_before - n_after,
                'removed_samples': removed_samples,
                'removed_reasons': removed_reasons,
                'genes_before': genes_before_filter,
                'genes_after': adata_filtered.n_vars,
                'genes_removed': genes_before_filter - adata_filtered.n_vars,
                'outlier_samples': outlier_samples if detect_outliers else [],
                'filter_strategy': filter_strategy,
                'median_lib_size': int(np.median(adata_filtered.obs['total_counts'])),
                'input_measurement': input_measurement,
                'median_genes': int(np.median(adata_filtered.obs['n_genes_by_counts'])),
                'median_ribo_pct': round(float(np.median(adata_filtered.obs['pct_counts_ribo'])), 2) if 'pct_counts_ribo' in adata_filtered.obs.columns else 0,
                'median_gini': round(float(np.median(gini_filtered)), 4),
                'correlation_off_diagonal': corr_metadata['off_diagonal'],
                'correlation_display_range': [
                    corr_metadata['display_vmin'], corr_metadata['display_vmax'],
                ],
            }
        }


# --- 辅助函数 ---


def _gini(values):
    """计算 Gini 系数。values 为原始 count 数组（非负）。"""
    vals = np.sort(values[values >= 0])
    n = len(vals)
    if n == 0:
        return 0.0
    if np.sum(vals) == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2.0 * np.sum(index * vals) / (n * np.sum(vals))) - (n + 1) / n


def _infer_groups(sample_names):
    """从样本名推断推荐分组，优先保留多因素联合组。"""
    from modules.io_utils import infer_sample_group_candidates

    candidates = infer_sample_group_candidates(sample_names)
    if candidates:
        mapping = candidates[0]['mapping']
        return [mapping[str(name)] for name in sample_names]

    groups = []
    for name in sample_names:
        name = str(name)
        assigned = False
        for sep in ['-', '_']:
            if sep in name:
                groups.append(name.split(sep)[0])
                assigned = True
                break
        if not assigned:
            groups.append(name)
    return groups


def _detect_outliers_mahal(pca_coords, sample_names):
    """基于 PCA 坐标的马氏距离检测离群样本，返回离群样本名列表。"""
    if pca_coords.shape[0] < 4 or pca_coords.shape[1] < 2:
        return []
    n_components = min(3, pca_coords.shape[1])
    coords = pca_coords[:, :n_components]
    mean = coords.mean(axis=0)
    cov = np.cov(coords.T)
    cov_inv = np.linalg.pinv(cov)
    distances = []
    for i in range(coords.shape[0]):
        diff = coords[i] - mean
        d = np.sqrt(diff @ cov_inv @ diff)
        distances.append(d)
    distances = np.array(distances)
    med = np.median(distances)
    mad = np.median(np.abs(distances - med))
    if mad < 1e-10:
        return []
    threshold = med + 3 * 1.4826 * mad
    return [sample_names[i] for i, d in enumerate(distances) if d > threshold]
