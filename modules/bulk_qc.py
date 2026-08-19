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
        from figure_engine import NatureFigureDirector, export_registered_figure
        director = NatureFigureDirector()
        nature_formats = ('svg', 'pdf', 'png')

        def _export_diagnostic(fig, stem, label, *, width='single', height_mm=None,
                               category='qc', plot_type='diagnostic'):
            spec = director.spec_from_params(
                plot_type, self.params, width=width, title=label,
            ).with_updates(formats=nature_formats, height_mm=height_mm)
            exported, report = export_registered_figure(
                fig, os.path.join(plots_dir, stem), spec,
                category=category, label=label,
                qa_path=os.path.join(results_dir, f'{stem}_nature_readiness.json'),
            )
            if not report.ready:
                self.progress(-1, f'{label} Nature readiness {report.score}/100；请查看 QA 报告。')
            plt.close(fig)
            return exported

        overview = [
            {'values': lib_sizes, 'title': 'Library size', 'ylabel': 'Library size'},
            {'values': n_genes_detected, 'title': 'Detected genes', 'ylabel': 'Detected genes'},
            {'values': mt_pct, 'title': 'Mitochondrial fraction', 'ylabel': 'MT%'},
            {'values': ribo_pct, 'title': 'Ribosomal fraction', 'ylabel': 'Ribo%'},
            {'values': gini_values, 'title': 'Gini coefficient', 'ylabel': 'Gini'},
        ]
        invariant = []
        for item in overview:
            finite_values = np.asarray(item['values'], dtype=float)
            finite_values = finite_values[np.isfinite(finite_values)]
            if finite_values.size == 0 or np.nanmax(finite_values) - np.nanmin(finite_values) <= 1e-12:
                invariant.append(item['title'])
        result_files.extend(_export_diagnostic(
            director.render(
                director.spec_from_params('diagnostic', self.params, width='double',
                                          title='Bulk RNA-seq QC overview').with_updates(
                                              extra={'kind': 'qc_overview'},
                                              formats=nature_formats, height_mm=130.0),
                {
                    'kind': 'qc_overview',
                    'sample_labels': sample_names,
                    'metrics': overview,
                    'library_size': lib_sizes,
                    'detected_genes': n_genes_detected,
                    'pass_colors': [('#4C78A8' if m else '#77808C') for m in mask],
                    'omitted_metrics': ', '.join(invariant),
                },
            ),
            'bulk_qc_overview', '质控总览', width='double', height_mm=130.0,
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
        from figure_engine import NatureFigureDirector, export_registered_figure
        director = NatureFigureDirector()
        corr_spec = director.spec_from_params(
            'correlation', self.params, width='double',
            title='样本相关性热图 (Pearson)', correlation_method='pearson',
        )
        fig_corr = director.render(corr_spec, {
            'correlation_matrix': corr_matrix,
            'sample_labels': sample_labels_corr,
            'groups': filtered_groups,
        })
        corr_metadata = fig_corr._nature_correlation_metadata
        exported, corr_readiness = export_registered_figure(
            fig_corr, os.path.join(plots_dir, 'bulk_qc_corr'), corr_spec,
            category='heatmap', label='样本相关性热图',
            qa_path=os.path.join(results_dir, 'bulk_qc_corr_nature_readiness.json'),
        )
        result_files.extend(exported)
        if not corr_readiness.ready:
            self.progress(-1, f'Correlation Nature readiness {corr_readiness.score}/100；请查看 QA 报告。')
        plt.close(fig_corr)
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
            fig_r = director.render(
                director.spec_from_params('diagnostic', self.params, width='single',
                                          title='Sample filtering').with_updates(
                                              extra={'kind': 'normalization_library'},
                                              formats=nature_formats, height_mm=60.0),
                {'kind': 'normalization_library', 'sample_labels': ['Before', 'After'],
                 'raw_values': [n_before, n_before], 'normalized_values': [n_after, n_after],
                 'raw_ylabel': 'Samples', 'normalized_ylabel': 'Samples'},
            )
            result_files.extend(_export_diagnostic(
                fig_r, 'bulk_qc_filter', '样本过滤结果', width='single', height_mm=60.0,
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

            pca_variance_qc = adata_normed.uns.get('pca', {}).get('variance_ratio', [])
            pca_spec = director.spec_from_params(
                'pca', self.params, width='single', title='QC-filtered sample PCA',
                show_legend=len(unique_groups_pca) <= 8,
            ).with_updates(formats=nature_formats, height_mm=82.0,
                           outlier_labels=tuple(outlier_samples))
            # Many QC projects are auto-grouped by sample ID, producing one
            # category per replicate.  Preserve that information with marker
            # shape while using one stable publication colour; this avoids a
            # misleading rainbow legend and color-vision collisions.
            pca_group_colors = ({group: NATURE_PALETTE[0] for group in unique_groups_pca}
                                if len(unique_groups_pca) > 8 else None)
            fig_pca = director.render(pca_spec, {
                'coordinates': pc[:, :2],
                'groups': filtered_groups_pca,
                'batches': filtered_groups_pca if len(unique_groups_pca) > 8 else None,
                'group_colors': pca_group_colors,
                'samples': filtered_sample_names,
                'explained_variance': pca_variance_qc[:2],
            })
            result_files.extend(_export_diagnostic(
                fig_pca, 'bulk_qc_pca', '样本 PCA', width='single', height_mm=82.0, category='pca',
                plot_type='pca',
            ))

        # PCA 方差解释 elbow 图（使用 scanpy 存储的 variance_ratio）
        pca_variance = adata_normed.uns.get('pca', {}).get('variance_ratio', None)
        if pca_variance is None and n_comps >= 2:
            pca_var = np.var(pc, axis=0)
            pca_variance = pca_var / pca_var.sum()
        if pca_variance is not None and len(pca_variance) > 0:
            n_pcs = len(pca_variance)
            pc_labels = [f'PC{i+1}' for i in range(n_pcs)]
            fig_elbow = director.render(
                director.spec_from_params('diagnostic', self.params, width='single',
                                          title='PCA variance explained').with_updates(
                                              extra={'kind': 'variance'},
                                              formats=nature_formats, height_mm=70.0),
                {'kind': 'variance', 'variance': pca_variance, 'labels': pc_labels},
            )
            result_files.extend(_export_diagnostic(
                fig_elbow, 'bulk_qc_pca_elbow', 'PCA 方差解释', width='single', height_mm=70.0,
                category='pca',
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
            box_data = [values for values in (intra_dists, inter_dists) if values]
            box_labels = [label for label, values in zip(['组内距离', '组间距离'],
                                                          (intra_dists, inter_dists)) if values]
            fig_dist = director.render(
                director.spec_from_params('diagnostic', self.params, width='single',
                                          title=f'Within vs between-group distance{title_suffix}').with_updates(
                                              extra={'kind': 'boxplot'},
                                              formats=nature_formats, height_mm=68.0),
                {'kind': 'boxplot', 'groups': box_labels, 'values': box_data,
                 'ylabel': '1 - Pearson r'},
            )
            result_files.extend(_export_diagnostic(
                fig_dist, 'bulk_qc_group_distance', '组内/组间距离', width='single', height_mm=68.0,
            ))

        # QC 指标散点矩阵 (Pairs Plot)
        obs_filtered = adata_filtered.obs
        pairs_groups = [groups[sample_names.index(s)] for s in obs_filtered.index.tolist()]
        unique_pg = sorted(set(pairs_groups))
        pg_color_map = ({g: NATURE_PALETTE[0] for g in unique_pg}
                        if len(unique_pg) > 3 else
                        {g: NATURE_PALETTE[i % len(NATURE_PALETTE)]
                         for i, g in enumerate(unique_pg)})
        from figure_engine.style import get_style
        pg_marker_map = get_style('nature').batch_markers(unique_pg)
        pg_colors = [pg_color_map[g] for g in pairs_groups]
        pair_values = [
            obs_filtered['total_counts'].values,
            obs_filtered['n_genes_by_counts'].values,
            obs_filtered['pct_counts_mt'].values if 'pct_counts_mt' in obs_filtered.columns else np.zeros(n_after),
            obs_filtered['pct_counts_ribo'].values if 'pct_counts_ribo' in obs_filtered.columns else np.zeros(n_after),
        ]
        pair_labels = ['Library Size', 'N Genes', 'MT%', 'Ribo%']
        fig_pairs = director.render(
            director.spec_from_params('diagnostic', self.params, width='double',
                                      title='QC metric pairs').with_updates(
                                          extra={'kind': 'pairs'}, formats=nature_formats,
                                          height_mm=112.0),
            {'kind': 'pairs', 'labels': pair_labels, 'values': pair_values,
             'groups': pairs_groups, 'group_colors': pg_color_map,
             'group_markers': pg_marker_map},
        )
        result_files.extend(_export_diagnostic(
            fig_pairs, 'bulk_qc_pairs_plot', 'QC 指标散点矩阵', width='double', height_mm=112.0,
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
            violin_metrics = []
            for metric in ('MT%', 'Ribo%', 'Library Size', 'N Genes'):
                vals_by_group = [violin_df[violin_df['group'] == g][metric].values
                                 for g in unique_groups]
                if any((lambda finite: finite.size > 0 and
                        np.nanmax(finite) - np.nanmin(finite) > 1e-12)
                       (np.asarray(v, dtype=float)[np.isfinite(v)]) for v in vals_by_group):
                    violin_metrics.append({'label': metric, 'values': vals_by_group})
            if violin_metrics:
                fig_violin = director.render(
                    director.spec_from_params('diagnostic', self.params, width='double',
                                              title='QC metrics by group').with_updates(
                                                  extra={'kind': 'violin'}, formats=nature_formats,
                                                  height_mm=90.0),
                    {'kind': 'violin', 'groups': unique_groups, 'metrics': violin_metrics},
                )
                result_files.extend(_export_diagnostic(
                    fig_violin, 'bulk_qc_violin_by_group', '各组 QC 指标分布', width='double', height_mm=90.0,
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
            hk_spec = director.spec_from_params(
                'heatmap', self.params, width='double', title='Housekeeping gene stability',
                zscore='row', row_cluster=True, col_cluster=True,
            ).with_updates(formats=nature_formats, height_mm=100.0, max_row_labels=20,
                           max_col_labels=18)
            fig_hk = director.render(hk_spec, {
                'matrix': hk_data.T,
                'gene_labels': cv_labels,
                'sample_labels': hk_samples,
                # Auto-grouping can create one legend entry per sample.  A
                # long annotation legend is not informative and collides with
                # the heatmap/colorbar at publication size, so only retain it
                # when it remains compact.
                'annotations': ({'Group': [groups[sample_names.index(s)] for s in hk_samples]}
                                if len(set(groups[sample_names.index(s)] for s in hk_samples)) <= 6
                                else {}),
                'colorbar_label': 'Gene-wise z-score',
            })
            result_files.extend(_export_diagnostic(
                fig_hk, 'bulk_qc_housekeeping', '管家基因稳定性', width='double', height_mm=100.0,
                category='heatmap', plot_type='heatmap',
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
