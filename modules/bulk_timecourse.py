import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from modules.base import BaseAnalysis


class BulkTimecourseAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_timecourse"
    DISPLAY_NAME = "时序分析"
    DESCRIPTION = "多时间点差异基因检测（spline + F-test）+ 模糊 c-means 轨迹聚类"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def _build_spline_basis(self, time_vals, df=3):
        """Use patsy to create natural cubic spline basis matrix."""
        from patsy import dmatrix
        basis = dmatrix(f'bs(x, df={df}, degree=3) - 1', {'x': time_vals}, return_type='dataframe')
        return basis.values

    def _run_temporal_f_test(self, counts, time_basis, n_spline_cols):
        """For each gene fit OLS full (intercept + spline) vs reduced (intercept),
        compute F-statistic and p-value."""
        from scipy import stats
        import statsmodels.api as sm

        n_genes = counts.shape[1]
        n_obs = counts.shape[0]
        F_stats = np.zeros(n_genes)
        pvalues = np.ones(n_genes)

        X_full = sm.add_constant(time_basis)
        X_red = sm.add_constant(np.ones(n_obs))
        df_diff = n_spline_cols
        df_resid = n_obs - n_spline_cols - 1

        if df_resid <= 0:
            return F_stats, pvalues

        for i in range(n_genes):
            y = counts[:, i]
            try:
                model_full = sm.OLS(y, X_full).fit()
                model_red = sm.OLS(y, X_red).fit()
                RSS_full = model_full.ssr
                RSS_red = model_red.ssr
                if RSS_full > 0 and df_resid > 0:
                    F = ((RSS_red - RSS_full) / df_diff) / (RSS_full / df_resid)
                    F_stats[i] = max(F, 0)
                    pvalues[i] = 1 - stats.f.cdf(F, df_diff, df_resid)
            except Exception:
                F_stats[i] = 0
                pvalues[i] = 1.0

        return F_stats, pvalues

    def _run_interaction_f_test(self, counts, full_basis, n_interaction_cols, n_full_model_cols):
        """F-test for interaction terms: full (time + group + interaction) vs
        reduced (time + group only)."""
        from scipy import stats
        import statsmodels.api as sm

        n_genes = counts.shape[1]
        n_obs = counts.shape[0]
        F_stats = np.zeros(n_genes)
        pvalues = np.ones(n_genes)

        n_reduced = n_full_model_cols - n_interaction_cols
        df_diff = n_interaction_cols
        df_resid = n_obs - n_full_model_cols

        if df_resid <= 0 or n_reduced <= 0:
            return F_stats, pvalues

        X_full = sm.add_constant(full_basis)
        X_red = X_full[:, :n_reduced + 1]  # const + first n_reduced columns

        for i in range(n_genes):
            y = counts[:, i]
            try:
                model_full = sm.OLS(y, X_full).fit()
                model_red = sm.OLS(y, X_red).fit()
                RSS_full = model_full.ssr
                RSS_red = model_red.ssr
                if RSS_full > 0:
                    F = ((RSS_red - RSS_full) / df_diff) / (RSS_full / df_resid)
                    F_stats[i] = max(F, 0)
                    pvalues[i] = 1 - stats.f.cdf(F, df_diff, df_resid)
            except Exception:
                F_stats[i] = 0
                pvalues[i] = 1.0

        return F_stats, pvalues

    def run(self, input_path):
        import matplotlib.pyplot as plt
        from statsmodels.stats.multitest import multipletests
        from modules.native_figures import heatmap_figure, line_figure
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT, NATURE_GRID

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)
        from modules.io_utils import infer_expression_measurement
        from modules.io_utils import obs_grouping_info

        time_column = self.params.get('time_column', 'minute')
        requested_group_column = str(self.params.get('group_column', '') or '').strip()
        group_column = requested_group_column
        spline_df = int(self.params.get('spline_df', 3))
        n_clusters = int(self.params.get('n_clusters', 6))
        fdr_threshold = float(self.params.get('fdr_threshold', 0.05))

        if group_column:
            grouping = obs_grouping_info(
                adata, group_column, max_categories=20,
                max_numeric_categories=20, require_multiple=True,
            )
            if not grouping['valid']:
                self.progress(-1, f"交互分组列已跳过：{grouping['reason']}")
                group_column = ''

        self.progress(15, "解析时间信息...")
        if time_column not in adata.obs.columns:
            raise ValueError(
                f"时间列 '{time_column}' 不存在于 obs 中。"
                f"可用列: {', '.join(adata.obs.columns.tolist())}"
            )
        time_vals = pd.to_numeric(adata.obs[time_column], errors='coerce').values
        if np.any(np.isnan(time_vals)):
            raise ValueError(f"时间列 '{time_column}' 包含非数值，请检查数据。")

        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        counts = counts.astype(float)
        gene_ids = adata.var_names.tolist()
        if 'gene_name' in adata.var.columns:
            gene_names = adata.var['gene_name'].tolist()
        else:
            gene_names = gene_ids

        # Filter low-expression genes (mean > 1)
        self.progress(20, "过滤低表达基因...")
        raw_mean = counts.mean(axis=0)
        expr_mask = raw_mean > 1
        counts = counts[:, expr_mask]
        gene_ids = [gene_ids[i] for i in range(len(gene_ids)) if expr_mask[i]]
        gene_names = [gene_names[i] for i in range(len(gene_names)) if expr_mask[i]]
        n_genes = counts.shape[1]

        if n_genes == 0:
            raise ValueError("过滤低表达基因后无剩余基因，请降低过滤阈值或检查数据。")

        # Preserve a completed normalization; only raw counts receive CPM.
        input_measurement = infer_expression_measurement(adata, input_path)
        self.progress(25, "准备时序表达矩阵...")
        if input_measurement == 'raw_counts':
            lib_sizes = counts.sum(axis=1, keepdims=True)
            lib_sizes = np.where(lib_sizes > 0, lib_sizes, 1)
            lognorm = np.log2(counts / lib_sizes * 1e6 + 1)
        elif input_measurement == 'continuous_expression':
            lognorm = np.log2(np.maximum(counts, 0) + 1)
        else:
            lognorm = counts

        n_obs = lognorm.shape[0]
        unique_times = np.unique(time_vals)

        if len(unique_times) < 3:
            raise ValueError(
                f"需要至少 3 个不同时间点进行时序分析，当前仅 {len(unique_times)} 个。"
            )

        # Build spline basis
        self.progress(35, "构建 spline 基函数...")
        try:
            time_basis = self._build_spline_basis(time_vals, df=spline_df)
            n_spline_cols = time_basis.shape[1]
        except Exception:
            # Fallback: use raw time as single column basis
            time_basis = time_vals.reshape(-1, 1)
            n_spline_cols = 1

        # Temporal F-test
        self.progress(45, "运行时序 F-test...")
        F_stats, pvalues = self._run_temporal_f_test(lognorm, time_basis, n_spline_cols)

        # BH FDR correction
        self.progress(55, "BH FDR 校正...")
        try:
            _, qvalues, _, _ = multipletests(pvalues, method='fdr_bh')
        except Exception:
            qvalues = pvalues

        sig = (qvalues < fdr_threshold).astype(int)
        n_sig = int(sig.sum())

        # Save timecourse_results.csv
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        tc_df = pd.DataFrame({
            'gene': gene_names,
            'F': np.round(F_stats, 4),
            'pvalue': pvalues,
            'qvalue': qvalues,
            'sig': sig
        })
        tc_df = tc_df.sort_values('qvalue')
        csv_path = os.path.join(results_dir, 'timecourse_results.csv')
        tc_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '时序差异基因列表'})

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)

        # Q-Q plot: observed vs theoretical F quantiles
        self.progress(60, "生成 Q-Q 图...")
        df1 = n_spline_cols
        df2 = max(n_obs - n_spline_cols - 1, 1)
        valid_f = F_stats[np.isfinite(F_stats) & (F_stats > 0)]
        if len(valid_f) < 10:
            self.progress(-1, "有效 F 统计量不足，跳过 Q-Q 图")
        else:
            from scipy.stats import f as f_dist
            sorted_obs = np.sort(valid_f)
            n_pts = len(sorted_obs)
            theoretical_q = np.array([(i + 0.5) / (n_pts + 1) for i in range(n_pts)])
            theoretical_f = f_dist.ppf(theoretical_q, df1, df2)
            fig_qq, ax_qq = plt.subplots(figsize=(6.8, 5.0), dpi=150)
            ax_qq.scatter(theoretical_f, sorted_obs, s=12, color=NATURE_PALETTE[0],
                          alpha=0.68, linewidths=0, rasterized=True)
            max_val = max(float(theoretical_f.max()), float(sorted_obs.max())) * 1.1
            ax_qq.plot([0, max_val], [0, max_val], color=NATURE_PALETTE[3],
                       linestyle='--', linewidth=1.0)
            ax_qq.set_title('Q-Q Plot (F-statistic)', loc='left', fontsize=10,
                            fontweight='semibold', color=NATURE_TEXT)
            ax_qq.set_xlabel('Theoretical F quantiles', fontsize=9)
            ax_qq.set_ylabel('Observed F statistics', fontsize=9)
            ax_qq.grid(False)
            result_files.extend(self.save_matplotlib_figure(
                fig_qq, plots_dir, 'timecourse_qq.png', 'qq', 'Q-Q 图',
                formats=('png', 'svg'), dpi=300,
            ))

        # Fuzzy c-means trajectory clustering
        cluster_df = None
        do_cluster = n_sig >= n_clusters

        if do_cluster:
            self.progress(70, "模糊 c-means 轨迹聚类...")
            # Per-timepoint mean trajectories
            time_unique = np.sort(unique_times)
            time_idx_map = {t: np.where(time_vals == t)[0] for t in time_unique}
            n_times = len(time_unique)
            sig_idx = np.where(sig == 1)[0]
            sig_expr = lognorm[:, sig_idx]

            traj = np.zeros((len(sig_idx), n_times))
            for ti, t in enumerate(time_unique):
                t_mask = time_idx_map[t]
                traj[:, ti] = sig_expr[t_mask, :].mean(axis=0)

            # Z-score standardize per gene
            traj_mean = traj.mean(axis=1, keepdims=True)
            traj_std = traj.std(axis=1, keepdims=True)
            traj_std = np.where(traj_std > 0, traj_std, 1)
            traj_z = (traj - traj_mean) / traj_std

            n_clust = min(n_clusters, len(sig_idx), n_times)
            if n_clust < 2:
                do_cluster = False
            else:
                try:
                    from fcmeans import FCM
                    fcm = FCM(n_clusters=n_clust, random_state=42, max_iter=300)
                    fcm.fit(traj_z)
                    cluster_labels = fcm.predict(traj_z)
                    membership = fcm.u
                except Exception:
                    from sklearn.cluster import KMeans
                    km = KMeans(n_clusters=n_clust, random_state=42, n_init=10)
                    km.fit(traj_z)
                    cluster_labels = km.labels_
                    # Soft membership from inverse distances
                    dists = np.linalg.norm(
                        traj_z[:, np.newaxis, :] - km.cluster_centers_[np.newaxis, :, :],
                        axis=2
                    )
                    inv_dists = 1.0 / (dists + 1e-8)
                    membership = inv_dists / inv_dists.sum(axis=1, keepdims=True)

                sig_gene_names = [gene_names[i] for i in sig_idx]
                sig_gene_ids = [gene_ids[i] for i in sig_idx]
                cluster_df = pd.DataFrame({
                    'gene': sig_gene_names,
                    'cluster': [f'C{c}' for c in cluster_labels + 1],
                })
                # Top membership score per gene
                cluster_df['membership'] = np.round(membership.max(axis=1), 4)
                cluster_df = cluster_df.sort_values(['cluster', 'membership'], ascending=[True, False])

                tc_csv = os.path.join(results_dir, 'timecourse_clusters.csv')
                cluster_df.to_csv(tc_csv, index=False)
                result_files.append({'file_path': tc_csv, 'file_type': 'csv', 'category': 'table', 'label': '轨迹聚类结果'})

                # Cluster centers line chart
                colors = [NATURE_PALETTE[i % len(NATURE_PALETTE)] for i in range(n_clust)]
                fig_centers, ax_centers = plt.subplots(figsize=(8.0, 5.0), dpi=150)
                for ci in range(n_clust):
                    mask_c = cluster_labels == ci
                    if mask_c.sum() == 0:
                        continue
                    center = traj_z[mask_c, :].mean(axis=0)
                    ax_centers.plot(time_unique, center, color=colors[ci], linewidth=1.8,
                                    marker='o', markersize=3.5,
                                    label=f'C{ci + 1} (n={mask_c.sum()})')
                ax_centers.set_title('Cluster Centers (z-scored)', loc='left', fontsize=10,
                                     fontweight='semibold', color=NATURE_TEXT)
                ax_centers.set_xlabel(time_column, fontsize=9)
                ax_centers.set_ylabel('Z-score', fontsize=9)
                ax_centers.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
                ax_centers.legend(frameon=False, fontsize=8)
                result_files.extend(self.save_matplotlib_figure(
                    fig_centers, plots_dir, 'timecourse_cluster_centers.png',
                    'cluster_centers', '聚类中心轨迹', formats=('png', 'svg'), dpi=300,
                ))

                # Gene x Time heatmap ordered by cluster
                order = np.argsort(cluster_labels)
                heatmap_z = traj_z[order, :]
                y_labels_cluster = [f"C{cluster_labels[i] + 1}_{sig_gene_names[i]}" for i in order]
                # Cap displayed genes at 200 for readability
                max_heat = 200
                if len(order) > max_heat:
                    step = len(order) // max_heat
                    sel = np.arange(0, len(order), step)[:max_heat]
                    heatmap_z = heatmap_z[sel, :]
                    y_labels_cluster = [y_labels_cluster[j] for j in sel]

                finite_heat = heatmap_z[np.isfinite(heatmap_z)]
                max_abs = float(np.nanmax(np.abs(finite_heat))) if finite_heat.size else 1.0
                fig_heat = heatmap_figure(
                    heatmap_z, x_labels=[str(t) for t in time_unique],
                    y_labels=y_labels_cluster,
                    title='Gene x Time Heatmap (ordered by cluster)',
                    x_label=time_column, y_label='Gene', colorbar_label='Z-score',
                    vmin=-max_abs, vmax=max_abs,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_heat, plots_dir, 'timecourse_heatmap.png', 'heatmap',
                    '基因x时间热图', formats=('png', 'svg'), dpi=300,
                ))

        # Pairwise group comparison per timepoint
        pairwise_groups_str = self.params.get('pairwise_groups', '').strip()
        n_pairwise_sig = 0
        if pairwise_groups_str and group_column and group_column in adata.obs.columns:
            pw_pair = None
            if '-vs-' in pairwise_groups_str:
                parts = pairwise_groups_str.split('-vs-')
                if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                    pw_pair = (parts[0].strip(), parts[1].strip())

            if pw_pair:
                self.progress(80, f"配对比较: {pw_pair[0]} vs {pw_pair[1]}...")
                pw_a, pw_b = pw_pair
                groups_series = adata.obs[group_column].astype(str)
                pairwise_rows = []

                for t in time_unique:
                    t_mask = time_vals == t
                    a_mask = t_mask & (groups_series == pw_a)
                    b_mask = t_mask & (groups_series == pw_b)
                    a_idx = np.where(a_mask)[0]
                    b_idx = np.where(b_mask)[0]

                    if len(a_idx) < 2 or len(b_idx) < 2:
                        continue

                    for gi in range(n_genes):
                        vals_a = lognorm[a_idx, gi]
                        vals_b = lognorm[b_idx, gi]
                        try:
                            tstat, pval = ttest_ind(vals_a, vals_b, equal_var=False)
                        except Exception:
                            tstat, pval = 0.0, 1.0
                        mean_a = vals_a.mean()
                        mean_b = vals_b.mean()
                        l2fc = np.log2((mean_b + 1) / (mean_a + 1))
                        pairwise_rows.append({
                            'gene': gene_names[gi],
                            'time': t,
                            'tstat': round(float(tstat), 4),
                            'pvalue': float(pval),
                            'log2FC': round(float(l2fc), 4),
                            'mean_a': round(float(mean_a), 4),
                            'mean_b': round(float(mean_b), 4),
                        })

                if pairwise_rows:
                    pw_df = pd.DataFrame(pairwise_rows)

                    # Per-timepoint BH FDR correction
                    for t in time_unique:
                        t_mask = pw_df['time'] == t
                        if t_mask.sum() > 0:
                            try:
                                _, qvals, _, _ = multipletests(
                                    pw_df.loc[t_mask, 'pvalue'].values, method='fdr_bh')
                                pw_df.loc[t_mask, 'qvalue'] = qvals
                            except Exception:
                                pw_df.loc[t_mask, 'qvalue'] = pw_df.loc[t_mask, 'pvalue']

                    pw_df = pw_df.sort_values('qvalue')

                    safe_name = f'{pw_a}_vs_{pw_b}'.replace(' ', '_')
                    pw_csv = os.path.join(results_dir, f'timecourse_pairwise_{safe_name}.csv')
                    pw_df.to_csv(pw_csv, index=False)
                    result_files.append({
                        'file_path': pw_csv, 'file_type': 'csv',
                        'category': 'table',
                        'label': f'配对比较结果 ({pw_a} vs {pw_b})'
                    })

                    # Heatmap: -log10(qvalue)
                    sig_pw = pw_df[pw_df['qvalue'] < fdr_threshold]
                    pw_genes_ordered = sig_pw['gene'].unique().tolist()
                    if not pw_genes_ordered:
                        pw_genes_ordered = pw_df.groupby('gene')['qvalue'].min().nsmallest(50).index.tolist()

                    time_list = sorted(pw_df['time'].unique())
                    gene_subset = pw_genes_ordered[:200]
                    pw_pivot = pw_df[pw_df['gene'].isin(gene_subset)].pivot_table(
                        index='gene', columns='time', values='qvalue', fill_value=1.0)
                    pw_pivot = pw_pivot.reindex(index=gene_subset, columns=time_list, fill_value=1.0)

                    neg_log_q = -np.log10(pw_pivot.values + 1e-300)
                    fig_pw = heatmap_figure(
                        neg_log_q, x_labels=[str(t) for t in time_list],
                        y_labels=pw_pivot.index.tolist(),
                        title=f'Pairwise -log10(q): {pw_a} vs {pw_b}',
                        x_label=time_column, y_label='Gene', colorbar_label='-log10(q)',
                        vmin=0, vmax=max(1.0, float(np.nanmax(neg_log_q)) if neg_log_q.size else 1.0),
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_pw, plots_dir, f'timecourse_pairwise_{safe_name}_heatmap.png',
                        'heatmap', f'配对比较热图 ({pw_a} vs {pw_b})',
                        formats=('png', 'svg'), dpi=300,
                    ))

                    # Heatmap: log2FC
                    pw_pivot_fc = pw_df[pw_df['gene'].isin(gene_subset)].pivot_table(
                        index='gene', columns='time', values='log2FC', fill_value=0.0)
                    pw_pivot_fc = pw_pivot_fc.reindex(index=gene_subset, columns=time_list, fill_value=0.0)

                    max_fc = float(np.nanmax(np.abs(pw_pivot_fc.values))) if pw_pivot_fc.size else 1.0
                    fig_pw_fc = heatmap_figure(
                        pw_pivot_fc.values, x_labels=[str(t) for t in time_list],
                        y_labels=pw_pivot_fc.index.tolist(),
                        title=f'Pairwise log2FC: {pw_a} vs {pw_b}',
                        x_label=time_column, y_label='Gene', colorbar_label='log2FC',
                        vmin=-max(1.0, max_fc), vmax=max(1.0, max_fc),
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_pw_fc, plots_dir, f'timecourse_pairwise_{safe_name}_log2fc.png',
                        'heatmap', f'配对比较 log2FC ({pw_a} vs {pw_b})',
                        formats=('png', 'svg'), dpi=300,
                    ))

                    n_pairwise_sig = int((pw_df['qvalue'] < fdr_threshold).sum())
                    del pw_df, pairwise_rows

        # Interaction test (if group_column specified)
        self.progress(85, "交互效应分析...")
        interaction_csv = None
        if group_column and group_column in adata.obs.columns:
            groups = adata.obs[group_column].astype(str)
            group_dummies = pd.get_dummies(groups, drop_first=True, dtype=float).values
            n_group_cols = group_dummies.shape[1]
            interaction = time_basis[:, :, np.newaxis] * group_dummies[:, np.newaxis, :]
            n_interaction_cols = time_basis.shape[1] * n_group_cols
            interaction_2d = interaction.reshape(n_obs, -1)
            full_basis = np.hstack([time_basis, group_dummies, interaction_2d])
            n_full_model_cols = time_basis.shape[1] + n_group_cols + n_interaction_cols

            F_int, p_int = self._run_interaction_f_test(
                lognorm, full_basis, n_interaction_cols, n_full_model_cols
            )
            try:
                _, q_int, _, _ = multipletests(p_int, method='fdr_bh')
            except Exception:
                q_int = p_int

            int_df = pd.DataFrame({
                'gene': gene_names,
                'F': np.round(F_int, 4),
                'pvalue': p_int,
                'qvalue': q_int,
                'sig': (q_int < fdr_threshold).astype(int)
            })
            int_df = int_df.sort_values('qvalue')
            interaction_csv = os.path.join(results_dir, 'timecourse_interaction.csv')
            int_df.to_csv(interaction_csv, index=False)
            result_files.append({
                'file_path': interaction_csv, 'file_type': 'csv',
                'category': 'table', 'label': '交互效应检验结果'
            })

        # Save output h5ad
        self.progress(92, "保存 h5ad...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_timecourse_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_genes_total': n_genes,
                'n_temporal_sig': n_sig,
                'fdr_threshold': fdr_threshold,
                'spline_df': spline_df,
                'n_timepoints': len(unique_times),
                'n_clusters': n_clusters if do_cluster else 0,
                'has_interaction': interaction_csv is not None,
                'n_pairwise_sig': n_pairwise_sig,
                'requested_group_column': requested_group_column,
                'group_column': group_column,
            }
        }
