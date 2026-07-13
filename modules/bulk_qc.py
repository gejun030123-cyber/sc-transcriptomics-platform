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
        from modules.visualization import scatter_plot, bar_plot
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

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

        fig = make_subplots(rows=3, cols=2,
            subplot_titles=['文库大小分布', '检测基因数',
                           '线粒体基因比例', '核糖体基因比例',
                           'Gini 系数', '文库大小 vs 检测基因数'],
            vertical_spacing=0.08)
        sample_idx = list(range(n_before))
        fig.add_trace(go.Bar(x=sample_idx, y=lib_sizes.tolist(), marker_color='#1a237e', name='文库大小'), row=1, col=1)
        fig.add_trace(go.Bar(x=sample_idx, y=n_genes_detected.tolist(), marker_color='#283593', name='基因数'), row=1, col=2)
        fig.add_trace(go.Bar(x=sample_idx, y=mt_pct.tolist(), marker_color='#e53935', name='MT%'), row=2, col=1)
        fig.add_trace(go.Bar(x=sample_idx, y=ribo_pct.tolist(), marker_color='#ff8f00', name='Ribo%'), row=2, col=2)
        fig.add_trace(go.Bar(x=sample_idx, y=gini_values.tolist(), marker_color='#6a1b9a', name='Gini'), row=3, col=1)
        colors = ['#4caf50' if m else '#e53935' for m in mask]
        fig.add_trace(go.Scattergl(x=lib_sizes.tolist(), y=n_genes_detected.tolist(), mode='markers',
            marker=dict(color=colors, size=6), name='样本'), row=3, col=2)
        fig.update_layout(height=900, width=800, showlegend=False, title='Bulk RNA-seq 质控总览')
        fpath = os.path.join(plots_dir, 'bulk_qc_overview.json')
        with open(fpath, 'w') as f: f.write(fig.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '质控总览'})

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

        # 分组颜色
        filtered_groups = [groups[sample_names.index(s)] for s in sample_labels_corr]
        unique_groups = sorted(set(filtered_groups))
        group_color_map = {g: f'hsl({i*360//max(1,len(unique_groups))},70%,50%)' for i, g in enumerate(unique_groups)}
        group_colors = [group_color_map[g] for g in filtered_groups]

        fig_corr = go.Figure()
        fig_corr.add_trace(go.Heatmap(
            z=corr_matrix.tolist(), x=sample_labels_corr, y=sample_labels_corr,
            colorscale='Blues', zmin=0, zmax=1,
            colorbar=dict(title='Pearson r'),
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
        ))
        # 分组 annotation bar (y-axis color strip via annotations)
        for i, (sample, color) in enumerate(zip(sample_labels_corr, group_colors)):
            fig_corr.add_annotation(
                x=-0.02, y=i, xref='paper', yref='y',
                text='', showarrow=False,
                xanchor='right',
                bgcolor=color, bordercolor=color, width=12, height=12)
        fig_corr.update_layout(
            title='样本相关性热图 (Pearson)',
            height=max(400, n_after * 30 + 100), width=max(500, n_after * 30 + 200),
            plot_bgcolor='white',
            margin=dict(l=80)
        )
        fpath = os.path.join(plots_dir, 'bulk_qc_corr.json')
        with open(fpath, 'w') as f: f.write(fig_corr.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '样本相关性热图'})

        if n_before > n_after:
            fig_r = go.Figure()
            fig_r.add_trace(go.Bar(x=['过滤前', '过滤后'], y=[n_before, n_after], marker_color=['#e53935', '#4caf50']))
            fig_r.update_layout(title='样本过滤结果', yaxis_title='样本数', width=400, height=300)
            fpath = os.path.join(plots_dir, 'bulk_qc_filter.json')
            with open(fpath, 'w') as f: f.write(fig_r.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '样本过滤结果'})

        self.progress(75, "PCA 离群检测...")
        outlier_samples = []
        n_comps = min(10, n_after - 1, adata_normed.n_vars - 1)
        if n_comps >= 2:
            sc.pp.pca(adata_normed, n_comps=n_comps)
            pc = adata_normed.obsm['X_pca']
            filtered_sample_names = adata_normed.obs.index.tolist()
            filtered_groups_pca = [groups[sample_names.index(s)] for s in filtered_sample_names]
            unique_groups_pca = sorted(set(filtered_groups_pca))
            group_color_map_pca = {g: f'hsl({i*360//max(1,len(unique_groups_pca))},70%,50%)' for i, g in enumerate(unique_groups_pca)}

            # 离群检测
            if detect_outliers:
                outlier_samples = _detect_outliers_mahal(pc, filtered_sample_names)

            fig_pca = go.Figure()
            for g in unique_groups_pca:
                grp_idx = [i for i in range(n_after) if filtered_groups_pca[i] == g]
                fig_pca.add_trace(go.Scattergl(
                    x=pc[grp_idx, 0].tolist(), y=pc[grp_idx, 1].tolist(),
                    mode='markers+text',
                    text=[filtered_sample_names[i] for i in grp_idx],
                    textposition='top center',
                    marker=dict(size=8, color=group_color_map_pca[g]),
                    name=g))
            # 离群点高亮
            if outlier_samples:
                out_idx = [filtered_sample_names.index(s) for s in outlier_samples if s in filtered_sample_names]
                if out_idx:
                    fig_pca.add_trace(go.Scattergl(
                        x=pc[out_idx, 0].tolist(), y=pc[out_idx, 1].tolist(),
                        mode='markers',
                        marker=dict(size=14, color='red', symbol='x', line=dict(width=2, color='darkred')),
                        name='离群样本'))
            fig_pca.update_layout(title='质控后样本 PCA', xaxis_title='PC1', yaxis_title='PC2',
                                 plot_bgcolor='white', width=600, height=500)
            fpath = os.path.join(plots_dir, 'bulk_qc_pca.json')
            with open(fpath, 'w') as f: f.write(fig_pca.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': '样本 PCA'})

        # PCA 方差解释 elbow 图（使用 scanpy 存储的 variance_ratio）
        pca_variance = adata_normed.uns.get('pca', {}).get('variance_ratio', None)
        if pca_variance is None and n_comps >= 2:
            pca_var = np.var(pc, axis=0)
            pca_variance = pca_var / pca_var.sum()
        if pca_variance is not None and len(pca_variance) > 0:
            n_pcs = len(pca_variance)
            pc_labels = [f'PC{i+1}' for i in range(n_pcs)]
            cumulative = np.cumsum(pca_variance).tolist()
            fig_elbow = go.Figure()
            fig_elbow.add_trace(go.Bar(x=pc_labels, y=pca_variance.tolist(),
                marker_color='#1a237e', name='方差比例'))
            fig_elbow.add_trace(go.Scatter(x=pc_labels, y=cumulative,
                mode='lines+markers', marker_color='#e53935', name='累积比例', yaxis='y2'))
            fig_elbow.update_layout(
                title='PCA 方差解释比例',
                xaxis_title='主成分', yaxis_title='方差解释比例',
                yaxis2=dict(title='累积比例', overlaying='y', side='right', range=[0, 1.05]),
                width=600, height=400, plot_bgcolor='white')
            fpath = os.path.join(plots_dir, 'bulk_qc_pca_elbow.json')
            with open(fpath, 'w') as f: f.write(fig_elbow.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': 'PCA 方差解释'})

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

            fig_dist = go.Figure()
            if intra_dists:
                fig_dist.add_trace(go.Box(y=intra_dists, name='组内距离', marker_color='#4caf50'))
            if inter_dists:
                fig_dist.add_trace(go.Box(y=inter_dists, name='组间距离', marker_color='#e53935'))
            try:
                from scipy.stats import ttest_ind
                _, pval = ttest_ind(intra_dists, inter_dists, equal_var=False)
                title_suffix = f' (p={pval:.2e})'
            except Exception:
                title_suffix = ''
            fig_dist.update_layout(
                title=f'组内 vs 组间距离{title_suffix}',
                yaxis_title='1 - Pearson r', width=500, height=400, plot_bgcolor='white')
            fpath = os.path.join(plots_dir, 'bulk_qc_group_distance.json')
            with open(fpath, 'w') as f: f.write(fig_dist.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '组内/组间距离'})

        # QC 指标散点矩阵 (Pairs Plot)
        obs_filtered = adata_filtered.obs
        pairs_groups = [groups[sample_names.index(s)] for s in obs_filtered.index.tolist()]
        unique_pg = sorted(set(pairs_groups))
        pg_color_map = {g: f'hsl({i*360//max(1,len(unique_pg))},70%,50%)' for i, g in enumerate(unique_pg)}
        pg_colors = [pg_color_map[g] for g in pairs_groups]

        fig_pairs = go.Figure(data=go.Splom(
            dimensions=[
                dict(label='Library Size', values=obs_filtered['total_counts'].tolist()),
                dict(label='N Genes', values=obs_filtered['n_genes_by_counts'].tolist()),
                dict(label='MT%', values=obs_filtered['pct_counts_mt'].tolist() if 'pct_counts_mt' in obs_filtered.columns else [0]*n_after),
                dict(label='Ribo%', values=obs_filtered['pct_counts_ribo'].tolist() if 'pct_counts_ribo' in obs_filtered.columns else [0]*n_after),
            ],
            marker=dict(color=pg_colors, size=5, line=dict(width=0.5, color='white')),
            text=obs_filtered.index.tolist(),
            showupperhalf=False,
        ))
        fig_pairs.update_layout(title='QC 指标散点矩阵', width=700, height=700)
        fpath = os.path.join(plots_dir, 'bulk_qc_pairs_plot.json')
        with open(fpath, 'w') as f: f.write(fig_pairs.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': 'QC 指标散点矩阵'})

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
            fig_violin = make_subplots(rows=2, cols=2,
                subplot_titles=['MT%', 'Ribo%', 'Library Size', 'N Genes'],
                shared_xaxes=True)
            metrics = [('MT%', 1, 1), ('Ribo%', 1, 2), ('Library Size', 2, 1), ('N Genes', 2, 2)]
            for metric, row, col in metrics:
                for g in unique_groups:
                    vals = violin_df[violin_df['group'] == g][metric].tolist()
                    fig_violin.add_trace(go.Violin(
                        y=vals, name=g, box_visible=True, meanline_visible=True,
                        legendgroup=g, showlegend=(row == 1 and col == 1)),
                        row=row, col=col)
            fig_violin.update_layout(title='各组 QC 指标分布', height=600, width=700)
            fpath = os.path.join(plots_dir, 'bulk_qc_violin_by_group.json')
            with open(fpath, 'w') as f: f.write(fig_violin.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '各组 QC 指标分布'})

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
            fig_hk = go.Figure(data=go.Heatmap(
                z=hk_data.T.tolist(), x=hk_samples, y=cv_labels,
                colorscale='YlOrRd', colorbar=dict(title='ln(CPM+1)')))
            fig_hk.update_layout(title='管家基因表达稳定性',
                width=max(400, len(hk_samples)*40+200),
                height=max(200, len(hk_genes)*30+100))
            fpath = os.path.join(plots_dir, 'bulk_qc_housekeeping.json')
            with open(fpath, 'w') as f: f.write(fig_hk.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '管家基因稳定性'})

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
