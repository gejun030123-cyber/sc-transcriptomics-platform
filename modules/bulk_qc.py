import os
import json
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

        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'bulk_qc')

        min_counts = int(self.params.get('min_counts', 100000))
        min_genes = int(self.params.get('min_genes', 5000))
        max_mt_pct = float(self.params.get('max_mt_pct', 20.0))
        max_ribo_pct = float(self.params.get('max_ribo_pct', 40.0))
        min_gini = float(self.params.get('min_gini', 0))
        min_sample_expr = int(self.params.get('min_sample_expr', 0))
        detect_outliers = self.params.get('detect_outliers', True)
        filter_strategy = self.params.get('filter_strategy', 'standard')

        # 过滤策略预设覆盖阈值
        if filter_strategy == 'conservative':
            min_counts = min(min_counts, 50000)
            min_genes = min(min_genes, 3000)
            max_mt_pct = max(max_mt_pct, 30.0)
            max_ribo_pct = max(max_ribo_pct, 60.0)

        self.progress(20, "计算质控指标...")
        # 优先用 gene_name 检测线粒体基因（Ensembl ID 不以 MT- 开头）
        if 'gene_name' in adata.var.columns:
            gene_names_for_mt = adata.var['gene_name'].fillna('').astype(str)
        else:
            gene_names_for_mt = adata.var_names.astype(str)
        adata.var['mt'] = gene_names_for_mt.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # 核糖体基因检测
        if 'gene_name' in adata.var.columns:
            gene_names_for_ribo = adata.var['gene_name'].fillna('').astype(str)
        else:
            gene_names_for_ribo = adata.var_names.astype(str)
        adata.var['ribo'] = gene_names_for_ribo.str.startswith(('RPL', 'RPS'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)

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
        if group_col and group_col in adata.obs.columns:
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
        n_after = adata_filtered.n_obs

        # 输出样本指标表和过滤日志
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        metrics_df = pd.DataFrame(filter_log_rows)
        metrics_csv = os.path.join(results_dir, 'bulk_qc_sample_metrics.csv')
        metrics_df.to_csv(metrics_csv, index=False)

        filter_log_csv = os.path.join(results_dir, 'bulk_qc_filter_log.csv')
        metrics_df.to_csv(filter_log_csv, index=False)

        result_files.append({'file_path': metrics_csv, 'file_type': 'csv', 'category': 'table', 'label': '样本 QC 指标'})
        result_files.append({'file_path': filter_log_csv, 'file_type': 'csv', 'category': 'table', 'label': '过滤日志'})

        # 基因层面过滤
        genes_before_filter = adata_filtered.n_vars
        gene_filter_rows = []
        if min_sample_expr > 0:
            raw_filt = adata_filtered.X.toarray() if hasattr(adata_filtered.X, 'toarray') else np.asarray(adata_filtered.X)
            lib_sizes_filt = raw_filt.sum(axis=1, keepdims=True)
            lib_sizes_filt[lib_sizes_filt == 0] = 1
            cpm = raw_filt / lib_sizes_filt * 1e6
            expr_count_per_gene = (cpm > 1).sum(axis=0)
            gene_mask = expr_count_per_gene >= min_sample_expr

            removed_gene_names = adata_filtered.var_names[~gene_mask].tolist()
            for g in removed_gene_names:
                idx = list(adata_filtered.var_names).index(g)
                gene_filter_rows.append({
                    'gene': g,
                    'expressed_in_n_samples': int(expr_count_per_gene[idx]),
                    'action': 'removed',
                })
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

        fig = make_subplots(rows=2, cols=2,
            subplot_titles=['文库大小分布', '检测基因数',
                           '线粒体基因比例', '文库大小 vs 检测基因数'])
        fig.add_trace(go.Bar(x=list(range(n_before)), y=lib_sizes, marker_color='#1a237e', name='文库大小'), row=1, col=1)
        fig.add_trace(go.Bar(x=list(range(n_before)), y=n_genes_detected, marker_color='#283593', name='基因数'), row=1, col=2)
        fig.add_trace(go.Bar(x=list(range(n_before)), y=mt_pct, marker_color='#e53935', name='MT%'), row=2, col=1)
        colors = ['#4caf50' if m else '#e53935' for m in mask]
        fig.add_trace(go.Scattergl(x=lib_sizes, y=n_genes_detected, mode='markers',
            marker=dict(color=colors, size=6), name='样本'), row=2, col=2)
        fig.update_layout(height=600, width=800, showlegend=False, title='Bulk RNA-seq 质控总览')
        fpath = os.path.join(plots_dir, 'bulk_qc_overview.json')
        with open(fpath, 'w') as f: f.write(fig.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '质控总览'})

        # 样本相关性热图
        norm_for_corr = adata_filtered.copy()
        sc.pp.normalize_total(norm_for_corr, target_sum=1e6)
        sc.pp.log1p(norm_for_corr)
        corr_data = norm_for_corr.X if not hasattr(norm_for_corr.X, 'toarray') else norm_for_corr.X.toarray()
        corr_matrix = np.corrcoef(corr_data)
        sample_labels_corr = norm_for_corr.obs.index.tolist()

        fig_corr = go.Figure()
        fig_corr.add_trace(go.Heatmap(
            z=corr_matrix.tolist(), x=sample_labels_corr, y=sample_labels_corr,
            colorscale='Blues', zmin=0, zmax=1,
            colorbar=dict(title='Pearson r'),
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
        ))
        fig_corr.update_layout(
            title='样本相关性热图 (Pearson)',
            height=max(400, n_after * 30 + 100), width=max(400, n_after * 30 + 100),
            plot_bgcolor='white'
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

        self.progress(80, "运行 PCA 离群检测...")
        sc.pp.normalize_total(adata_filtered, target_sum=1e6)
        sc.pp.log1p(adata_filtered)
        sc.pp.pca(adata_filtered, n_comps=min(10, n_after - 1))
        fig_pca = go.Figure()
        pc = adata_filtered.obsm['X_pca']
        fig_pca.add_trace(go.Scattergl(x=pc[:, 0], y=pc[:, 1], mode='markers+text',
            text=adata_filtered.obs.index.tolist(), textposition='top center',
            marker=dict(size=8, color='#1a237e')))
        fig_pca.update_layout(title='质控后样本 PCA', xaxis_title='PC1', yaxis_title='PC2',
                             plot_bgcolor='white', width=600, height=500)
        fpath = os.path.join(plots_dir, 'bulk_qc_pca.json')
        with open(fpath, 'w') as f: f.write(fig_pca.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': '样本 PCA'})

        self.progress(90, "保存输出...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_qc_output.h5ad')
        adata_filtered.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'samples_before': n_before,
                'samples_after': n_after,
                'samples_removed': n_before - n_after,
                'genes_total': adata_filtered.n_vars,
                'median_lib_size': int(np.median(adata_filtered.obs['total_counts'])),
                'median_genes': int(np.median(adata_filtered.obs['n_genes_by_counts'])),
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
    """从样本名推断分组，取第一个分隔符前的前缀。"""
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
    try:
        cov_inv = np.linalg.pinv(cov)
    except np.linalg.LinAlgError:
        return []
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
