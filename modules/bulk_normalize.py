import os
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis

class BulkNormalizeAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_normalize"
    DISPLAY_NAME = "Bulk 数据标准化"
    DESCRIPTION = "计数矩阵标准化：DESeq2 size factors、CPM、分位数标准化"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import plotly.graph_objects as go

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)

        method = self.params.get('method', 'deseq2')
        self.progress(20, f"标准化方法: {method}...")

        # 前置过滤参数
        min_expr_value = float(self.params.get('min_expr_value', 1))
        min_expr_samples = int(self.params.get('min_expr_samples', 3))
        max_zero_pct = float(self.params.get('max_zero_pct', 0))

        # 低表达基因前置过滤
        n_genes_before = adata.n_vars
        if min_expr_samples > 0:
            from scipy import sparse as _sp
            raw_for_filter = adata.X.toarray() if _sp.issparse(adata.X) else np.asarray(adata.X)
            lib_for_filter = raw_for_filter.sum(axis=1, keepdims=True)
            lib_for_filter[lib_for_filter == 0] = 1
            cpm_check = raw_for_filter / lib_for_filter * 1e6
            n_expr = np.array((cpm_check >= min_expr_value).sum(axis=0)).flatten()
            gene_mask = n_expr >= min_expr_samples
            if max_zero_pct > 0:
                zero_pct = np.array((raw_for_filter == 0).sum(axis=0)).flatten() / adata.n_obs * 100
                gene_mask = gene_mask & (zero_pct <= max_zero_pct)
            adata = adata[:, gene_mask].copy()

        if adata.n_vars == 0:
            raise ValueError("所有基因均被过滤掉（低表达），请降低 min_expr_samples 参数")

        from scipy import sparse as _sp
        raw_counts = np.asarray(adata.X.toarray()) if _sp.issparse(adata.X) else np.asarray(adata.X)
        raw_counts = raw_counts.astype(float)
        adata.layers['raw'] = raw_counts.copy()
        raw_lib = raw_counts.sum(axis=1)

        if method == 'deseq2':
            size_factors = _estimate_size_factors(raw_counts)
            adata.obs['size_factor'] = size_factors
            norm_counts = raw_counts / size_factors[:, None]
            adata.layers['normalized'] = norm_counts
            adata.X = np.log2(norm_counts + 1)
            adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            adata.layers['normalized'] = np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)

        elif method == 'tmm':
            tmm_factors = _tmm_normalize(raw_counts)
            tmm_factors = np.where(tmm_factors > 0, tmm_factors, 1.0)
            adata.obs['size_factor'] = tmm_factors
            norm_counts = raw_counts / tmm_factors[:, None]
            adata.layers['normalized'] = norm_counts
            adata.X = np.log2(norm_counts + 1)
            adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            adata.layers['normalized'] = np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)

        elif method == 'cpm':
            cpm_target = float(self.params.get('cpm_target', 1e6))
            lib_sizes = raw_counts.sum(axis=1, keepdims=True)
            lib_sizes[lib_sizes == 0] = 1  # 避免除零
            cpm = raw_counts / lib_sizes * cpm_target
            adata.layers['normalized'] = cpm
            adata.X = np.log2(cpm + 1)

            # 清理 inf/NaN 值
            adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            adata.layers['normalized'] = np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)

        elif method == 'log2_quantile':
            log_counts = np.log2(raw_counts + 1)
            from scipy.stats import rankdata
            ranked = np.apply_along_axis(rankdata, 0, log_counts)
            ref_distribution = np.mean(np.sort(log_counts, axis=0), axis=1)
            norm = np.zeros_like(log_counts)
            for i in range(log_counts.shape[1]):
                sorted_idx = np.argsort(ranked[:, i])
                norm[sorted_idx, i] = ref_distribution
            adata.layers['normalized'] = 2**norm - 1
            adata.X = norm

            # 清理 inf/NaN 值
            adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            adata.layers['normalized'] = np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)

        elif method == 'vst':
            size_factors = _estimate_size_factors(raw_counts)
            adata.obs['size_factor'] = size_factors
            adata.X = _vst_transform(raw_counts, size_factors)

        elif method == 'rlog':
            size_factors = _estimate_size_factors(raw_counts)
            adata.obs['size_factor'] = size_factors
            adata.X = _rlog_transform(raw_counts, size_factors)

        else:
            raise ValueError(f"未知标准化方法: {method}，支持: deseq2/tmm/cpm/log2_quantile/vst/rlog")

        # 统一输出标记
        adata.uns['normalization'] = {
            'method': method,
            'is_log_transformed': method in ('deseq2', 'tmm', 'cpm', 'vst', 'log2_quantile', 'rlog'),
        }

        self.progress(60, "生成标准化前后对比图...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        # 文库大小对比图：仅对线性尺度方法有意义
        if method in ('vst', 'rlog'):
            norm_lib = raw_lib  # vst/rlog 无 normalized 层，用原始文库大小展示
        else:
            norm_layer = adata.layers.get('normalized', adata.X)
            norm_lib = norm_layer.sum(axis=1) if hasattr(norm_layer, 'sum') else np.ones(adata.n_obs)

        from plotly.subplots import make_subplots
        fig = make_subplots(rows=1, cols=2, subplot_titles=['标准化前 (Raw)', '标准化后 (Normalized)'])
        fig.add_trace(go.Bar(y=raw_lib, marker_color='#e53935', name='Raw'), row=1, col=1)
        fig.add_trace(go.Bar(y=norm_lib, marker_color='#4caf50', name='Normalized'), row=1, col=2)
        fig.update_layout(height=350, width=700, showlegend=False, title='文库大小对比')
        fpath = os.path.join(plots_dir, 'bulk_norm_libsize.json')
        with open(fpath, 'w') as f: f.write(fig.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'bar', 'label': '文库大小对比'})

        if method in ('deseq2', 'tmm') and 'size_factor' in adata.obs.columns:
            fig_sf = go.Figure()
            fig_sf.add_trace(go.Bar(x=adata.obs.index.tolist(), y=adata.obs['size_factor'].values,
                                   marker_color='#1a237e'))
            fig_sf.update_layout(title=f'{method.upper()} Size Factors', yaxis_title='Size Factor',
                                plot_bgcolor='white', width=600, height=300)
            fpath = os.path.join(plots_dir, 'bulk_norm_sizefactors.json')
            with open(fpath, 'w') as f: f.write(fig_sf.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'bar', 'label': 'Size Factors'})

        # 标准化前后表达分布对比
        sample_labels = adata.obs.index.tolist()
        fig_box = make_subplots(rows=1, cols=2, subplot_titles=['标准化前 (log2 raw)', '标准化后'])
        # 使用单个 violin trace 展示所有样本的分布
        raw_log2 = np.log2(raw_counts + 1)
        fig_box.add_trace(go.Violin(y=raw_log2.flatten(), name='Raw', box_visible=True,
                                     meanline_visible=True, marker_color='#e53935', showlegend=False), row=1, col=1)
        norm_flat = adata.X.flatten()
        fig_box.add_trace(go.Violin(y=norm_flat, name='Normalized', box_visible=True,
                                     meanline_visible=True, marker_color='#4caf50', showlegend=False), row=1, col=2)
        fig_box.update_layout(height=400, width=800, title='标准化前后表达分布对比')
        fpath = os.path.join(plots_dir, 'bulk_norm_boxplot_compare.json')
        with open(fpath, 'w') as f: f.write(fig_box.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '表达分布对比'})

        # PCA 前后对比
        n_pcs = min(10, adata.n_obs - 1)
        if n_pcs >= 2:
            fig_pca = make_subplots(rows=1, cols=2, subplot_titles=['原始数据 PCA', '标准化后 PCA'])
            # 原始数据 PCA
            adata_raw_pca = sc.AnnData(X=np.log2(raw_counts + 1), obs=adata.obs.copy())
            sc.pp.scale(adata_raw_pca, max_value=10)
            sc.pp.pca(adata_raw_pca, n_comps=n_pcs)
            pc_raw = adata_raw_pca.obsm['X_pca']
            fig_pca.add_trace(go.Scattergl(x=pc_raw[:, 0].tolist(), y=pc_raw[:, 1].tolist(),
                mode='markers+text', text=sample_labels, textposition='top center',
                marker=dict(size=8, color='#e53935'), showlegend=False), row=1, col=1)
            # 标准化后 PCA
            adata_norm_pca = sc.AnnData(X=adata.X.copy(), obs=adata.obs.copy())
            sc.pp.scale(adata_norm_pca, max_value=10)
            sc.pp.pca(adata_norm_pca, n_comps=n_pcs)
            pc_norm = adata_norm_pca.obsm['X_pca']
            fig_pca.add_trace(go.Scattergl(x=pc_norm[:, 0].tolist(), y=pc_norm[:, 1].tolist(),
                mode='markers+text', text=sample_labels, textposition='top center',
                marker=dict(size=8, color='#4caf50'), showlegend=False), row=1, col=2)
            fig_pca.update_layout(height=400, width=900, title='标准化前后 PCA 对比',
                                  plot_bgcolor='white')
            fpath = os.path.join(plots_dir, 'bulk_norm_pca_compare.json')
            with open(fpath, 'w') as f: f.write(fig_pca.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': 'PCA 前后对比'})

        self.progress(85, "保存结果...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_normalize_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'method': method,
                'n_samples': adata.n_obs,
                'n_genes_before_filter': n_genes_before,
                'n_genes_after_filter': adata.n_vars,
                'genes_filtered': n_genes_before - adata.n_vars,
                'median_size_factor': round(float(adata.obs['size_factor'].median()), 3) if 'size_factor' in adata.obs.columns else None,
                'is_log_transformed': adata.uns.get('normalization', {}).get('is_log_transformed', True),
            }
        }


# --- 辅助函数 ---


def _tmm_normalize(counts):
    """TMM 标准化（edgeR 风格）。counts: (n_samples, n_genes) ndarray of raw counts.
    返回 per-sample TMM 因子。"""
    n_samples, n_genes = counts.shape
    lib_sizes = counts.sum(axis=1)
    lib_sizes[lib_sizes == 0] = 1

    # 参考样本：文库大小最接近中位数
    median_lib = np.median(lib_sizes)
    ref_idx = int(np.argmin(np.abs(lib_sizes - median_lib)))

    factors = np.ones(n_samples)
    for i in range(n_samples):
        if i == ref_idx:
            continue
        mask = (counts[i] > 0) & (counts[ref_idx] > 0)
        if mask.sum() < 10:
            continue
        fi = lib_sizes[i]
        fr = lib_sizes[ref_idx]
        Mi = np.log2((counts[i, mask] / fi) / (counts[ref_idx, mask] / fr))
        Ai = 0.5 * (np.log2(counts[i, mask] / fi) + np.log2(counts[ref_idx, mask] / fr))

        m_lo, m_hi = np.percentile(Mi, [30, 70])
        a_hi = np.percentile(Ai, 95)
        keep = (Mi >= m_lo) & (Mi <= m_hi) & (Ai <= a_hi)

        if keep.sum() > 0:
            wi = (fi - counts[i, mask][keep]) / (fi * counts[i, mask][keep])
            wr = (fr - counts[ref_idx, mask][keep]) / (fr * counts[ref_idx, mask][keep])
            weights = 1.0 / (wi + wr + 1e-30)
            factors[i] = 2 ** (-np.average(Mi[keep], weights=weights))

    # 归一化使几何均值为 1
    log_factors = np.log(factors[factors > 0])
    if len(log_factors) > 0:
        geo_mean = np.exp(np.mean(log_factors))
        factors = factors / geo_mean

    return factors


def _vst_transform(counts, size_factors):
    """近似 VST。counts: raw counts (n_samples, n_genes), size_factors: (n_samples,)."""
    normed = counts / size_factors[:, None]
    normed = np.maximum(normed, 0)
    vst = np.log2(normed + 0.5)
    return np.nan_to_num(vst, nan=0.0, posinf=0.0, neginf=0.0)


def _rlog_transform(counts, size_factors, prior_mean=None):
    """近似 rlog。小样本时通过正则化收缩基因效应。"""
    normed = counts / size_factors[:, None]
    normed = np.maximum(normed, 0)
    log_normed = np.log2(normed + 0.5)

    n_samples = counts.shape[0]
    if prior_mean is None:
        prior_mean = float(np.mean(log_normed))

    gene_effects = log_normed.mean(axis=0) - prior_mean
    shrinkage = min(1.0, n_samples / 30.0)

    sample_residuals = log_normed - log_normed.mean(axis=0, keepdims=True)
    rlog = prior_mean + shrinkage * gene_effects[None, :] + sample_residuals

    return np.nan_to_num(rlog, nan=0.0, posinf=0.0, neginf=0.0)


def _estimate_size_factors(raw_counts):
    """DESeq2 中位比率法计算 size factors。raw_counts: (n_samples, n_genes)."""
    from scipy.stats import gmean
    # 用所有样本（不过滤全零样本）
    nonzero_all = (raw_counts > 0).all(axis=0)
    if nonzero_all.sum() == 0:
        geo_means = np.exp(np.log(raw_counts + 1).mean(axis=0))
    else:
        geo_means = np.ones(raw_counts.shape[1])
        geo_means[nonzero_all] = gmean(raw_counts[:, nonzero_all], axis=0)
    ratios = raw_counts / (geo_means + 1e-10)
    size_factors = np.median(ratios, axis=1)
    size_factors = np.where(size_factors > 0, size_factors, 1.0)
    return size_factors
