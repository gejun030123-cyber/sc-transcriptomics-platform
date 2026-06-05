import os
import json
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
        from modules.visualization import scatter_plot

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)

        method = self.params.get('method', 'deseq2')
        self.progress(20, f"标准化方法: {method}...")

        raw_counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        raw_counts = raw_counts.astype(float)
        raw_lib = raw_counts.sum(axis=1)

        if method == 'deseq2':
            from scipy.stats import gmean
            counts = raw_counts[raw_counts.sum(axis=1) > 0]
            # 只用所有样本都 >0 的基因计算 geometric mean（DESeq2 标准做法）
            nonzero_mask = (counts > 0).all(axis=0)
            if nonzero_mask.sum() == 0:
                # Fallback: 没有全非零基因，用 log-based gmean
                geo_means = np.exp(np.log(counts + 1).mean(axis=0))
            else:
                geo_means = np.ones(counts.shape[1])
                geo_means[nonzero_mask] = gmean(counts[:, nonzero_mask], axis=0)
            ratios = counts / (geo_means + 1e-10)
            size_factors = np.median(ratios, axis=1)
            # 防止 size_factor 为 0
            size_factors = np.where(size_factors > 0, size_factors, 1.0)
            adata.obs['size_factor'] = size_factors
            norm_counts = counts / size_factors[:, None]
            adata.layers['normalized'] = norm_counts
            adata.X = np.log2(norm_counts + 1)

            # 清理 inf/NaN 值
            import numpy as _np
            adata.X = _np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            if 'normalized' in adata.layers:
                adata.layers['normalized'] = _np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)

        elif method == 'cpm':
            lib_sizes = raw_counts.sum(axis=1, keepdims=True)
            cpm = raw_counts / lib_sizes * 1e6
            adata.layers['normalized'] = cpm
            adata.X = np.log2(cpm + 1)

            # 清理 inf/NaN 值
            import numpy as _np
            adata.X = _np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            if 'normalized' in adata.layers:
                adata.layers['normalized'] = _np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)

        elif method == 'log2_quantile':
            from scipy.stats import gmean
            log_counts = np.log2(raw_counts + 1)
            from scipy.stats import rankdata
            ranked = np.apply_along_axis(rankdata, 0, log_counts)
            ref_distribution = np.sort(np.mean(log_counts, axis=1))
            norm = np.zeros_like(log_counts)
            for i in range(log_counts.shape[1]):
                sorted_idx = np.argsort(ranked[:, i])
                norm[sorted_idx, i] = ref_distribution
            adata.layers['normalized'] = 2**norm - 1
            adata.X = norm

            # 清理 inf/NaN 值
            import numpy as _np
            adata.X = _np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            if 'normalized' in adata.layers:
                adata.layers['normalized'] = _np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)

        self.progress(60, "生成标准化前后对比图...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

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

        if method == 'deseq2':
            fig_sf = go.Figure()
            fig_sf.add_trace(go.Bar(x=adata.obs.index.tolist(), y=adata.obs['size_factor'].values,
                                   marker_color='#1a237e'))
            fig_sf.update_layout(title='DESeq2 Size Factors', yaxis_title='Size Factor',
                                plot_bgcolor='white', width=600, height=300)
            fpath = os.path.join(plots_dir, 'bulk_norm_sizefactors.json')
            with open(fpath, 'w') as f: f.write(fig_sf.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'bar', 'label': 'Size Factors'})

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
                'n_genes': adata.n_vars,
                'median_size_factor': round(float(adata.obs['size_factor'].median()), 3) if 'size_factor' in adata.obs.columns else None,
            }
        }
