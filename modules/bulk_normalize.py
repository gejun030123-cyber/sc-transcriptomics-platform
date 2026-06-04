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
        if input_path.endswith(('.xlsx', '.xls')):
            df = pd.read_excel(input_path, index_col=0)
            adata = sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
        elif input_path.endswith('.csv') or input_path.endswith('.txt'):
            df = pd.read_csv(input_path, sep=None if input_path.endswith('.csv') else '\t', index_col=0)
            adata = sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
        else:
            adata = sc.read_h5ad(input_path)

        method = self.params.get('method', 'deseq2')
        self.progress(20, f"标准化方法: {method}...")

        if method == 'deseq2':
            from scipy.stats import gmean
            counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
            counts = counts.astype(float)
            counts = counts[counts.sum(axis=1) > 0]
            geo_means = gmean(counts + 1, axis=0)
            ratios = counts / (geo_means + 1e-10)
            size_factors = np.median(ratios, axis=1)
            adata.obs['size_factor'] = size_factors
            norm_counts = counts / size_factors[:, None]
            adata.layers['normalized'] = norm_counts
            adata.X = np.log2(norm_counts + 1)

        elif method == 'cpm':
            counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
            lib_sizes = counts.sum(axis=1, keepdims=True)
            cpm = counts / lib_sizes * 1e6
            adata.layers['normalized'] = cpm
            adata.X = np.log2(cpm + 1)

        elif method == 'log2_quantile':
            from scipy.stats import gmean
            counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
            log_counts = np.log2(counts + 1)
            from scipy.stats import rankdata
            ranked = np.apply_along_axis(rankdata, 0, log_counts)
            ref_distribution = np.sort(np.mean(log_counts, axis=1))
            norm = np.zeros_like(log_counts)
            for i in range(log_counts.shape[1]):
                sorted_idx = np.argsort(ranked[:, i])
                norm[sorted_idx, i] = ref_distribution
            adata.layers['normalized'] = 2**norm - 1
            adata.X = norm

        self.progress(60, "生成标准化前后对比图...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        raw_lib = counts.sum(axis=1) if 'counts' in dir() else np.ones(adata.n_obs)
        norm_layer = adata.layers.get('normalized', adata.X)
        norm_lib = norm_layer.sum(axis=1) if hasattr(norm_layer, 'sum') else np.ones(adata.n_obs)

        from plotly.subplots import make_subplots
        fig = make_subplots(rows=1, cols=2, subplot_titles=['标准化前 (Raw)', '标准化后 (Normalized)'])
        fig.add_trace(go.Bar(y=raw_lib, marker_color='#e53935', name='Raw'), row=1, col=1)
        fig.add_trace(go.Bar(y=norm_lib, marker_color='#4caf50', name='Normalized'), row=1, col=2)
        fig.update_layout(height=350, width=700, showlegend=False, title='文库大小对比')
        fpath = os.path.join(plots_dir, 'bulk_norm_libsize.json')
        with open(fpath, 'w') as f: json.dump(json.loads(fig.to_json()), f)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'bar', 'label': '文库大小对比'})

        if method == 'deseq2':
            fig_sf = go.Figure()
            fig_sf.add_trace(go.Bar(x=adata.obs.index.tolist(), y=adata.obs['size_factor'].values,
                                   marker_color='#1a237e'))
            fig_sf.update_layout(title='DESeq2 Size Factors', yaxis_title='Size Factor',
                                plot_bgcolor='white', width=600, height=300)
            fpath = os.path.join(plots_dir, 'bulk_norm_sizefactors.json')
            with open(fpath, 'w') as f: json.dump(json.loads(fig_sf.to_json()), f)
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
