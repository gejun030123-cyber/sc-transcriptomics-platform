import os
from modules.base import BaseAnalysis

class NormalizeAnalysis(BaseAnalysis):
    MODULE_NAME = "normalize"
    DISPLAY_NAME = "标准化"
    DESCRIPTION = "数据标准化（log1p / Pearson 残差）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import omicverse as ov
        import json
        import numpy as np

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        from modules.io_utils import remap_var_names
        adata = remap_var_names(adata)
        adata.layers["counts"] = adata.X.copy()

        method = self.params.get('method', 'log1p')
        target_sum = float(self.params.get('target_sum', 10000))

        self.progress(20, f"Normalizing ({method})...")
        if method == 'pearson_residuals':
            self.progress(30, "Computing Pearson residuals...")
            clip_values = self.params.get('clip_values', True)
            adata.layers["normalized"] = sc.experimental.pp.normalize_pearson_residuals(
                adata, clip=clip_values, copy=True
            ).X if hasattr(sc.experimental.pp, 'normalize_pearson_residuals') else adata.X
            # Fallback: use omicverse pearson mode
            try:
                adata = ov.pp.preprocess(adata, mode='pearson', target_sum=target_sum)
            except Exception:
                sc.pp.normalize_total(adata, target_sum=target_sum)
                sc.pp.log1p(adata)
        else:
            adata = ov.pp.preprocess(adata, mode='shiftlog', target_sum=target_sum)

        self.progress(70, "Generating normalization plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        # Library size distribution before/after
        import plotly.graph_objects as go
        fig = go.Figure()
        if 'counts' in adata.layers:
            counts_lib = np.array(adata.layers["counts"].sum(axis=1)).flatten()
            fig.add_trace(go.Histogram(x=counts_lib, name='Before', opacity=0.6, nbinsx=50))
        norm_lib = np.array(adata.X.sum(axis=1)).flatten() if hasattr(adata.X, 'sum') else None
        if norm_lib is not None:
            fig.add_trace(go.Histogram(x=norm_lib, name='After', opacity=0.6, nbinsx=50))
        fig.update_layout(title='Library Size Distribution', xaxis_title='Total Counts',
                         yaxis_title='Frequency', barmode='overlay',
                         plot_bgcolor='white', width=600, height=400)
        fpath = os.path.join(plots_dir, 'normalize_libsize.json')
        with open(fpath, 'w') as f: json.dump(json.loads(fig.to_json()), f)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'histogram', 'label': 'Library Size Distribution'})

        self.progress(90, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'normalize_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'method': method,
                'target_sum': target_sum,
            }
        }
