import os
from modules.base import BaseAnalysis

class PreprocessAnalysis(BaseAnalysis):
    MODULE_NAME = "preprocess"
    DISPLAY_NAME = "预处理"
    DESCRIPTION = "标准化，选择高变异基因"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import omicverse as ov
        from modules.visualization import umap_scatter
        import json

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)

        n_hvg = int(self.params.get('n_top_genes', 2000))
        target_sum = float(self.params.get('target_sum', 10000))

        self.progress(20, "Normalizing and selecting HVGs (shiftlog|pearson)...")
        adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', target_sum=target_sum)

        self.progress(50, f"Selecting top {n_hvg} HVGs...")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor='seurat_v3', layer='counts')
        adata_hvg = adata[:, adata.var['highly_variable']].copy()

        self.progress(70, "Generating HVG plot...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        import plotly.graph_objects as go
        var_df = adata.var.copy()
        var_df = var_df.sort_values('variances_norm', ascending=False).head(3000)
        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            x=var_df['means'], y=var_df['variances_norm'],
            mode='markers',
            marker=dict(size=2, color=var_df['highly_variable'].map({True: '#e53935', False: '#9e9e9e'}))
        ))
        fig.update_layout(title='Highly Variable Genes', xaxis_title='Mean', yaxis_title='Normalized Variance',
                         plot_bgcolor='white', width=600, height=400)
        fpath = os.path.join(plots_dir, 'preprocess_hvg.json')
        with open(fpath, 'w') as f: json.dump(json.loads(fig.to_json()), f)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'scatter', 'label': 'Highly Variable Genes'})

        self.progress(85, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'preprocess_output.h5ad')
        adata.write_h5ad(output_path)

        n_hvg_actual = int(adata.var['highly_variable'].sum()) if 'highly_variable' in adata.var.columns else n_hvg
        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'n_genes_total': adata.n_vars,
                'n_hvgs': n_hvg_actual,
                'target_sum': target_sum,
            }
        }
