import os
import json
from modules.base import BaseAnalysis

class TrajectoryAnalysis(BaseAnalysis):
    MODULE_NAME = "trajectory"
    DISPLAY_NAME = "轨迹分析"
    DESCRIPTION = "基于扩散图或 Slingshot 的拟时序分析"
    INPUT_REQUIRES = ['X_umap']

    def validate_input(self, adata):
        if 'X_umap' not in adata.obsm:
            return "UMAP not found. Run dimensionality reduction first."
        return None

    def run(self, input_path):
        import scanpy as sc
        from modules.visualization import umap_scatter
        import plotly.graph_objects as go

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        from modules.io_utils import remap_var_names
        adata = remap_var_names(adata)
        method = self.params.get('method', 'diffusion_map')
        cluster_key = self.params.get('cluster_key', 'leiden')

        self.progress(20, f"Computing diffusion map...")
        sc.tl.diffmap(adata)

        self.progress(40, "Computing diffusion pseudotime...")
        sc.tl.dpt(adata, n_branchings=0, n_dcs=10)

        self.progress(60, "Generating trajectory plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        fig_json = json.dumps(umap_scatter(adata, 'dpt_pseudotime', title='Diffusion Pseudotime'))
        fpath = os.path.join(plots_dir, 'trajectory_pseudotime.json')
        with open(fpath, 'w') as f: f.write(fig_json)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': 'Pseudotime UMAP'})

        if 'X_diffmap' in adata.obsm:
            dc = adata.obsm['X_diffmap'][:, :2]
            fig = go.Figure()
            color_vals = adata.obs['dpt_pseudotime'].values if 'dpt_pseudotime' in adata.obs.columns else None
            fig.add_trace(go.Scattergl(x=dc[:, 0], y=dc[:, 1], mode='markers',
                                       marker=dict(size=2, color=color_vals, colorscale='Viridis', colorbar=dict(title='Pseudotime'))))
            fig.update_layout(title='Diffusion Map', xaxis_title='DC1', yaxis_title='DC2',
                             plot_bgcolor='white', width=600, height=500)
            fpath = os.path.join(plots_dir, 'trajectory_diffmap.json')
            with open(fpath, 'w') as f: json.dump(json.loads(fig.to_json()), f)
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'diffusion_map', 'label': 'Diffusion Map'})

        self.progress(85, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'trajectory_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'method': method,
                'max_pseudotime': round(float(adata.obs['dpt_pseudotime'].max()), 3) if 'dpt_pseudotime' in adata.obs.columns else None,
                'n_cells': adata.n_obs,
            }
        }
