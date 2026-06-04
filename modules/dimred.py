import os
from modules.base import BaseAnalysis

class DimredAnalysis(BaseAnalysis):
    MODULE_NAME = "dimred"
    DISPLAY_NAME = "Dimensionality Reduction"
    DESCRIPTION = "Scale, PCA, UMAP"
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
        n_comps = int(self.params.get('n_comps', 50))

        self.progress(20, "Scaling data...")
        ov.pp.scale(adata, max_value=10)

        self.progress(35, f"Running PCA ({n_comps} components)...")
        ov.pp.pca(adata, n_comps=n_comps)

        self.progress(55, "Computing neighbors...")
        sc.pp.neighbors(adata, n_pcs=n_comps)

        self.progress(70, "Computing UMAP...")
        sc.tl.umap(adata)

        self.progress(80, "Generating UMAP plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        for color_key in ['batch', 'leiden', 'n_genes_by_counts']:
            if color_key in adata.obs.columns:
                fig_json = json.dumps(umap_scatter(adata, color_key, title=f'UMAP colored by {color_key}'))
                fpath = os.path.join(plots_dir, f'dimred_umap_{color_key}.json')
                with open(fpath, 'w') as f: f.write(fig_json)
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'UMAP by {color_key}'})

        self.progress(90, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'dimred_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'n_pcs': n_comps,
                'pca_variance_ratio_top5': round(float(adata.uns['pca']['variance_ratio'][:5].sum()), 3) if 'pca' in adata.uns else None,
            }
        }
