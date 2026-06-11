import os
from modules.base import BaseAnalysis

class DimredAnalysis(BaseAnalysis):
    MODULE_NAME = "dimred"
    DISPLAY_NAME = "降维分析"
    DESCRIPTION = "标准化、PCA、UMAP 降维"
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
        from modules.io_utils import remap_var_names
        adata = remap_var_names(adata)
        n_comps = int(self.params.get('n_comps', 50))
        use_mde = self.params.get('use_mde', False)

        self.progress(20, "Scaling data...")
        ov.pp.scale(adata, max_value=10)

        self.progress(35, f"Running PCA ({n_comps} components)...")
        sc.pp.pca(adata, n_comps=n_comps, layer='scaled')

        self.progress(55, "Computing neighbors...")
        sc.pp.neighbors(adata, n_pcs=n_comps)

        self.progress(70, "Computing 2D embedding...")
        embedding_method = 'umap'
        if use_mde:
            try:
                import pymde
                mde = pymde.preserve_neighbors(adata.obsm['X_pca'], embedding_dim=2, device='cpu')
                embedding = mde.embed(verbose=False)
                adata.obsm['X_mde'] = embedding.cpu().numpy() if hasattr(embedding, 'cpu') else embedding.numpy()
                adata.obsm['X_umap'] = adata.obsm['X_mde']
                embedding_method = 'mde'
            except ImportError:
                self.progress(71, "pymde not installed, falling back to UMAP...")
                sc.tl.umap(adata)
        else:
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
                'embedding_method': embedding_method,
            }
        }
