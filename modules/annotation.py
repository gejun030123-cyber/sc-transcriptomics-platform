import os
from modules.base import BaseAnalysis

DEFAULT_TME_MARKERS = {
    'Tumor Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
    'CAF': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'FAP', 'ACTA2'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'ENG', 'CLDN5'],
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'TRAC'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1', 'NCAM1', 'PRF1'],
    'B cells': ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'PAX5'],
    'Plasma cells': ['JCHAIN', 'MZB1', 'SDC1', 'IGHG1', 'IGKC'],
    'Monocyte/Macrophage': ['CD14', 'CD68', 'CSF1R', 'LYZ', 'S100A8', 'S100A9'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA'],
    'Neutrophils': ['CSF3R', 'CXCR2', 'FCGR3B', 'S100A12'],
    'Mast cells': ['KIT', 'TPSAB1', 'TPSB2', 'HDC', 'MS4A2'],
    'Pericytes': ['RGS5', 'PDGFRB', 'NOTCH3', 'MCAM', 'ACTA2'],
    'Proliferating': ['MKI67', 'TOP2A', 'PCNA', 'STMN1', 'CDK1'],
}

class AnnotationAnalysis(BaseAnalysis):
    MODULE_NAME = "annotation"
    DISPLAY_NAME = "细胞注释"
    DESCRIPTION = "基于 Marker 基因的细胞类型自动注释"
    INPUT_REQUIRES = ['leiden']

    def validate_input(self, adata):
        cluster_key = self.params.get('cluster_key', 'leiden')
        if cluster_key not in adata.obs.columns:
            return f"Column '{cluster_key}' not found in adata.obs."
        return None

    def run(self, input_path):
        import scanpy as sc
        from modules.visualization import umap_scatter
        import json

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        cluster_key = self.params.get('cluster_key', 'leiden')
        resolution = self.params.get('resolution', '0.8')
        leiden_key = f'leiden_{resolution}' if f'leiden_{resolution}' in adata.obs.columns else cluster_key

        self.progress(20, "Scoring cell type markers...")
        markers = DEFAULT_TME_MARKERS
        for ct, genes in markers.items():
            available_genes = [g for g in genes if g in adata.var_names]
            if available_genes:
                sc.tl.score_genes(adata, available_genes, score_name=f'score_{ct}', use_raw=False)

        self.progress(50, "Assigning cell types to clusters...")
        score_cols = [f'score_{ct}' for ct in markers if f'score_{ct}' in adata.obs.columns]
        if score_cols:
            cluster_annotations = {}
            for cluster in adata.obs[leiden_key].cat.categories:
                mask = adata.obs[leiden_key] == cluster
                mean_scores = {col: adata.obs.loc[mask, col].mean() for col in score_cols}
                best_col = max(mean_scores, key=mean_scores.get)
                best_ct = best_col.replace('score_', '')
                cluster_annotations[cluster] = best_ct
            adata.obs['celltype'] = adata.obs[leiden_key].map(cluster_annotations).astype('category')
        else:
            adata.obs['celltype'] = adata.obs[leiden_key].astype(str)

        self.progress(70, "Generating dotplot and UMAP...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        import plotly.graph_objects as go
        dotplot_genes = []
        for genes in markers.values():
            dotplot_genes.extend([g for g in genes[:2] if g in adata.var_names])
        dotplot_genes = list(dict.fromkeys(dotplot_genes))[:20]

        if dotplot_genes and 'celltype' in adata.obs.columns:
            sc.tl.dendrogram(adata, groupby='celltype')
            fig_dotplot = sc.pl.dotplot(adata, var_names=dotplot_genes, groupby='celltype', return_fig=True)
            import matplotlib.pyplot as plt
            import io, base64
            buf = io.BytesIO()
            fig_dotplot.savefig(buf, format='png', dpi=100, bbox_inches='tight')
            plt.close('all')
            buf.seek(0)
            img_b64 = base64.b64encode(buf.read()).decode()
            fpath = os.path.join(plots_dir, 'annotation_dotplot.json')
            with open(fpath, 'w') as f:
                json.dump({'data': [{'type': 'image', 'source': f'data:image/png;base64,{img_b64}', 'xref': 'paper', 'yref': 'paper', 'x': 0, 'y': 1, 'sizex': 1, 'sizey': 1, 'sizing': 'stretch'}], 'layout': {'width': 800, 'height': 500, 'title': 'Cell Type Marker Dotplot'}}, f)
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'dotplot', 'label': 'Cell Type Dotplot'})

        fig_json = json.dumps(umap_scatter(adata, 'celltype', title='UMAP by Cell Type'))
        fpath = os.path.join(plots_dir, 'annotation_umap_celltype.json')
        with open(fpath, 'w') as f: f.write(fig_json)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': 'UMAP by Cell Type'})

        self.progress(90, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'annotation_output.h5ad')
        adata.write_h5ad(output_path)

        ct_counts = adata.obs['celltype'].value_counts().to_dict()
        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_celltypes': adata.obs['celltype'].nunique(),
                'celltype_counts': {str(k): int(v) for k, v in ct_counts.items()},
                'cluster_column': leiden_key,
            }
        }
