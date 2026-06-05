from modules.base import BaseAnalysis

class QCAnalysis(BaseAnalysis):
    MODULE_NAME = "qc"
    DISPLAY_NAME = "质控"
    DESCRIPTION = "按线粒体比例、基因数、UMI数过滤细胞 + Scrublet 去除双细胞"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import omicverse as ov
        from modules.visualization import umap_scatter, violin_plot
        import plotly.graph_objects as go
        import os, json

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        adata.layers["counts"] = adata.X.copy()

        self.progress(15, "Flagging MT/ribo/hb genes...")
        # 优先用 gene_name 检测（Ensembl ID 不以 MT-/RPS/RPL 开头）
        if 'gene_name' in adata.var.columns:
            gene_names_str = adata.var['gene_name'].fillna('').astype(str)
        else:
            gene_names_str = adata.var_names.astype(str)
        adata.var["mt"] = gene_names_str.str.startswith("MT-")
        adata.var["ribo"] = gene_names_str.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = gene_names_str.str.contains("^HB[^(P)]")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"],
            inplace=True, percent_top=[20], log1p=True
        )

        self.progress(30, "Running Scrublet doublet detection + QC filtering...")
        n_before = adata.shape[0]
        adata = ov.pp.qc(
            adata,
            tresh={
                'mito_perc': float(self.params.get('mito_perc', 0.2)),
                'nUMIs': int(self.params.get('nUMIs', 500)),
                'detected_genes': int(self.params.get('detected_genes', 250)),
            },
            doublets_method='scrublet',
            batch_key=self.params.get('batch_key', 'batch'),
            filter_doublets=True,
        )
        n_after = adata.shape[0]

        self.progress(60, "Generating QC plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        violin_keys = [k for k in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'] if k in adata.obs.columns]
        if violin_keys:
            fig_json = json.dumps(violin_plot(adata, violin_keys, title='QC Metrics'))
            fpath = os.path.join(plots_dir, 'qc_violin.json')
            with open(fpath, 'w') as f: f.write(fig_json)
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'violin', 'label': 'QC Violin Plots'})

        # QC 散点图
        if 'total_counts' in adata.obs.columns and 'n_genes_by_counts' in adata.obs.columns:
            fig_scatter = go.Figure()
            color_vals = adata.obs['pct_counts_mt'].values if 'pct_counts_mt' in adata.obs.columns else None
            fig_scatter.add_trace(go.Scattergl(
                x=adata.obs['total_counts'], y=adata.obs['n_genes_by_counts'],
                mode='markers', marker=dict(size=3, color=color_vals, colorscale='Reds',
                                            colorbar=dict(title='MT%'), opacity=0.6),
                text=adata.obs.index.tolist(),
                hovertemplate='%{text}<br>Counts: %{x:.0f}<br>Genes: %{y:.0f}<br>MT%: %{marker.color:.1f}'
            ))
            fig_scatter.update_layout(title='QC: Counts vs Genes', xaxis_title='Total Counts',
                                     yaxis_title='Detected Genes', plot_bgcolor='white', width=600, height=400)
            fpath = os.path.join(plots_dir, 'qc_scatter.json')
            with open(fpath, 'w') as f: f.write(json.dumps(json.loads(fig_scatter.to_json())))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'scatter', 'label': 'QC Scatter'})

        if 'X_umap' in adata.obsm:
            for color_key in ['batch', 'leiden']:
                if color_key in adata.obs.columns:
                    fig_json = json.dumps(umap_scatter(adata, color_key, title=f'UMAP by {color_key}'))
                    fpath = os.path.join(plots_dir, f'qc_umap_{color_key}.json')
                    with open(fpath, 'w') as f: f.write(fig_json)
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'UMAP by {color_key}'})

        self.progress(85, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'qc_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'cells_before': n_before,
                'cells_after': n_after,
                'cells_removed': n_before - n_after,
                'pct_removed': round((n_before - n_after) / max(n_before, 1) * 100, 1),
                'n_genes': adata.shape[1],
            }
        }
