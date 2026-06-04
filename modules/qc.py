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
        import os, json

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)

        self.progress(15, "Flagging MT/ribo/hb genes...")
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
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
