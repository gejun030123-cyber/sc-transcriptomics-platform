from modules.base import BaseAnalysis
from modules.constants import S_GENES, G2M_GENES


class HVGAnalysis(BaseAnalysis):
    MODULE_NAME = "hvg"
    DISPLAY_NAME = "高变异基因选择"
    DESCRIPTION = "选择高变异基因（HVG），支持批次感知和基因过滤"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import scanpy as sc
        import numpy as np
        import json

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)

        n_hvg = int(self.params.get('n_top_genes', 2000))
        batch_key = self.params.get('batch_key', '').strip()
        hvg_flavor = self.params.get('hvg_flavor', 'seurat_v3')
        batch_hvg_strategy = self.params.get('batch_hvg_strategy', 'intersection')
        exclude_mt = self.params.get('exclude_mt_genes', False)
        exclude_cc = self.params.get('exclude_cc_genes', False)
        force_genes_str = self.params.get('force_include_genes', '').strip()
        cc_scoring = self.params.get('cc_scoring', False)
        regress_cc = self.params.get('regress_cc', False)

        self.progress(20, f"Selecting HVGs (flavor={hvg_flavor})...")
        hvg_kwargs = dict(n_top_genes=n_hvg, flavor=hvg_flavor, layer='counts')
        if batch_key and batch_key in adata.obs.columns:
            hvg_kwargs['batch_key'] = batch_key
        sc.pp.highly_variable_genes(adata, **hvg_kwargs)

        # Batch-aware HVG merging (per-batch selection)
        if batch_key and batch_key in adata.obs.columns and 'highly_variable' in adata.var.columns:
            batch_hvgs = set()
            for batch_val in adata.obs[batch_key].unique():
                batch_mask = adata.obs[batch_key] == batch_val
                adata_batch = adata[batch_mask].copy()
                try:
                    sc.pp.highly_variable_genes(adata_batch, n_top_genes=n_hvg, flavor=hvg_flavor, layer='counts')
                    batch_hvgs_batch = set(adata_batch.var_names[adata_batch.var['highly_variable']])
                    if not batch_hvgs:
                        batch_hvgs = batch_hvgs_batch
                    elif batch_hvg_strategy == 'union':
                        batch_hvgs = batch_hvgs | batch_hvgs_batch
                    else:  # intersection
                        batch_hvgs = batch_hvgs & batch_hvgs_batch
                except Exception:
                    pass

            if batch_hvgs:
                adata.var['highly_variable'] = adata.var_names.isin(batch_hvgs)

        # Gene exclusion filters
        if exclude_mt:
            mt_genes = adata.var_names[adata.var_names.str.upper().str.startswith('MT-')]
            adata.var.loc[mt_genes, 'highly_variable'] = False

        if exclude_cc:
            var_names_set = set(adata.var_names.astype(str))
            cc_genes = [g for g in S_GENES + G2M_GENES if g in var_names_set]
            if cc_genes:
                adata.var.loc[cc_genes, 'highly_variable'] = False

        # Force include genes
        force_genes = []
        if force_genes_str:
            force_genes = [g.strip() for g in force_genes_str.replace('\n', ',').split(',') if g.strip()]
            force_found = [g for g in force_genes if g in adata.var_names]
            if force_found:
                adata.var.loc[force_found, 'highly_variable'] = True

        # Cell cycle scoring
        if cc_scoring:
            self.progress(50, "Scoring cell cycle...")
            var_names_set = set(adata.var_names.astype(str))
            s_in = [g for g in S_GENES if g in var_names_set]
            g2m_in = [g for g in G2M_GENES if g in var_names_set]
            if len(s_in) >= 5 and len(g2m_in) >= 5:
                adata_cc = adata.copy()
                if 'log1p' not in adata_cc.uns:
                    sc.pp.normalize_total(adata_cc, target_sum=1e4)
                    sc.pp.log1p(adata_cc)
                sc.tl.score_genes_cell_cycle(adata_cc, s_genes=s_in, g2m_genes=g2m_in)
                adata.obs['S_score'] = adata_cc.obs['S_score']
                adata.obs['G2M_score'] = adata_cc.obs['G2M_score']
                adata.obs['phase'] = adata_cc.obs['phase']
                del adata_cc

        # Regress out cell cycle
        if regress_cc and 'S_score' in adata.obs.columns and 'G2M_score' in adata.obs.columns:
            self.progress(65, "Regressing out cell cycle...")
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'])

        # Subset to HVGs
        adata_hvg = adata[:, adata.var['highly_variable']].copy()

        self.progress(75, "Generating HVG plot...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        import plotly.graph_objects as go
        var_df = adata.var.copy()
        plot_cols = [c for c in ['variances_norm', 'variances', 'dispersions_norm', 'dispersions'] if c in var_df.columns]
        y_col = plot_cols[0] if plot_cols else None
        if y_col:
            var_df = var_df.sort_values(y_col, ascending=False).head(3000)
            fig = go.Figure()
            fig.add_trace(go.Scattergl(
                x=var_df['means'] if 'means' in var_df.columns else range(len(var_df)),
                y=var_df[y_col],
                mode='markers',
                marker=dict(size=2, color=var_df['highly_variable'].map({True: '#e53935', False: '#9e9e9e'}))
            ))
            fig.update_layout(title='Highly Variable Genes', xaxis_title='Mean', yaxis_title=y_col,
                             plot_bgcolor='white', width=600, height=400)
            result_files.append(self.save_plotly_json(fig, plots_dir, 'hvg_scatter.json', 'scatter', 'Highly Variable Genes'))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'hvg')

        n_hvg_actual = int(adata.var['highly_variable'].sum()) if 'highly_variable' in adata.var.columns else n_hvg
        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'n_genes_total': adata.n_vars,
                'n_hvgs': n_hvg_actual,
                'hvg_flavor': hvg_flavor,
                'force_include_count': len(force_genes),
            }
        }
