import os
import json
from modules.base import BaseAnalysis

class DEGAnalysis(BaseAnalysis):
    MODULE_NAME = "deg"
    DISPLAY_NAME = "差异表达"
    DESCRIPTION = "差异表达基因分析（Wilcoxon 检验）"
    INPUT_REQUIRES = ['leiden']

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import pandas as pd
        from modules.visualization import umap_scatter
        import plotly.graph_objects as go

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        groupby = self.params.get('groupby', 'celltype')
        method = self.params.get('method', 'wilcoxon')
        n_genes = int(self.params.get('n_genes', 20))

        if groupby not in adata.obs.columns:
            groupby = 'leiden'

        self.progress(20, f"Running DEG analysis ({method})...")
        sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=100)

        self.progress(50, "Extracting results...")
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        deg_data = []
        for g in groups:
            for i in range(min(n_genes, len(result['names'][g]))):
                deg_data.append({
                    'cluster': g,
                    'gene': result['names'][g][i],
                    'logfc': round(float(result['logfoldchanges'][g][i]), 3),
                    'pval': float(result['pvals'][g][i]),
                    'pval_adj': float(result['pvals_adj'][g][i]),
                    'score': round(float(result['scores'][g][i]), 3),
                })

        self.progress(70, "Generating volcano plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        deg_df = pd.DataFrame(deg_data)
        csv_path = os.path.join(results_dir, 'deg_results.csv')
        deg_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': 'DEG Results Table'})

        first_group = groups[0] if groups else None
        if first_group:
            g_mask = deg_df['cluster'] == first_group
            g_df = deg_df[g_mask].copy()
            g_df['-log10(pval_adj)'] = -g_df['pval_adj'].apply(lambda x: __import__('math').log10(max(x, 1e-300)))
            fig = go.Figure()
            sig = (g_df['pval_adj'] < 0.05) & (g_df['logfc'].abs() > 1)
            fig.add_trace(go.Scattergl(x=g_df.loc[sig, 'logfc'], y=g_df.loc[sig, '-log10(pval_adj)'],
                                       mode='markers', marker=dict(color='#e53935', size=4), name='Significant'))
            fig.add_trace(go.Scattergl(x=g_df.loc[~sig, 'logfc'], y=g_df.loc[~sig, '-log10(pval_adj)'],
                                       mode='markers', marker=dict(color='#9e9e9e', size=3), name='Not significant'))
            fig.update_layout(title=f'Volcano Plot: {first_group}', xaxis_title='Log2 FC', yaxis_title='-log10(padj)',
                             plot_bgcolor='white', width=600, height=400)
            fpath = os.path.join(plots_dir, f'deg_volcano_{first_group}.json')
            with open(fpath, 'w') as f: json.dump(json.loads(fig.to_json()), f)
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'volcano', 'label': f'Volcano: {first_group}'})

        self.progress(90, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'deg_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_groups': len(groups),
                'groups': list(groups),
                'method': method,
                'total_deg_genes': len(deg_data),
            }
        }
