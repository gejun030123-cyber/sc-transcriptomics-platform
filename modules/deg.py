from modules.base import BaseAnalysis

class DEGAnalysis(BaseAnalysis):
    MODULE_NAME = "deg"
    DISPLAY_NAME = "差异表达"
    DESCRIPTION = "差异表达基因分析（Wilcoxon 检验）"
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        import scanpy as sc
        import pandas as pd
        from modules.visualization import umap_scatter
        import plotly.graph_objects as go
        import os, json

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        groupby = self.params.get('groupby', 'celltype')
        method = self.params.get('method', 'wilcoxon')
        n_genes = int(self.params.get('n_genes', 20))
        reference = self.params.get('reference', 'rest')
        pval_cutoff = float(self.params.get('pval_cutoff', 0.05))
        logfc_cutoff = float(self.params.get('logfc_cutoff', 1.0))
        min_pct = float(self.params.get('min_pct', 0.1))
        correction_method = self.params.get('correction_method', 'benjamini_hochberg')
        volcano_top_n = int(self.params.get('volcano_top_n', 10))
        volcano_genes_str = self.params.get('volcano_genes', '').strip()

        if groupby not in adata.obs.columns:
            groupby = 'leiden'

        self.progress(20, f"Running DEG analysis ({method})...")
        ref_kwarg = {} if reference == 'rest' else {'reference': str(reference)}
        sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=100, **ref_kwarg)

        self.progress(50, "Extracting results...")
        from statsmodels.stats.multitest import multipletests
        import numpy as np

        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names

        # 预计算每组每基因的表达比例（用于 min_pct 过滤）
        pct_expr = {}
        if min_pct > 0:
            for g in groups:
                mask = adata.obs[groupby] == g
                n_cells = mask.sum()
                if n_cells > 0:
                    expr = (adata[mask].X > 0).sum(axis=0)
                    if hasattr(expr, 'A1'):
                        expr = expr.A1
                    pct_expr[g] = dict(zip(adata.var_names, expr / n_cells))

        deg_data = []
        correction_map = {'benjamini_hochberg': 'fdr_bh', 'bonferroni': 'bonferroni', 'BY': 'fdr_by'}
        meth = correction_map.get(correction_method, 'fdr_bh')

        for g in groups:
            raw_pvals = []
            gene_info = []
            for i in range(min(n_genes, len(result['names'][g]))):
                gene = result['names'][g][i]
                raw_p = float(result['pvals'][g][i])
                if min_pct > 0 and g in pct_expr:
                    if pct_expr[g].get(gene, 0) < min_pct:
                        continue
                raw_pvals.append(max(raw_p, 1e-300))
                gene_info.append({
                    'gene': gene,
                    'logfc': round(float(result['logfoldchanges'][g][i]), 3),
                    'pval': raw_p,
                    'score': round(float(result['scores'][g][i]), 3),
                })

            if raw_pvals:
                reject, padj, _, _ = multipletests(raw_pvals, method=meth)
                for j, info in enumerate(gene_info):
                    info['pval_adj'] = float(padj[j])
                    info['cluster'] = g
                    deg_data.append(info)

        self.progress(70, "Generating volcano plots...")
        plots_dir = self.ensure_plots_dir()
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
            g_df['-log10(pval_adj)'] = -np.log10(g_df['pval_adj'].clip(lower=1e-300))
            fig = go.Figure()
            sig = (g_df['pval_adj'] < pval_cutoff) & (g_df['logfc'].abs() > logfc_cutoff)
            fig.add_trace(go.Scattergl(x=g_df.loc[sig, 'logfc'], y=g_df.loc[sig, '-log10(pval_adj)'],
                                       mode='markers', marker=dict(color='#e53935', size=4), name='Significant'))
            fig.add_trace(go.Scattergl(x=g_df.loc[~sig, 'logfc'], y=g_df.loc[~sig, '-log10(pval_adj)'],
                                       mode='markers', marker=dict(color='#9e9e9e', size=3), name='Not significant'))
            # Gene annotations on volcano
            annotate_genes = set()
            if volcano_genes_str:
                annotate_genes = {g.strip() for g in volcano_genes_str.replace('\n', ',').split(',') if g.strip()}
            if volcano_top_n > 0:
                top_sig = g_df[sig].nsmallest(volcano_top_n, 'pval_adj')
                annotate_genes.update(top_sig['gene'].tolist())
            for gene_name in annotate_genes:
                gene_row = g_df[g_df['gene'] == gene_name]
                if not gene_row.empty:
                    row = gene_row.iloc[0]
                    fig.add_annotation(x=row['logfc'], y=row['-log10(pval_adj)'],
                                      text=gene_name, showarrow=True, arrowhead=2, ax=20, ay=-20,
                                      font=dict(size=9))
            fig.update_layout(title=f'Volcano Plot: {first_group}', xaxis_title='Log2 FC', yaxis_title='-log10(padj)',
                             plot_bgcolor='white', width=600, height=400)
            result_files.append(self.save_plotly_json(fig, plots_dir, f'deg_volcano_{first_group}.json', 'volcano', f'Volcano: {first_group}'))

        # DEG Dotplot
        if self.params.get('show_dotplot', True):
            try:
                custom_dotplot_str = self.params.get('custom_dotplot_genes', '').strip()
                if custom_dotplot_str:
                    top_genes_list = [g.strip() for g in custom_dotplot_str.replace('\n', ',').split(',') if g.strip()]
                    not_found = [g for g in top_genes_list if g not in adata.var_names]
                    top_genes_list = [g for g in top_genes_list if g in adata.var_names]
                    dotplot_label = '自定义基因 Dotplot'
                else:
                    top_genes_list = []
                    not_found = []
                    for g in groups:
                        group_df = sc.get.rank_genes_groups_df(adata, group=g)
                        top_genes_list.extend(group_df.head(5)['names'].tolist())
                    top_genes_list = list(dict.fromkeys(top_genes_list))[:30]
                    dotplot_label = 'DEG Dotplot'

                if top_genes_list:
                    sc.tl.dendrogram(adata, groupby=groupby)
                    fig_dot = sc.pl.dotplot(adata, var_names=top_genes_list, groupby=groupby, return_fig=True)
                    import io, base64
                    buf = io.BytesIO()
                    fig_dot.savefig(buf, format='png', dpi=100, bbox_inches='tight')
                    import matplotlib.pyplot as plt
                    plt.close('all')
                    buf.seek(0)
                    img_b64 = base64.b64encode(buf.read()).decode()
                    fpath = os.path.join(plots_dir, 'deg_dotplot.json')
                    with open(fpath, 'w') as f:
                        json.dump({'data': [{'type': 'image', 'source': f'data:image/png;base64,{img_b64}', 'xref': 'paper', 'yref': 'paper', 'x': 0, 'y': 1, 'sizex': 1, 'sizey': 1, 'sizing': 'stretch'}], 'layout': {'width': 900, 'height': 500, 'title': 'DEG Dotplot'}}, f)
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'dotplot', 'label': dotplot_label})
            except Exception:
                pass

        # Gene expression UMAP
        plot_genes_umap = self.params.get('plot_genes_umap', '').strip()
        if plot_genes_umap and 'X_umap' in adata.obsm:
            gene_list = [g.strip() for g in plot_genes_umap.split(',') if g.strip() and g.strip() in adata.var_names]
            for gene in gene_list[:5]:  # Limit to 5 genes
                fig_gene = umap_scatter(adata, gene, title=f'{gene} Expression')
                fpath = os.path.join(plots_dir, f'deg_gene_umap_{gene}.json')
                with open(fpath, 'w') as f: f.write(json.dumps(fig_gene))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'{gene} Expression'})

        # Export full DEG results
        all_deg = []
        for g in groups:
            df_g = sc.get.rank_genes_groups_df(adata, group=g)
            df_g['cluster'] = g
            all_deg.append(df_g)
        if all_deg:
            all_deg_df = pd.concat(all_deg)
            full_csv = os.path.join(results_dir, 'sc_deg_full_results.csv')
            all_deg_df.to_csv(full_csv, index=False)
            result_files.append({'file_path': full_csv, 'file_type': 'csv', 'category': 'table', 'label': '完整 DEG 结果'})

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'deg')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_groups': len(groups),
                'groups': list(groups),
                'method': method,
                'reference': reference,
                'total_deg_genes': len(deg_data),
                'custom_dotplot_genes': self.params.get('custom_dotplot_genes', '').strip() or None,
            }
        }
