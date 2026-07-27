from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping


VOLCANO_COLOR_MAP = {
    'Up': '#B64342',       # Nature 风格红色：上调
    'Down': '#0F4D92',     # Nature 风格深蓝：下调
    'NS': '#98A2B3',       # 中性灰：不显著
}


def classify_volcano_regulation(logfc, padj, logfc_cutoff=1.0, pval_cutoff=0.05):
    """Classify genes for the three-colour volcano plot.

    The same adjusted-p-value and absolute-logFC thresholds are used by the
    DEG count plot, CSV interpretation, and volcano colours:
    ``Up`` / ``Down`` / ``NS``.
    """
    import numpy as np

    logfc = np.asarray(logfc, dtype=float)
    padj = np.asarray(padj, dtype=float)
    labels = np.full(logfc.shape, 'NS', dtype=object)
    significant = np.isfinite(logfc) & np.isfinite(padj) & (padj < pval_cutoff)
    labels[significant & (logfc > logfc_cutoff)] = 'Up'
    labels[significant & (logfc < -logfc_cutoff)] = 'Down'
    return labels

class DEGAnalysis(BaseAnalysis):
    MODULE_NAME = "deg"
    DISPLAY_NAME = "差异表达"
    DESCRIPTION = "差异表达基因分析（Wilcoxon 检验）"
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        import scanpy as sc
        import pandas as pd
        import os
        from modules.native_figures import (
            heatmap_figure, umap_panel_figure, scatter_figure, grouped_bar_figure,
        )

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        requested_groupby = str(self.params.get('groupby', 'celltype') or '').strip()
        method = self.params.get('method', 'wilcoxon')
        n_genes = int(self.params.get('n_genes', 20))
        reference = self.params.get('reference', 'rest')
        pval_cutoff = float(self.params.get('pval_cutoff', 0.05))
        logfc_cutoff = float(self.params.get('logfc_cutoff', 1.0))
        min_pct = float(self.params.get('min_pct', 0.1))
        correction_method = self.params.get('correction_method', 'benjamini_hochberg')
        volcano_top_n = int(self.params.get('volcano_top_n', 10))
        volcano_genes_str = self.params.get('volcano_genes', '').strip()

        groupby, group_info = resolve_obs_grouping(
            adata, requested_groupby, fallbacks=['leiden'],
            max_categories=50, max_numeric_categories=20,
            require_multiple=True,
        )
        if groupby is None:
            raise ValueError(
                f"groupby '{requested_groupby}' 不是有效的分类分组列："
                f"{group_info.get('reason', '')}"
            )
        if groupby != requested_groupby:
            self.progress(-1, f"分组列已改用 '{groupby}'：{group_info.get('requested_reason', '')}")
        if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
            adata.obs[groupby] = adata.obs[groupby].astype(str).astype('category')

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
        if not deg_df.empty:
            deg_df['regulation'] = classify_volcano_regulation(
                deg_df['logfc'].to_numpy(),
                deg_df['pval_adj'].to_numpy(),
                logfc_cutoff=logfc_cutoff,
                pval_cutoff=pval_cutoff,
            )
        csv_path = os.path.join(results_dir, 'deg_results.csv')
        deg_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': f'Top {n_genes} DEG per Cluster'})

        marker_selection = {
            'cluster_markers': {}, 'genes': [], 'warnings': [],
            'n_data_driven': 0, 'n_classic_anchor': 0, 'n_relaxed': 0,
            'parameters': {},
        }
        custom_dotplot_str = self.params.get('custom_dotplot_genes', '').strip()
        if not custom_dotplot_str:
            try:
                from modules.annotation import select_cluster_marker_genes
                marker_selection = select_cluster_marker_genes(
                    adata,
                    groupby,
                    method=method,
                    n_rank_genes=max(1, min(100, int(n_genes))),
                    min_markers=1,
                    max_markers=5,
                    padj_cutoff=pval_cutoff,
                    min_pct=min_pct,
                    min_delta_pct=0.05,
                )
                if marker_selection.get('warnings'):
                    self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
            except Exception as exc:
                self.progress(-1, f'DEG marker dotplot 选择失败，回退到排名结果：{exc}')

        first_group = groups[0] if groups else None
        if first_group:
            g_mask = deg_df['cluster'] == first_group
            g_df = deg_df[g_mask].copy()
            g_df['-log10(pval_adj)'] = -np.log10(g_df['pval_adj'].clip(lower=1e-300))
            if 'regulation' not in g_df:
                g_df['regulation'] = classify_volcano_regulation(
                    g_df['logfc'].to_numpy(),
                    g_df['pval_adj'].to_numpy(),
                    logfc_cutoff=logfc_cutoff,
                    pval_cutoff=pval_cutoff,
                )
            sig = g_df['regulation'].isin(['Up', 'Down'])
            volcano_colors = g_df['regulation'].map(VOLCANO_COLOR_MAP).fillna(
                VOLCANO_COLOR_MAP['NS']
            ).to_numpy()
            fig = scatter_figure(
                g_df['logfc'].values, g_df['-log10(pval_adj)'].values,
                title=f'Volcano Plot: {first_group}', x_label='Log2 FC',
                y_label='-log10(padj)', colors=volcano_colors, size=11,
                alpha=0.78,
            )
            ax = fig.axes[0]
            ax.axvline(logfc_cutoff, color='#98A2B3', linestyle='--', linewidth=0.7)
            ax.axvline(-logfc_cutoff, color='#98A2B3', linestyle='--', linewidth=0.7)
            ax.axhline(-np.log10(pval_cutoff), color='#98A2B3', linestyle='--', linewidth=0.7)
            from matplotlib.lines import Line2D
            from matplotlib.patches import Patch

            regulation_counts = g_df['regulation'].value_counts()
            legend_handles = [
                Patch(facecolor=VOLCANO_COLOR_MAP['Up'], edgecolor='none',
                      label=f"Up-regulated (n={int(regulation_counts.get('Up', 0))})"),
                Patch(facecolor=VOLCANO_COLOR_MAP['Down'], edgecolor='none',
                      label=f"Down-regulated (n={int(regulation_counts.get('Down', 0))})"),
                Patch(facecolor=VOLCANO_COLOR_MAP['NS'], edgecolor='none',
                      label=f"Not significant (n={int(regulation_counts.get('NS', 0))})"),
                Line2D([0], [0], color='#98A2B3', linestyle='--', linewidth=0.8,
                       label='Thresholds'),
            ]
            ax.legend(handles=legend_handles, loc='upper right', frameon=False,
                      fontsize=7.5, handlelength=1.5, borderaxespad=0.5)
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
                    ax.annotate(gene_name, (row['logfc'], row['-log10(pval_adj)']),
                                xytext=(4, 4), textcoords='offset points', fontsize=7,
                                color='#111827')
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, f'deg_volcano_{first_group}.png', 'volcano',
                f'Volcano: {first_group}', formats=('png', 'svg'), dpi=300,
            ))

        if self.params.get('show_deg_counts_bar', True) and not deg_df.empty:
            try:
                group_order = [str(g) for g in groups]
                up_counts = []
                down_counts = []
                for g in group_order:
                    sub = deg_df[deg_df['cluster'].astype(str) == g]
                    sig = sub['pval_adj'] < pval_cutoff
                    up_counts.append(int((sig & (sub['logfc'] > logfc_cutoff)).sum()))
                    down_counts.append(int((sig & (sub['logfc'] < -logfc_cutoff)).sum()))
                fig_counts = grouped_bar_figure(
                    group_order, [('Up', up_counts), ('Down', down_counts)],
                    title=f'Significant DEG Counts (padj<{pval_cutoff}, |logFC|>{logfc_cutoff})',
                    x_label=groupby, y_label='Gene count', rotation=35,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_counts, plots_dir, 'deg_significant_counts_bar.png',
                    'bar', 'Significant DEG Counts', formats=('png', 'svg'), dpi=300,
                ))
            except Exception as e:
                self.progress(-1, f"DEG count summary plot failed: {e}")

        # DEG Dotplot
        if self.params.get('show_dotplot', True):
            try:
                if custom_dotplot_str:
                    top_genes_list = [g.strip() for g in custom_dotplot_str.replace('\n', ',').split(',') if g.strip()]
                    not_found = [g for g in top_genes_list if g not in adata.var_names]
                    top_genes_list = [g for g in top_genes_list if g in adata.var_names]
                    dotplot_label = '自定义基因 Dotplot'
                else:
                    top_genes_list = marker_selection.get('genes', [])[:30]
                    not_found = []
                    if not top_genes_list:
                        for g in groups:
                            group_df = sc.get.rank_genes_groups_df(adata, group=g)
                            top_genes_list.extend(group_df.head(5)['names'].tolist())
                        top_genes_list = list(dict.fromkeys(top_genes_list))[:30]
                    dotplot_label = 'DEG Dotplot'

                if top_genes_list:
                    sc.tl.dendrogram(adata, groupby=groupby)
                    fig_dot = sc.pl.dotplot(adata, var_names=top_genes_list, groupby=groupby, return_fig=True)
                    result_files.extend(self.save_matplotlib_figure(
                        fig_dot, plots_dir, 'deg_dotplot.png', 'dotplot', dotplot_label
                    ))
                    import matplotlib.pyplot as plt
                    plt.close('all')
            except Exception as e:
                self.progress(-1, f"DEG dotplot generation failed: {e}")

        # Gene expression UMAP
        plot_genes_umap = self.params.get('plot_genes_umap', '').strip()
        if plot_genes_umap and 'X_umap' in adata.obsm:
            gene_list = [g.strip() for g in plot_genes_umap.split(',') if g.strip() and g.strip() in adata.var_names]
            for gene in gene_list[:5]:  # Limit to 5 genes
                try:
                    fig_static = self.build_publication_umap(
                        adata, gene, title=f'{gene} Expression'
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_static, plots_dir, f'deg_gene_umap_{gene}.png', 'umap',
                        f'{gene} Expression'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_static)
                except Exception as exc:
                    self.progress(-1, f'{gene} 表达 UMAP 导出失败：{exc}')

        if self.params.get('show_top_marker_umap_panel', True) and 'X_umap' in adata.obsm and not deg_df.empty:
            try:
                marker_limit = int(self.params.get('top_marker_umap_genes', 6))
                marker_candidates = []
                for g in groups:
                    sub = deg_df[deg_df['cluster'].astype(str) == str(g)]
                    sig_sub = sub[(sub['pval_adj'] < pval_cutoff) & (sub['logfc'] > logfc_cutoff)]
                    genes_for_group = sig_sub['gene'].astype(str).tolist() or sub['gene'].astype(str).tolist()
                    marker_candidates.extend(genes_for_group[:1])
                marker_genes = [g for g in dict.fromkeys(marker_candidates) if g in adata.var_names][:max(1, marker_limit)]
                if marker_genes:
                    fig_marker_umap = umap_panel_figure(
                        adata, marker_genes, titles=marker_genes,
                        point_size=4, opacity=0.70, ncols=min(3, len(marker_genes)),
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_marker_umap, plots_dir, 'deg_top_marker_umap_panel.png',
                        'umap', 'Top Marker Feature UMAPs', formats=('png', 'svg'), dpi=300,
                    ))
            except Exception as e:
                self.progress(-1, f"Top marker UMAP panel failed: {e}")

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

        # Cluster marker heatmap
        if self.params.get('show_marker_heatmap', True) and not deg_df.empty:
            try:
                heatmap_top_n = int(self.params.get('marker_heatmap_top_n', 3))
                heatmap_genes = []
                for g in groups:
                    cluster_genes = deg_df.loc[deg_df['cluster'] == g, 'gene'].astype(str).tolist()
                    heatmap_genes.extend(cluster_genes[:max(1, heatmap_top_n)])
                heatmap_genes = [g for g in dict.fromkeys(heatmap_genes) if g in adata.var_names][:60]
                if heatmap_genes:
                    expr = adata[:, heatmap_genes].X
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray()
                    expr = np.asarray(expr, dtype=float)
                    group_labels = adata.obs[groupby].astype(str)
                    group_order = [str(g) for g in groups]
                    mean_matrix = []
                    for g in group_order:
                        mask = (group_labels == g).values
                        if mask.sum() == 0:
                            mean_matrix.append(np.zeros(len(heatmap_genes)))
                        else:
                            mean_matrix.append(expr[mask, :].mean(axis=0))
                    mean_matrix = np.asarray(mean_matrix).T
                    row_mean = mean_matrix.mean(axis=1, keepdims=True)
                    row_std = mean_matrix.std(axis=1, keepdims=True) + 1e-10
                    z = np.clip((mean_matrix - row_mean) / row_std, -3, 3)
                    fig_heat = heatmap_figure(
                        z, x_labels=group_order, y_labels=heatmap_genes,
                        title=f'Top Marker Heatmap by {groupby}',
                        x_label=groupby, y_label='Marker genes',
                        colorbar_label='Row z-score',
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_heat, plots_dir, 'deg_marker_heatmap.png',
                        'heatmap', 'Cluster Marker Heatmap', formats=('png', 'svg'), dpi=300,
                    ))
            except Exception as e:
                self.progress(-1, f"Marker heatmap generation failed: {e}")

        self.progress(90, "Saving output...")
        adata.uns['marker_selection'] = {
            'cluster_key': groupby,
            'used_in_dotplot': not bool(custom_dotplot_str),
            'cluster_markers': marker_selection.get('cluster_markers', {}),
            'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
            'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
            'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
            'warnings': marker_selection.get('warnings', []),
            'parameters': marker_selection.get('parameters', {}),
        }
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
                'requested_groupby': requested_groupby,
                'groupby': groupby,
                'marker_selection': {
                    'used_in_dotplot': not bool(custom_dotplot_str),
                    'cluster_markers': marker_selection.get('cluster_markers', {}),
                    'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
                    'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
                    'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
                    'warnings': marker_selection.get('warnings', []),
                    'parameters': marker_selection.get('parameters', {}),
                },
            }
        }
