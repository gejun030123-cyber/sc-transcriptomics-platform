import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkEnrichmentAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_enrichment"
    DISPLAY_NAME = "通路富集分析"
    DESCRIPTION = "GO/KEGG/WikiPathways 通路富集分析（ORA / GSEA）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import omicverse as ov
        import plotly.graph_objects as go

        method = self.params.get('method', 'ORA')
        database = self.params.get('database', 'GO_BP')
        organism = self.params.get('organism', 'Human')
        pvalue_cutoff = float(self.params.get('pvalue_cutoff', 0.05))
        top_n = int(self.params.get('top_n', 20))
        input_source = self.params.get('input_source', '')
        split_direction = self.params.get('split_direction', False) in (True, 'true', 'on', '1')
        custom_genes_str = self.params.get('custom_genes', '').strip()

        self.progress(5, "加载基因列表...")

        # Load gene list from DEG CSV result
        deg_genes = []
        gene_rnk = None

        if custom_genes_str:
            deg_genes = [g.strip() for g in custom_genes_str.replace('\n', ',').split(',') if g.strip()]
            gene_rnk = None
        elif input_source and os.path.exists(input_source):
            deg_df = pd.read_csv(input_source)
            if 'regulation' in deg_df.columns:
                deg_genes = deg_df[deg_df['regulation'] != 'NS']['gene'].tolist()
                if 'log2FC' in deg_df.columns:
                    gene_rnk = deg_df[['gene', 'log2FC']].dropna().sort_values('log2FC', ascending=False)
                    gene_rnk.columns = ['gene_name', 'rank']
            elif 'gene' in deg_df.columns:
                deg_genes = deg_df['gene'].tolist()

        if not deg_genes and method == 'ORA':
            raise ValueError("未找到差异基因列表。请先运行 DEG 分析，并在 input_source 中指定结果 CSV 文件路径。")

        self.progress(15, f"加载 {database} 基因集数据库...")

        # Map database name to file
        db_map = {
            'GO_BP': 'GO_Biological_Process_2021',
            'GO_MF': 'GO_Molecular_Function_2021',
            'GO_CC': 'GO_Cellular_Component_2021',
            'KEGG': 'KEGG_2021_Human' if organism == 'Human' else 'KEGG_2019_Mouse',
            'WikiPathways': 'WikiPathways_2019_Human' if organism == 'Human' else 'WikiPathways_2019_Mouse',
            'Reactome': 'Reactome_2022',
        }

        db_filename = db_map.get(database, 'GO_Biological_Process_2021')
        organism_lower = organism.lower()

        # Download pathway databases if not present
        try:
            ov.utils.download_pathway_database()
        except Exception:
            pass

        db_path = f'genesets/{db_filename}.txt'
        if not os.path.exists(db_path):
            alt_paths = [
                f'genesets/{db_filename}_{organism}.txt',
                f'genesets/{db_filename}_{organism_lower}.txt',
            ]
            for p in alt_paths:
                if os.path.exists(p):
                    db_path = p
                    break

        pathways_dict = ov.utils.geneset_prepare(db_path, organism=organism)

        self.progress(30, f"运行 {method} 富集分析...")

        result_files = []
        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        n_sig = 0

        if method == 'ORA':
            if split_direction and not custom_genes_str and input_source and os.path.exists(input_source):
                # Split direction: run ORA separately for Up and Down genes
                deg_df_split = pd.read_csv(input_source)
                all_split_results = []

                for direction in ['Up', 'Down']:
                    direction_genes = deg_df_split[deg_df_split['regulation'] == direction]['gene'].tolist()
                    if not direction_genes:
                        continue

                    self.progress(45 if direction == 'Up' else 60, f"运行 {direction} 基因 ORA...")
                    enr_dir = ov.bulk.geneset_enrichment(
                        gene_list=direction_genes,
                        pathways_dict=pathways_dict,
                        pvalue_type='adjust',
                        pvalue_threshold=pvalue_cutoff,
                        organism=organism_lower,
                        outdir=os.path.join(self.project_dir, f'enrichr_{direction.lower()}_tmp')
                    )

                    enr_dir['direction'] = direction
                    all_split_results.append(enr_dir)

                    suffix = '_up' if direction == 'Up' else '_down'
                    top_enr_dir = enr_dir.head(top_n)
                    if len(top_enr_dir) > 0:
                        x_col = 'Fractions' if 'Fractions' in top_enr_dir.columns else ('Odds Ratio' if 'Odds Ratio' in top_enr_dir.columns else None)
                        term_col = 'Term' if 'Term' in top_enr_dir.columns else None

                        fig_bubble = go.Figure()
                        x_vals = top_enr_dir[x_col].values if x_col else list(range(len(top_enr_dir)))
                        y_vals = top_enr_dir[term_col].tolist() if term_col else top_enr_dir.index.tolist()

                        if 'Overlap' in top_enr_dir.columns:
                            sizes = top_enr_dir['Overlap'].apply(lambda x: int(str(x).split('/')[0]) * 3 + 5).values
                        else:
                            sizes = [10] * len(top_enr_dir)

                        if 'P-value' in top_enr_dir.columns:
                            colors = -np.log10(top_enr_dir['P-value'].clip(lower=1e-300).values)
                        else:
                            colors = [1] * len(top_enr_dir)

                        fig_bubble.add_trace(go.Scatter(
                            x=x_vals, y=y_vals, mode='markers',
                            marker=dict(size=sizes, color=colors, colorscale='YlOrRd', showscale=True,
                                       colorbar=dict(title='-log10(p)')),
                            hovertemplate='%{y}<br>-log10(p): %{marker.color:.1f}<extra></extra>'
                        ))
                        fig_bubble.update_layout(
                            title=f'{database} ORA 富集分析 ({organism}) - {direction} genes',
                            xaxis_title='Gene Fraction' if x_col else 'Index',
                            yaxis=dict(autorange='reversed'),
                            plot_bgcolor='white', width=800, height=max(400, top_n * 25 + 100)
                        )
                        fpath = os.path.join(plots_dir, f'enrichment_ora_bubble{suffix}.json')
                        with open(fpath, 'w') as f: f.write(fig_bubble.to_json(engine="json"))
                        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': f'ORA 气泡图 ({direction})'})

                if all_split_results:
                    combined = pd.concat(all_split_results, ignore_index=True)
                    csv_path = os.path.join(results_dir, 'enrichment_split_results.csv')
                    combined.to_csv(csv_path, index=False)
                    result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '分方向 ORA 富集结果'})
                    n_sig = len(combined[combined['P-value'] < pvalue_cutoff]) if 'P-value' in combined.columns else 0

            else:
                enr = ov.bulk.geneset_enrichment(
                    gene_list=deg_genes,
                    pathways_dict=pathways_dict,
                    pvalue_type='adjust',
                    pvalue_threshold=pvalue_cutoff,
                    organism=organism_lower,
                    outdir=os.path.join(self.project_dir, 'enrichr_tmp')
                )

                self.progress(60, "保存结果表...")
                csv_path = os.path.join(results_dir, 'enrichment_ora_results.csv')
                enr.to_csv(csv_path, index=False)
                result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': 'ORA 富集结果'})

                self.progress(70, "生成气泡图...")
                top_enr = enr.head(top_n)
                if len(top_enr) > 0:
                    # Determine x-axis column
                    x_col = 'Fractions' if 'Fractions' in top_enr.columns else ('Odds Ratio' if 'Odds Ratio' in top_enr.columns else None)
                    term_col = 'Term' if 'Term' in top_enr.columns else None

                    fig_bubble = go.Figure()
                    x_vals = top_enr[x_col].values if x_col else list(range(len(top_enr)))
                    y_vals = top_enr[term_col].tolist() if term_col else top_enr.index.tolist()

                    # Size by overlap count
                    if 'Overlap' in top_enr.columns:
                        sizes = top_enr['Overlap'].apply(lambda x: int(str(x).split('/')[0]) * 3 + 5).values
                    else:
                        sizes = [10] * len(top_enr)

                    # Color by p-value
                    if 'P-value' in top_enr.columns:
                        colors = -np.log10(top_enr['P-value'].clip(lower=1e-300).values)
                    else:
                        colors = [1] * len(top_enr)

                    fig_bubble.add_trace(go.Scatter(
                        x=x_vals, y=y_vals, mode='markers',
                        marker=dict(size=sizes, color=colors, colorscale='YlOrRd', showscale=True,
                                   colorbar=dict(title='-log10(p)')),
                        hovertemplate='%{y}<br>-log10(p): %{marker.color:.1f}<extra></extra>'
                    ))
                    fig_bubble.update_layout(
                        title=f'{database} ORA 富集分析 ({organism})',
                        xaxis_title='Gene Fraction' if x_col else 'Index',
                        yaxis=dict(autorange='reversed'),
                        plot_bgcolor='white', width=800, height=max(400, top_n * 25 + 100)
                    )
                    fpath = os.path.join(plots_dir, 'enrichment_ora_bubble.json')
                    with open(fpath, 'w') as f: f.write(fig_bubble.to_json(engine="json"))
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': 'ORA 气泡图'})

                    # Bar chart
                    if 'P-value' in top_enr.columns:
                        fig_bar = go.Figure()
                        fig_bar.add_trace(go.Bar(
                            x=-np.log10(top_enr['P-value'].clip(lower=1e-300).values),
                            y=y_vals, orientation='h', marker_color='#e53935'
                        ))
                        fig_bar.update_layout(
                            title=f'Top {top_n} 富集通路',
                            xaxis_title='-log10(P-value)',
                            yaxis=dict(autorange='reversed'),
                            plot_bgcolor='white', width=700, height=max(400, top_n * 25 + 100)
                        )
                        fpath = os.path.join(plots_dir, 'enrichment_ora_bar.json')
                        with open(fpath, 'w') as f: f.write(fig_bar.to_json(engine="json"))
                        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': 'ORA 条形图'})

                n_sig = len(enr[enr['P-value'] < pvalue_cutoff]) if 'P-value' in enr.columns else len(enr)

        elif method == 'GSEA':
            if gene_rnk is None or len(gene_rnk) == 0:
                raise ValueError("GSEA 需要排序的基因列表。请确保 DEG 结果包含 log2FC 列。")

            pre_res = ov.bulk.geneset_enrichment_GSEA(
                gene_rnk=gene_rnk,
                pathways_dict=pathways_dict,
                processes=4,
                permutation_num=1000,
                outdir=os.path.join(self.project_dir, 'enrichr_gsea_tmp')
            )

            self.progress(60, "处理 GSEA 结果...")
            enr = pre_res.res2d
            enr_sig = enr[enr['fdr'] < pvalue_cutoff].copy()

            csv_path = os.path.join(results_dir, 'enrichment_gsea_results.csv')
            enr.to_csv(csv_path, index=False)
            result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': 'GSEA 富集结果'})

            self.progress(70, "生成 GSEA 图表...")
            if len(enr_sig) > 0:
                top_gsea = enr_sig.head(top_n)
                fig_bar = go.Figure()
                fig_bar.add_trace(go.Bar(
                    x=top_gsea['nes'].values,
                    y=top_gsea.index.tolist(),
                    orientation='h',
                    marker_color=['#e53935' if v > 0 else '#1a237e' for v in top_gsea['nes'].values]
                ))
                fig_bar.update_layout(
                    title=f'Top {top_n} GSEA 通路 (NES)',
                    xaxis_title='Normalized Enrichment Score',
                    yaxis=dict(autorange='reversed'),
                    plot_bgcolor='white', width=700, height=max(400, top_n * 25 + 100)
                )
                fpath = os.path.join(plots_dir, 'enrichment_gsea_nes.json')
                with open(fpath, 'w') as f: f.write(fig_bar.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': 'GSEA NES 条形图'})

            n_sig = len(enr_sig)
        else:
            raise ValueError(f"不支持的富集方法: {method}")

        self.progress(90, "保存输出...")
        self.progress(100, "完成")

        return {
            'output_adata': None,
            'result_files': result_files,
            'summary': {
                'method': method,
                'database': database,
                'organism': organism,
                'n_input_genes': len(deg_genes) if method == 'ORA' else (len(gene_rnk) if gene_rnk is not None else 0),
                'n_significant': n_sig,
                'pvalue_cutoff': pvalue_cutoff,
                'input_mode': 'custom_genes' if custom_genes_str else ('split_direction' if split_direction else 'standard'),
            }
        }
