import os
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


def _ensure_gsea_term_column(result_df):
    """Preserve pathway names stored in a GSEA result index."""
    result = result_df.copy()
    if 'Term' not in result.columns:
        index_name = result.index.name or 'index'
        result = result.reset_index().rename(columns={index_name: 'Term'})
    return result


ONTOLOGY_COLORS = {
    'BP': '#FDBE85',
    'MF': '#B8A9D1',
    'KEGG': '#7BC77B',
    'OTHER': '#9CB8D8',
}


def _enrichment_ontology(value, database=''):
    """Normalize database labels to the compact ontology legend used in plots."""
    text = str(value or database or '').upper()
    if 'KEGG' in text:
        return 'KEGG'
    if 'MF' in text or 'MOLECULAR FUNCTION' in text:
        return 'MF'
    if 'BP' in text or 'BIOLOGICAL PROCESS' in text:
        return 'BP'
    return 'OTHER'


def _enrichment_count(value):
    """Extract the leading overlap count from values such as ``9/132``."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return 1
    text = str(value).strip()
    try:
        return max(1, int(float(text.split('/', 1)[0])))
    except (TypeError, ValueError):
        return 1


def _enrichment_figure(result_df, title, database='GO_BP', score_column=None):
    """Create a clean horizontal enrichment chart similar to the supplied PDF.

    Bars encode ``-log10(adjusted p-value)``, the left bubbles encode overlap
    count, and fill colors encode BP/MF/KEGG.  The figure is deliberately native
    Matplotlib so it exports identically to high-DPI PNG and SVG.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    df = result_df.copy()
    if df.empty:
        return None
    term_col = next((c for c in ('Term', 'Description', 'term', 'description') if c in df.columns), None)
    if term_col is None:
        df = df.reset_index().rename(columns={df.index.name or 'index': 'Term'})
        term_col = 'Term'
    p_col = next((c for c in ('Adjusted P-value', 'Adjusted p-value', 'p.adjust', 'FDR', 'fdr', 'P-value', 'pvalue') if c in df.columns), None)
    if score_column and score_column in df.columns:
        score = pd.to_numeric(df[score_column], errors='coerce').abs()
        x_label = 'Enrichment score'
    elif p_col:
        score = -np.log10(pd.to_numeric(df[p_col], errors='coerce').clip(lower=1e-300))
        x_label = '-log10(adjusted P-value)'
    else:
        score = pd.Series(np.arange(len(df), 0, -1), index=df.index, dtype=float)
        x_label = 'Enrichment score'
    df = df.assign(_score=score).dropna(subset=['_score']).sort_values('_score', ascending=False)
    if df.empty:
        return None

    ontology_col = next((c for c in ('Ontology', 'ontology', 'Gene_set', 'gene_set', 'database') if c in df.columns), None)
    ontologies = [
        _enrichment_ontology(row[ontology_col] if ontology_col else database, database)
        for _, row in df.iterrows()
    ]
    count_col = next((c for c in ('Overlap', 'Count', 'count', 'Gene Count', 'gene_count', 'setSize', 'size') if c in df.columns), None)
    counts = [_enrichment_count(row[count_col]) if count_col else 1 for _, row in df.iterrows()]
    terms = [str(value) for value in df[term_col].tolist()]
    values = df['_score'].astype(float).to_numpy()
    y = np.arange(len(df))
    colors = [ONTOLOGY_COLORS.get(item, ONTOLOGY_COLORS['OTHER']) for item in ontologies]

    height = max(4.8, min(13.0, 1.25 + 0.43 * len(df)))
    fig, ax = plt.subplots(figsize=(9.2, height), dpi=150)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.barh(y, values, color=colors, edgecolor='none', height=0.82, alpha=0.96, zorder=1)
    ax.set_yticks(y, [''] * len(y))
    ax.invert_yaxis()
    ax.set_xlabel(x_label, fontsize=12, color='#142a8b', labelpad=12)
    ax.set_ylabel('Description', fontsize=12, color='#142a8b', labelpad=18)
    fig.suptitle(title, x=0.03, y=0.98, ha='left', va='top',
                 fontsize=17, fontweight='semibold', color='#111827')
    ax.grid(axis='x', color='#d8dee9', linewidth=0.65, alpha=0.55, zorder=0)
    ax.tick_params(axis='y', labelsize=10, length=0, pad=7, colors='#111827')
    ax.tick_params(axis='x', labelsize=10, colors='#374151', width=0.6)
    for spine in ax.spines.values():
        spine.set_visible(False)

    max_score = max(float(values.max()), 1.0)
    ax.set_xlim(-max(0.9, max_score * 0.12), max_score * 1.08)
    bubble_sizes = [65 + 30 * np.sqrt(max(1, count)) for count in counts]
    ax.scatter(np.full(len(y), -max_score * 0.045), y, s=bubble_sizes,
               c=colors, edgecolors='#111827', linewidths=0.8, zorder=3)
    for yi, count in zip(y, counts):
        ax.text(-max_score * 0.045, yi, str(count), ha='center', va='center',
                fontsize=8.5, color='#111827', zorder=4)
    for yi, term in zip(y, terms):
        ax.text(max_score * 0.012, yi, term, ha='left', va='center',
                fontsize=10, color='#111827', zorder=2)

    present_ontologies = [item for item in ('BP', 'MF', 'KEGG', 'OTHER') if item in ontologies]
    ontology_handles = [Patch(facecolor=ONTOLOGY_COLORS[item], edgecolor='none', label=item)
                        for item in present_ontologies]
    ontology_legend = fig.legend(handles=ontology_handles, title='ONTOLOGY',
                                 loc='upper left', bbox_to_anchor=(0.82, 0.88),
                                 bbox_transform=fig.transFigure,
                                 frameon=False, fontsize=10, title_fontsize=11)
    legend_counts = sorted(set(counts))[:4]
    if legend_counts:
        count_handles = [Line2D([0], [0], marker='o', linestyle='none',
                                markerfacecolor='#111827', markeredgecolor='#111827',
                                markersize=5 + 2.2 * np.sqrt(max(1, value)),
                                label=str(value)) for value in legend_counts]
        fig.legend(handles=count_handles, title='Count', loc='upper left',
                   bbox_to_anchor=(0.82, 0.54), bbox_transform=fig.transFigure,
                   frameon=False, fontsize=10, title_fontsize=11)
    fig.subplots_adjust(left=0.07, right=0.75, top=0.88, bottom=0.13)
    return fig


class BulkEnrichmentAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_enrichment"
    DISPLAY_NAME = "通路富集分析"
    DESCRIPTION = "GO/KEGG/WikiPathways 通路富集分析（ORA / GSEA）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import omicverse as ov

        method = self.params.get('method', 'ORA')
        database = self.params.get('database', 'GO_BP')
        organism = self.params.get('organism', 'Human')
        pvalue_cutoff = float(self.params.get('pvalue_cutoff', 0.05))
        top_n = int(self.params.get('top_n', 20))
        input_source = self.params.get('input_source', '')
        split_direction = self.params.get('split_direction', False) in (True, 'true', 'on', '1')
        custom_genes_str = self.params.get('custom_genes', '').strip()

        self.progress(5, "加载基因列表...")

        # Validate input_source path
        if input_source:
            real_input = os.path.realpath(input_source)
            real_project = os.path.realpath(self.project_dir)
            if not real_input.startswith(real_project + os.sep) and real_input != real_project:
                raise ValueError("input_source 必须在项目目录内")

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
        except Exception as e:
            import logging
            logging.getLogger(__name__).warning(f"通路数据库下载失败: {e}")

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

        if not os.path.exists(db_path):
            raise FileNotFoundError(f"基因集数据库文件不存在: {db_path}")

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
                        fig_enrichment = _enrichment_figure(
                            top_enr_dir,
                            title=f'{direction}_enrich_result',
                            database=database,
                        )
                        if fig_enrichment is not None:
                            result_files.extend(self.save_matplotlib_figure(
                                fig_enrichment, plots_dir,
                                f'enrichment_ora{suffix}.png', 'enrichment',
                                f'ORA 富集图 ({direction})',
                            ))
                            import matplotlib.pyplot as plt
                            plt.close(fig_enrichment)

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
                    fig_enrichment = _enrichment_figure(
                        top_enr,
                        title=f'{database} enrich result',
                        database=database,
                    )
                    if fig_enrichment is not None:
                        result_files.extend(self.save_matplotlib_figure(
                            fig_enrichment, plots_dir, 'enrichment_ora.png',
                            'enrichment', 'ORA 富集图',
                        ))
                        import matplotlib.pyplot as plt
                        plt.close(fig_enrichment)

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

            self.progress(71, "处理 GSEA 结果...")
            enr = _ensure_gsea_term_column(pre_res.res2d)
            # omicverse 将通路名保存在索引中；导出 index=False 前显式保留为 Term 列。
            enr_sig = enr[enr['fdr'] < pvalue_cutoff].copy()

            csv_path = os.path.join(results_dir, 'enrichment_gsea_results.csv')
            enr.to_csv(csv_path, index=False)
            result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': 'GSEA 富集结果'})

            self.progress(75, "生成 GSEA 图表...")
            if len(enr_sig) > 0:
                top_gsea = enr_sig.head(top_n)
                fig_enrichment = _enrichment_figure(
                    top_gsea,
                    title='GSEA enrich result',
                    database=database,
                    score_column='nes' if 'nes' in top_gsea.columns else None,
                )
                if fig_enrichment is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        fig_enrichment, plots_dir, 'enrichment_gsea.png',
                        'enrichment', 'GSEA 富集图',
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_enrichment)

            n_sig = len(enr_sig)
        else:
            raise ValueError(f"不支持的富集方法: {method}")

        self.progress(90, "保存输出...")

        # Clean up temporary enrichment directories
        import shutil
        for tmp_suffix in ['_tmp', '_up_tmp', '_down_tmp', '_gsea_tmp']:
            tmp_dir = os.path.join(self.project_dir, f'enrichr{tmp_suffix}')
            if os.path.isdir(tmp_dir):
                shutil.rmtree(tmp_dir, ignore_errors=True)

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
