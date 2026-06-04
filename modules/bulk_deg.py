import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkDEGAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_deg"
    DISPLAY_NAME = "Bulk 差异表达分析"
    DESCRIPTION = "组间差异表达基因检测：t-test / Mann-Whitney / DESeq2 风格"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from scipy import stats
        from modules.visualization import scatter_plot

        self.progress(5, "加载数据...")
        if input_path.endswith('.csv') or input_path.endswith('.txt'):
            df = pd.read_csv(input_path, sep=None if input_path.endswith('.csv') else '\t', index_col=0)
            adata = sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
        else:
            adata = sc.read_h5ad(input_path)

        groupby = self.params.get('groupby', 'condition')
        group1 = self.params.get('group1', '')
        group2 = self.params.get('group2', '')
        method = self.params.get('method', 't-test')
        fc_threshold = float(self.params.get('fc_threshold', 2.0))
        pval_threshold = float(self.params.get('pval_threshold', 0.05))
        top_n = int(self.params.get('top_n', 20))

        self.progress(15, "解析分组信息...")

        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        counts = counts.astype(float)
        gene_names = adata.var_names.tolist()

        if groupby in adata.obs.columns:
            groups = adata.obs[groupby].unique().tolist()
            if not group1 or group1 not in groups:
                group1 = groups[0]
            if not group2 or group2 not in groups:
                group2 = groups[1] if len(groups) > 1 else groups[0]
            mask1 = adata.obs[groupby] == group1
            mask2 = adata.obs[groupby] == group2
            data1 = counts[mask1.values]
            data2 = counts[mask2.values]
        else:
            n = counts.shape[0]
            half = n // 2
            group1, group2 = "Group1", "Group2"
            data1, data2 = counts[:half], counts[half:]
            mask1 = pd.Series([True]*half + [False]*(n-half))
            mask2 = pd.Series([False]*half + [True]*(n-half))

        self.progress(30, f"差异分析：{method}...")
        n_genes = counts.shape[1]
        log2fc = np.zeros(n_genes)
        pvalues = np.zeros(n_genes)

        mean1 = np.mean(data1, axis=0)
        mean2 = np.mean(data2, axis=0)
        mean1_safe = np.where(mean1 > 0, mean1, 1e-10)
        mean2_safe = np.where(mean2 > 0, mean2, 1e-10)
        log2fc = np.log2(mean1_safe / mean2_safe)

        for i in range(n_genes):
            if method == 't-test':
                _, p = stats.ttest_ind(data1[:, i], data2[:, i], equal_var=False)
            else:
                _, p = stats.mannwhitneyu(data1[:, i], data2[:, i], alternative='two-sided')
            pvalues[i] = p if not np.isnan(p) else 1.0

        from statsmodels.stats.multitest import multipletests
        try:
            _, padj, _, _ = multipletests(pvalues, method='fdr_bh')
        except Exception:
            padj = pvalues

        self.progress(60, "生成结果表...")

        def classify(lfc, padj_val):
            if padj_val < pval_threshold and lfc >= np.log2(fc_threshold):
                return 'Up'
            elif padj_val < pval_threshold and lfc <= -np.log2(fc_threshold):
                return 'Down'
            return 'NS'

        regulation = [classify(log2fc[i], padj[i]) for i in range(n_genes)]
        deg_df = pd.DataFrame({
            'gene': gene_names,
            'log2FC': np.round(log2fc, 4),
            'pvalue': pvalues,
            'padj': padj,
            'mean_group1': np.round(mean1, 2),
            'mean_group2': np.round(mean2, 2),
            'regulation': regulation
        })
        deg_df = deg_df.sort_values('padj')

        n_up = sum(1 for r in regulation if r == 'Up')
        n_down = sum(1 for r in regulation if r == 'Down')

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        self.progress(70, "生成火山图...")
        neg_log_padj = -np.log10(padj + 1e-300)
        color_map = {'Up': '#e53935', 'Down': '#1a237e', 'NS': '#bdbdbd'}
        fig_vol = go.Figure()
        for reg in ['NS', 'Up', 'Down']:
            idx = [i for i in range(n_genes) if regulation[i] == reg]
            fig_vol.add_trace(go.Scattergl(
                x=log2fc[idx], y=neg_log_padj[idx], mode='markers',
                marker=dict(color=color_map[reg], size=5, opacity=0.7),
                name=f'{reg} ({len(idx)})',
                text=[gene_names[i] for i in idx], hovertemplate='%{text}<br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}'
            ))
        fig_vol.add_hline(y=-np.log10(pval_threshold), line_dash='dash', line_color='gray')
        fig_vol.add_vline(x=np.log2(fc_threshold), line_dash='dash', line_color='gray')
        fig_vol.add_vline(x=-np.log2(fc_threshold), line_dash='dash', line_color='gray')
        fig_vol.update_layout(
            title=f'火山图 ({group1} vs {group2})',
            xaxis_title='log2(Fold Change)', yaxis_title='-log10(padj)',
            plot_bgcolor='white', width=700, height=500
        )
        fpath = os.path.join(plots_dir, 'bulk_deg_volcano.json')
        with open(fpath, 'w') as f: json.dump(json.loads(fig_vol.to_json()), f)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'volcano', 'label': '火山图'})

        self.progress(80, "生成 MA 图...")
        avg_expr = (mean1 + mean2) / 2
        fig_ma = go.Figure()
        for reg in ['NS', 'Up', 'Down']:
            idx = [i for i in range(n_genes) if regulation[i] == reg]
            fig_ma.add_trace(go.Scattergl(
                x=np.log2(avg_expr[idx] + 1), y=log2fc[idx], mode='markers',
                marker=dict(color=color_map[reg], size=5, opacity=0.7),
                name=f'{reg}', text=[gene_names[i] for i in idx],
                hovertemplate='%{text}<br>AvgExpr: %{x:.2f}<br>log2FC: %{y:.2f}'
            ))
        fig_ma.add_hline(y=0, line_dash='solid', line_color='gray', line_width=0.5)
        fig_ma.update_layout(
            title='MA 图', xaxis_title='log2(Average Expression)', yaxis_title='log2(Fold Change)',
            plot_bgcolor='white', width=700, height=500
        )
        fpath = os.path.join(plots_dir, 'bulk_deg_ma.json')
        with open(fpath, 'w') as f: json.dump(json.loads(fig_ma.to_json()), f)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'ma', 'label': 'MA 图'})

        self.progress(88, "保存差异基因 CSV...")
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        csv_path = os.path.join(results_dir, 'bulk_deg_results.csv')
        deg_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '差异表达基因列表'})

        top_genes = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
        top_csv = os.path.join(results_dir, 'bulk_deg_top_genes.csv')
        top_genes.to_csv(top_csv, index=False)
        result_files.append({'file_path': top_csv, 'file_type': 'csv', 'category': 'table', 'label': f'Top {top_n} 差异基因'})

        self.progress(95, "保存 h5ad...")
        adata.obs['group'] = adata.obs[groupby] if groupby in adata.obs.columns else 'unknown'
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_deg_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'comparison': f'{group1} vs {group2}',
                'method': method,
                'n_genes_total': n_genes,
                'n_up': n_up,
                'n_down': n_down,
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
            }
        }
