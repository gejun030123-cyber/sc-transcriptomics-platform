import os
import json
from modules.base import BaseAnalysis

def _run_stat_test(ct_abs, test_type, n_permutations=1000):
    """对列联表运行指定统计检验，返回 (statistic, p_value)。"""
    from scipy.stats import chi2_contingency, fisher_exact
    import numpy as np

    if test_type == 'chi_square':
        chi2, pval, _, _ = chi2_contingency(ct_abs)
        return chi2, pval
    elif test_type == 'fisher_exact':
        if ct_abs.shape == (2, 2):
            stat, pval = fisher_exact(ct_abs.values)
            return stat, pval
        else:
            chi2, pval, _, _ = chi2_contingency(ct_abs)
            return chi2, pval
    elif test_type == 'permutation':
        observed_chi2, _, _, _ = chi2_contingency(ct_abs)
        count = 0
        import pandas as pd
        for _ in range(n_permutations):
            shuffled = ct_abs.copy()
            total = shuffled.values.sum()
            col_sums = shuffled.sum(axis=0).values
            row_sums = shuffled.sum(axis=1).values
            pvals = col_sums / total
            pvals = pvals / pvals.sum()
            perm_table = np.random.multinomial(row_sums[0], pvals).reshape(1, -1)
            for rs in row_sums[1:]:
                row = np.random.multinomial(rs, pvals)
                perm_table = np.vstack([perm_table, row])
            perm_df = pd.DataFrame(perm_table, index=ct_abs.index, columns=ct_abs.columns)
            perm_chi2, _, _, _ = chi2_contingency(perm_df)
            if perm_chi2 >= observed_chi2:
                count += 1
        pval = (count + 1) / (n_permutations + 1)
        return observed_chi2, pval
    else:
        chi2, pval, _, _ = chi2_contingency(ct_abs)
        return chi2, pval

class ProportionAnalysis(BaseAnalysis):
    MODULE_NAME = "proportion"
    DISPLAY_NAME = "细胞比例分析"
    DESCRIPTION = "各分组间的细胞比例差异分析"
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        import scanpy as sc
        import pandas as pd
        from scipy.stats import chi2_contingency
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        groupby = self.params.get('groupby', 'celltype')
        batch_key = self.params.get('batch_key', 'batch')
        stat_test = self.params.get('stat_test', 'chi_square')
        n_permutations = int(self.params.get('n_permutations', 1000))
        min_cells_per_group = int(self.params.get('min_cells_per_group', 10))

        if groupby not in adata.obs.columns:
            groupby = 'leiden'

        self.progress(30, "Computing cell proportions...")
        if min_cells_per_group > 0:
            group_counts = adata.obs[groupby].value_counts()
            valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()
            if len(valid_groups) < len(group_counts):
                removed = set(group_counts.index) - set(valid_groups)
                self.progress(-1, f"移除 {len(removed)} 个低细胞数组: {removed}")
                adata = adata[adata.obs[groupby].isin(valid_groups)].copy()
        ct = pd.crosstab(adata.obs[batch_key], adata.obs[groupby], normalize='index')
        ct_abs = pd.crosstab(adata.obs[batch_key], adata.obs[groupby])

        self.progress(50, "Running chi-squared test...")
        chi2, pval = _run_stat_test(ct_abs, stat_test, n_permutations)

        self.progress(65, "Generating proportion plots...")
        plots_dir = self.ensure_plots_dir()
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        fig = go.Figure()
        for col in ct.columns:
            fig.add_trace(go.Bar(name=str(col), x=ct.index.tolist(), y=ct[col].values))
        fig.update_layout(**self.get_plotly_layout(
            title='Cell Proportions by Group',
            barmode='stack',
            xaxis_title=batch_key, yaxis_title='Proportion',
        ))
        result_files.append(self.save_plotly_json(fig, plots_dir, 'proportion_stacked.json', 'bar', 'Cell Proportions (Stacked)'))

        fig2 = make_subplots(rows=1, cols=len(ct.index), subplot_titles=[str(x) for x in ct.index])
        for i, idx in enumerate(ct.index, 1):
            fig2.add_trace(go.Pie(labels=[str(x) for x in ct.columns], values=ct.loc[idx].values, hole=0.3),
                          row=1, col=i)
        fig2.update_layout(title='Cell Type Distribution per Batch', width=300 * len(ct.index), height=400)
        result_files.append(self.save_plotly_json(fig2, plots_dir, 'proportion_pie.json', 'pie', 'Cell Type Distribution'))

        ct_abs.to_csv(os.path.join(results_dir, 'cell_counts.csv'))
        ct.to_csv(os.path.join(results_dir, 'cell_proportions.csv'))
        result_files.append({'file_path': os.path.join(results_dir, 'cell_counts.csv'), 'file_type': 'csv', 'category': 'table', 'label': 'Cell Counts'})
        result_files.append({'file_path': os.path.join(results_dir, 'cell_proportions.csv'), 'file_type': 'csv', 'category': 'table', 'label': 'Cell Proportions'})

        # Pairwise group comparison
        compare_groups_str = self.params.get('compare_groups', '').strip()
        compare_pairs = []
        if compare_groups_str:
            for item in compare_groups_str.replace('\n', ';').split(';'):
                item = item.strip()
                if '-vs-' in item:
                    parts = item.split('-vs-')
                    if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                        compare_pairs.append((parts[0].strip(), parts[1].strip()))

        for group_a, group_b in compare_pairs:
            if group_a in ct.index and group_b in ct.index:
                self.progress(80, f"Comparing {group_a} vs {group_b}...")
                mask = adata.obs[batch_key].isin([group_a, group_b])
                adata_sub = adata[mask]
                ct_sub_abs = pd.crosstab(adata_sub.obs[batch_key], adata_sub.obs[groupby])
                ct_sub = pd.crosstab(adata_sub.obs[batch_key], adata_sub.obs[groupby], normalize='index')

                chi2_sub, pval_sub = _run_stat_test(ct_sub_abs, stat_test, n_permutations)

                fig_sub = go.Figure()
                for col in ct_sub.columns:
                    fig_sub.add_trace(go.Bar(name=str(col), x=ct_sub.index.tolist(), y=ct_sub[col].values))
                fig_sub.update_layout(**self.get_plotly_layout(
                    title=f'Cell Proportions: {group_a} vs {group_b} (p={pval_sub:.4f})',
                    barmode='stack',
                    xaxis_title=batch_key, yaxis_title='Proportion',
                    width=600, height=400,
                ))
                safe_name = f'{group_a}_vs_{group_b}'.replace(' ', '_')
                result_files.append(self.save_plotly_json(fig_sub, plots_dir, f'proportion_compare_{safe_name}.json', 'bar', f'{group_a} vs {group_b} 比例比较'))

                ct_sub_abs.to_csv(os.path.join(results_dir, f'cell_counts_{safe_name}.csv'))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'proportion')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'chi2': round(float(chi2), 2),
                'p_value': float(pval),
                'n_batches': len(ct.index),
                'n_groups': len(ct.columns),
                'stat_test': stat_test,
            }
        }
