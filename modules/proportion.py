import os
import json
import re
import numpy as np
from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping

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
        valid_permutations = 0
        import pandas as pd
        # Conditional label permutation keeps both row totals (group cell
        # counts) and column totals (cell-type totals) fixed.  Independent
        # multinomial draws used to generate a different null table with
        # random margins, so the reported p value was not a permutation test
        # for the observed contingency table.
        rng = np.random.default_rng(0)
        row_sums = ct_abs.sum(axis=1).to_numpy(dtype=int)
        col_sums = ct_abs.sum(axis=0).to_numpy(dtype=int)
        labels = np.repeat(np.arange(len(col_sums), dtype=int), col_sums)
        for _ in range(max(1, int(n_permutations))):
            shuffled_labels = rng.permutation(labels)
            perm_table = np.zeros_like(ct_abs.values, dtype=int)
            offset = 0
            for row_index, row_size in enumerate(row_sums):
                assigned = shuffled_labels[offset:offset + row_size]
                perm_table[row_index, :] = np.bincount(
                    assigned, minlength=len(col_sums),
                )
                offset += row_size
            perm_df = pd.DataFrame(perm_table, index=ct_abs.index, columns=ct_abs.columns)
            try:
                perm_chi2, _, _, _ = chi2_contingency(perm_df)
            except ValueError:
                continue
            valid_permutations += 1
            if perm_chi2 >= observed_chi2:
                count += 1
        pval = ((count + 1) / (valid_permutations + 1)
                if valid_permutations else float('nan'))
        return observed_chi2, pval
    else:
        chi2, pval, _, _ = chi2_contingency(ct_abs)
        return chi2, pval


_TECHNICAL_BATCH_NAMES = {
    'batch', 'technical_batch', 'sequencing_batch', 'library_batch',
    'lane', 'run', 'sequencing_run',
}


def _bh_adjust(pvalues):
    """BH-adjust finite p values without adding a heavy dependency."""
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted
    indices = np.where(valid)[0]
    ordered = indices[np.argsort(values[indices])]
    running = 1.0
    for rank in range(len(ordered), 0, -1):
        index = ordered[rank - 1]
        running = min(running, values[index] * len(ordered) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _sample_level_composition(adata, groupby, sample_key, condition_key,
                              min_samples_per_condition=2,
                              min_cells_per_sample=10,
                              allow_technical_sample=False,
                              allow_technical_condition=False,
                              condition_pairs=None):
    """Build sample × celltype proportions and safe per-celltype tests.

    Cell counts are useful descriptive evidence but are not biological
    replicates.  This helper deliberately requires a valid biological sample
    column before emitting p values/FDR.
    """
    import pandas as pd
    from scipy.stats import kruskal, mannwhitneyu

    warnings = []
    base = {
        'available': False,
        'inference_ready': False,
        'warnings': warnings,
        'proportions': pd.DataFrame(),
        'tests': pd.DataFrame(),
        'design': pd.DataFrame(),
        'condition_counts': {},
    }
    sample_key = str(sample_key or '').strip()
    condition_key = str(condition_key or '').strip()
    if not sample_key:
        warnings.append('未指定生物学样本列；不报告样本级 p 值/FDR。')
        return base
    if not condition_key:
        warnings.append('未指定条件列；不报告样本级组间 p 值/FDR。')
        return base
    if sample_key not in adata.obs.columns or condition_key not in adata.obs.columns:
        warnings.append('样本列或条件列不存在；不报告样本级 p 值/FDR。')
        return base
    if sample_key == condition_key:
        warnings.append('样本列与条件列不能相同；不报告样本级 p 值/FDR。')
        return base
    if sample_key.lower() in _TECHNICAL_BATCH_NAMES and not allow_technical_sample:
        warnings.append(f"样本列 '{sample_key}' 看起来是技术 batch；未确认其代表生物学样本，不报告 p 值/FDR。")
        return base
    if condition_key.lower() in _TECHNICAL_BATCH_NAMES and not allow_technical_condition:
        warnings.append(f"条件列 '{condition_key}' 看起来是技术 batch；未确认其代表生物学条件，不报告 p 值/FDR。")
        return base

    obs = adata.obs
    sample_values = obs[sample_key]
    condition_values = obs[condition_key]
    missing = (
        sample_values.isna() | condition_values.isna()
        | sample_values.astype(str).str.strip().eq('')
        | condition_values.astype(str).str.strip().eq('')
    )
    if bool(missing.any()):
        warnings.append('样本列或条件列有缺失值；不报告样本级 p 值/FDR。')
        return base

    work = pd.DataFrame({
        'sample': sample_values.astype(str),
        'condition': condition_values.astype(str),
        'celltype': obs[groupby].astype(str),
    })
    sample_condition_count = work[['sample', 'condition']].drop_duplicates().groupby(
        'sample', observed=True
    )['condition'].nunique()
    if bool((sample_condition_count > 1).any()):
        warnings.append('同一个样本属于多个条件；当前模块不支持配对/交互模型，已关闭 p 值/FDR。')
        return base

    sample_design = work.groupby(['sample', 'condition'], observed=True).size().rename(
        'total_cells'
    ).reset_index()
    sample_design['eligible_for_inference'] = sample_design['total_cells'] >= int(min_cells_per_sample)
    celltypes = sorted(work['celltype'].unique().tolist())
    if not celltypes:
        warnings.append('没有可用的细胞类型，无法构建样本级组成表。')
        return base

    counts = work.groupby(['sample', 'condition', 'celltype'], observed=True).size()
    full_index = pd.MultiIndex.from_product(
        [sample_design['sample'].tolist(), celltypes], names=['sample', 'celltype']
    )
    # Each sample maps to exactly one condition above, so construct the full
    # table from (sample, celltype) and merge the design rather than dropping
    # zero-cell types.
    by_sample_celltype = counts.groupby(level=['sample', 'celltype']).sum().reindex(full_index, fill_value=0)
    proportions = by_sample_celltype.rename('n_cells').reset_index().merge(
        sample_design, on='sample', how='left', validate='many_to_one'
    )
    proportions['proportion'] = proportions['n_cells'] / proportions['total_cells'].clip(lower=1)
    proportions = proportions[[
        'sample', 'condition', 'celltype', 'n_cells', 'total_cells',
        'eligible_for_inference', 'proportion',
    ]]

    eligible = proportions[proportions['eligible_for_inference']].copy()
    condition_counts = eligible[['sample', 'condition']].drop_duplicates().groupby(
        'condition', observed=True
    )['sample'].nunique()
    condition_counts = {str(key): int(value) for key, value in condition_counts.items()}
    base.update({
        'available': True,
        'proportions': proportions,
        'design': sample_design,
        'condition_counts': condition_counts,
    })
    if len(condition_counts) < 2:
        warnings.append('有效样本只覆盖一个条件；仅输出样本级描述统计。')
        return base
    undersized = {key: value for key, value in condition_counts.items()
                  if value < int(min_samples_per_condition)}
    if undersized:
        detail = ', '.join(f'{key}=n{value}' for key, value in undersized.items())
        warnings.append(f'条件独立样本数不足 {int(min_samples_per_condition)}（{detail}）；仅输出描述统计。')
        return base

    conditions = sorted(condition_counts)
    rows = []
    if condition_pairs:
        # 指定比较组时，样本级检验只运行列出的对（与细胞级 pairwise 一致），
        # 避免“指定了 A-vs-B 却仍输出全部条件 Kruskal”的口径漂移。
        for experimental, control in condition_pairs:
            if experimental not in condition_counts or control not in condition_counts:
                continue
            for celltype in celltypes:
                values_exp = eligible.loc[
                    (eligible['celltype'] == celltype) & (eligible['condition'] == experimental),
                    'proportion',
                ].to_numpy(dtype=float)
                values_ctrl = eligible.loc[
                    (eligible['celltype'] == celltype) & (eligible['condition'] == control),
                    'proportion',
                ].to_numpy(dtype=float)
                try:
                    statistic, pvalue = mannwhitneyu(
                        values_exp, values_ctrl, alternative='two-sided'
                    )
                    test_name = 'mann_whitney_u'
                except Exception:
                    statistic, pvalue, test_name = np.nan, np.nan, 'unavailable'
                row = {
                    'celltype': str(celltype),
                    "comparison": f"{experimental} vs {control}",
                    'test': test_name,
                    'statistic': float(statistic) if np.isfinite(statistic) else np.nan,
                    'p_value': float(pvalue) if np.isfinite(pvalue) else np.nan,
                    'n_samples': int(len(values_exp) + len(values_ctrl)),
                }
                row[f'n_{experimental}'] = int(len(values_exp))
                row[f'n_{control}'] = int(len(values_ctrl))
                row[f'mean_{experimental}'] = float(np.mean(values_exp)) if len(values_exp) else np.nan
                row[f'mean_{control}'] = float(np.mean(values_ctrl)) if len(values_ctrl) else np.nan
                rows.append(row)
    else:
        for celltype in celltypes:
            values_by_condition = [
                eligible.loc[
                    (eligible['celltype'] == celltype) & (eligible['condition'] == condition),
                    'proportion',
                ].to_numpy(dtype=float)
                for condition in conditions
            ]
            try:
                if len(conditions) == 2:
                    statistic, pvalue = mannwhitneyu(
                        values_by_condition[0], values_by_condition[1], alternative='two-sided'
                    )
                    test_name = 'mann_whitney_u'
                else:
                    statistic, pvalue = kruskal(*values_by_condition)
                    test_name = 'kruskal'
            except Exception:
                statistic, pvalue, test_name = np.nan, np.nan, 'unavailable'
            row = {
                'celltype': str(celltype), 'test': test_name,
                'statistic': float(statistic) if np.isfinite(statistic) else np.nan,
                'p_value': float(pvalue) if np.isfinite(pvalue) else np.nan,
                'n_samples': int(sum(len(values) for values in values_by_condition)),
            }
            for condition, values in zip(conditions, values_by_condition):
                row[f'n_{condition}'] = int(len(values))
                row[f'mean_{condition}'] = float(np.mean(values)) if len(values) else np.nan
            rows.append(row)
    tests = pd.DataFrame(rows)
    if not tests.empty:
        tests['fdr_bh'] = _bh_adjust(tests['p_value'].to_numpy(dtype=float))
    base.update({'inference_ready': True, 'tests': tests})
    return base

class ProportionAnalysis(BaseAnalysis):
    MODULE_NAME = "proportion"
    DISPLAY_NAME = "细胞比例分析"
    DESCRIPTION = "各分组间的细胞比例差异分析"
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        import scanpy as sc
        import pandas as pd
        from scipy.stats import chi2_contingency
        import matplotlib.pyplot as plt
        from modules.native_figures import grouped_bar_figure, heatmap_figure
        from modules.figure_style import NATURE_PALETTE

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata)
        requested_groupby = str(self.params.get('groupby', 'celltype') or '').strip()
        requested_batch_key = str(self.params.get('batch_key', 'batch') or '').strip()
        requested_condition_key = str(self.params.get('condition_key', '') or '').strip()
        requested_sample_key = str(self.params.get('sample_key', '') or '').strip()
        analysis_unit = str(self.params.get('analysis_unit', 'auto') or 'auto').strip().lower()
        stat_test = self.params.get('stat_test', 'chi_square')
        n_permutations = int(self.params.get('n_permutations', 1000))
        min_cells_per_group = int(self.params.get('min_cells_per_group', 10))
        min_cells_per_sample = int(self.params.get('min_cells_per_sample', 10))
        min_samples_per_condition = int(self.params.get('min_samples_per_condition', 2))

        # 指定比较组（A-vs-B）同时约束细胞级与样本级检验，避免口径不一致。
        compare_groups_str = self.params.get('compare_groups', '').strip()
        compare_pairs = []
        if compare_groups_str:
            for item in compare_groups_str.replace('\n', ';').split(';'):
                item = item.strip()
                if '-vs-' in item:
                    parts = item.split('-vs-')
                    if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                        compare_pairs.append((parts[0].strip(), parts[1].strip()))

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
        # ``batch_key`` remains the compatible display/comparison fallback.
        # A supplied condition_key takes precedence for the new sample-level
        # analysis so technical batch is not silently treated as biology.
        requested_comparison_key = requested_condition_key or requested_batch_key
        comparison_key, comparison_info = resolve_obs_grouping(
            adata, requested_comparison_key, max_categories=50,
            max_numeric_categories=20, require_multiple=True,
        )
        if comparison_key is None:
            comparison_label = 'condition_key' if requested_condition_key else 'batch_key'
            if requested_comparison_key not in adata.obs.columns:
                raise ValueError(f"{comparison_label} '{requested_comparison_key}' 不在 adata.obs 中")
            raise ValueError(
                f"{comparison_label} '{requested_comparison_key}' 不是有效的分类分组列："
                f"{comparison_info.get('reason', '')}"
            )
        if groupby != requested_groupby:
            self.progress(-1, f"分组列已改用 '{groupby}'：{group_info.get('requested_reason', '')}")
        if comparison_key != requested_comparison_key:
            self.progress(-1, f"比较列已改用 '{comparison_key}'：{comparison_info.get('requested_reason', '')}")

        # Compute sample-level composition before removing rare cell types for
        # the cell-level table.  Filtering globally first changes every
        # sample's denominator and can manufacture condition differences when
        # a cell type is rare in only one condition.
        if analysis_unit == 'cell':
            sample_level = {
                'available': False, 'inference_ready': False,
                'warnings': ['已选择细胞级描述模式；不报告样本级 p 值/FDR。'],
                'proportions': pd.DataFrame(), 'tests': pd.DataFrame(),
                'design': pd.DataFrame(), 'condition_counts': {},
            }
        else:
            sample_level = _sample_level_composition(
                adata, groupby, requested_sample_key, requested_condition_key,
                min_samples_per_condition=min_samples_per_condition,
                min_cells_per_sample=min_cells_per_sample,
                allow_technical_sample=bool(self.params.get('confirm_batch_is_biological_sample', False)),
                allow_technical_condition=bool(self.params.get('confirm_batch_is_biological_condition', False)),
                condition_pairs=compare_pairs,
            )

        self.progress(30, "Computing cell proportions...")
        removed_groups = []
        if min_cells_per_group > 0:
            group_counts = adata.obs[groupby].value_counts()
            valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()
            if len(valid_groups) < len(group_counts):
                removed = set(group_counts.index) - set(valid_groups)
                removed_groups = sorted(map(str, removed))
                self.progress(-1, f"移除 {len(removed)} 个低细胞数组: {removed}")
                adata = adata[adata.obs[groupby].isin(valid_groups)].copy()
        ct = pd.crosstab(adata.obs[comparison_key], adata.obs[groupby], normalize='index')
        ct_abs = pd.crosstab(adata.obs[comparison_key], adata.obs[groupby])
        if ct.shape[0] < 2 or ct.shape[1] == 0:
            raise ValueError(
                f'过滤后没有可比较的细胞类型/分组（列联表 {ct.shape}）。'
                '请放宽 min_cells_per_group 或检查分组列。'
            )

        # 过滤后列联表可能为空或退化为单组：chi2_contingency 在空表上抛
        # ValueError、单行/单列表返回误导性 p=1.0。不足两组时跳过检验并在
        # summary 中明确记录，而不是崩溃或输出无意义 p 值。
        stat_skipped = False
        if ct_abs.shape[0] < 2 or ct_abs.shape[1] < 2:
            stat_skipped = True
            self.progress(-1, '过滤后有效分组不足两组，跳过列联表检验。')
            chi2, pval = float('nan'), float('nan')
        else:
            self.progress(50, "Running chi-squared test...")
            chi2, pval = _run_stat_test(ct_abs, stat_test, n_permutations)

        self.progress(65, "Generating proportion plots...")
        plots_dir = self.ensure_plots_dir()
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        fig = grouped_bar_figure(
            ct.index.tolist(), [(str(col), ct[col].values) for col in ct.columns],
            title='Cell Proportions by Group', x_label=comparison_key,
            y_label='Proportion', rotation=35, stacked=True,
        )
        result_files.extend(self.save_matplotlib_figure(
            fig, plots_dir, 'proportion_stacked.png', 'bar',
            'Cell Proportions (Stacked)', formats=('png', 'svg'), dpi=300,
        ))

        if self.params.get('show_proportion_heatmap', True):
            fig_heat = heatmap_figure(
                ct.values, x_labels=[str(x) for x in ct.columns],
                y_labels=[str(x) for x in ct.index], title='Cell Proportion Heatmap',
                x_label=groupby, y_label=comparison_key, colorbar_label='Proportion',
                vmin=0, vmax=max(1.0, float(ct.values.max()) if ct.size else 1.0),
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_heat, plots_dir, 'proportion_heatmap.png', 'heatmap',
                'Cell Proportion Heatmap', formats=('png', 'svg'), dpi=300,
            ))

        n_batches = len(ct.index)
        fig2, pie_axes = plt.subplots(1, max(1, n_batches),
                                      figsize=(5.0 * max(1, n_batches), 4.8), dpi=150,
                                      squeeze=False)
        for i, idx in enumerate(ct.index):
            ax = pie_axes[0, i]
            ax.pie(ct.loc[idx].values,
                   labels=[str(x) for x in ct.columns],
                   colors=[NATURE_PALETTE[j % len(NATURE_PALETTE)] for j in range(len(ct.columns))],
                   startangle=90, counterclock=False,
                   autopct='%1.1f%%', pctdistance=0.78,
                   wedgeprops={'width': 0.42, 'edgecolor': 'white', 'linewidth': 0.7},
                   textprops={'fontsize': 7})
            ax.set_title(str(idx), loc='left', fontsize=9)
        fig2.suptitle(f'Cell Type Distribution per {comparison_key}', x=0.05, ha='left',
                      fontsize=11, fontweight='semibold')
        result_files.extend(self.save_matplotlib_figure(
            fig2, plots_dir, 'proportion_pie.png', 'pie',
            'Cell Type Distribution', formats=('png', 'svg'), dpi=300,
        ))

        ct_abs.to_csv(os.path.join(results_dir, 'cell_counts.csv'))
        ct.to_csv(os.path.join(results_dir, 'cell_proportions.csv'))
        result_files.append({'file_path': os.path.join(results_dir, 'cell_counts.csv'), 'file_type': 'csv', 'category': 'table', 'label': 'Cell Counts'})
        result_files.append({'file_path': os.path.join(results_dir, 'cell_proportions.csv'), 'file_type': 'csv', 'category': 'table', 'label': 'Cell Proportions'})

        if sample_level['available']:
            sample_prop_path = os.path.join(results_dir, 'sample_level_cell_proportions.csv')
            sample_tests_path = os.path.join(results_dir, 'sample_level_proportion_tests.csv')
            sample_design_path = os.path.join(results_dir, 'sample_level_design_summary.csv')
            sample_level['proportions'].to_csv(sample_prop_path, index=False)
            sample_level['tests'].to_csv(sample_tests_path, index=False)
            sample_level['design'].to_csv(sample_design_path, index=False)
            result_files.extend([
                {'file_path': sample_prop_path, 'file_type': 'csv', 'category': 'table',
                 'label': 'Sample-level Cell Proportions'},
                {'file_path': sample_tests_path, 'file_type': 'csv', 'category': 'table',
                 'label': 'Sample-level Proportion Tests (BH-FDR)'},
                {'file_path': sample_design_path, 'file_type': 'csv', 'category': 'table',
                 'label': 'Sample-level Design Summary'},
            ])
            sample_means = sample_level['proportions'].groupby(
                ['condition', 'celltype'], observed=True
            )['proportion'].mean().unstack(fill_value=0)
            if not sample_means.empty:
                fig_sample = grouped_bar_figure(
                    sample_means.index.tolist(),
                    [(str(column), sample_means[column].values) for column in sample_means.columns],
                    title='Mean Cell Proportions across Biological Samples',
                    x_label=requested_condition_key, y_label='Mean sample proportion',
                    rotation=35, stacked=True,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_sample, plots_dir, 'proportion_sample_level.png', 'bar',
                    'Sample-level Mean Cell Proportions', formats=('png', 'svg'), dpi=300,
                ))

        # Pairwise group comparison（compare_pairs 已在开头统一解析）

        for group_a, group_b in compare_pairs:
            if group_a in ct.index and group_b in ct.index:
                self.progress(80, f"Comparing {group_a} vs {group_b}...")
                mask = adata.obs[comparison_key].isin([group_a, group_b])
                adata_sub = adata[mask]
                ct_sub_abs = pd.crosstab(adata_sub.obs[comparison_key], adata_sub.obs[groupby])
                ct_sub = pd.crosstab(adata_sub.obs[comparison_key], adata_sub.obs[groupby], normalize='index')

                chi2_sub, pval_sub = _run_stat_test(ct_sub_abs, stat_test, n_permutations) \
                    if (ct_sub_abs.shape[0] >= 2 and ct_sub_abs.shape[1] >= 2) \
                    else (float('nan'), float('nan'))

                fig_sub = grouped_bar_figure(
                    ct_sub.index.tolist(), [(str(col), ct_sub[col].values) for col in ct_sub.columns],
                    title=f'Cell Proportions: {group_a} vs {group_b} (p={pval_sub:.4f})',
                    x_label=comparison_key, y_label='Proportion', rotation=35, stacked=True,
                )
                safe_name = re.sub(r'[^A-Za-z0-9_.-]', '_', f'{group_a}_vs_{group_b}')
                result_files.extend(self.save_matplotlib_figure(
                    fig_sub, plots_dir, f'proportion_compare_{safe_name}.png', 'bar',
                    f'{group_a} vs {group_b} 比例比较', formats=('png', 'svg'), dpi=300,
                ))

                ct_sub_path = os.path.join(results_dir, f'cell_counts_{safe_name}.csv')
                ct_sub_abs.to_csv(ct_sub_path)
                result_files.append({'file_path': ct_sub_path, 'file_type': 'csv', 'category': 'table', 'label': f'{group_a} vs {group_b} 细胞计数'})

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'proportion')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                # 语义化统计量名：fisher_exact 返回的是 OR 而不是 chi2。
                'statistic': round(float(chi2), 2) if not stat_skipped else None,
                'chi2': (
                    round(float(chi2), 2)
                    if not stat_skipped and stat_test in ('chi_square', 'permutation')
                    else None
                ),
                'p_value': float(pval) if not stat_skipped else None,
                'n_batches': len(ct.index),
                'n_groups': len(ct.columns),
                'stat_test': stat_test,
                'stat_skipped': bool(stat_skipped),
                'removed_low_cell_groups': removed_groups,
                'cell_level_association_p_value': float(pval) if not stat_skipped else None,
                'inference_unit': ('sample'
                                   if sample_level['inference_ready']
                                   else ('sample (downgraded to descriptive)' if requested_sample_key else 'cell (descriptive only)')),
                'sample_key': requested_sample_key,
                'condition_key': requested_condition_key,
                'sample_level_inference_ready': bool(sample_level['inference_ready']),
                'sample_level_condition_counts': sample_level['condition_counts'],
                'sample_level_n_tests': int(len(sample_level['tests'])),
                'sample_level_warnings': sample_level['warnings'],
                'requested_groupby': requested_groupby,
                'groupby': groupby,
                'scope_key': str(self.params.get('scope_key', '') or '').strip() or None,
                'scope_values': ([v.strip() for v in str(self.params.get('scope_values', '') or '').split(',') if v.strip()] or None),
                'requested_batch_key': requested_batch_key,
                'batch_key': comparison_key,
                'comparison_key': comparison_key,
            }
        }
