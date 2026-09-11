"""Regression tests for auditable Bulk DEG threshold diagnostics."""

import pandas as pd
import pytest


def test_comparison_filter_diagnostics_separates_fdr_and_fc_failures():
    from modules.bulk_deg import _comparison_filter_diagnostics

    deg_df = pd.DataFrame({
        'pvalue': [0.01, 0.001, 0.20, 0.03],
        'padj': [0.20, 0.04, 0.04, 0.80],
        'log2FC': [0.4, 1.2, -0.9, -1.1],
    })

    diagnostics = _comparison_filter_diagnostics(
        deg_df, fc_threshold=2.0, padj_threshold=0.05,
    )

    assert diagnostics['n_tested_genes'] == 4
    assert diagnostics['n_raw_p_lt_0_05'] == 3
    assert diagnostics['n_padj_lt_threshold'] == 2
    assert diagnostics['n_abs_log2fc_ge_threshold'] == 2
    assert diagnostics['n_both_padj_and_fc'] == 1
    assert diagnostics['effective_log2fc_threshold'] == pytest.approx(1.0)
    assert diagnostics['significance_metric'] == 'padj'
    assert diagnostics['multiple_testing_scope'] == 'per_comparison'
    assert diagnostics['min_pvalue'] == pytest.approx(0.001)
    assert diagnostics['min_padj'] == pytest.approx(0.04)


def test_bulk_deg_summary_records_per_comparison_diagnostics_and_effective_cutoff(tmp_path):
    import anndata
    import numpy as np
    from modules.bulk_deg import BulkDEGAnalysis

    rng = np.random.default_rng(9)
    counts = rng.poisson(80, size=(6, 16)).astype(float)
    counts[3:, :3] += 150
    input_path = tmp_path / 'counts.h5ad'
    anndata.AnnData(
        counts,
        obs=pd.DataFrame(
            {'condition': ['Ctrl'] * 3 + ['Treat'] * 3},
            index=[f'S{i}' for i in range(6)],
        ),
        var=pd.DataFrame(index=[f'G{i}' for i in range(16)]),
    ).write_h5ad(input_path)

    result = BulkDEGAnalysis(
        str(tmp_path), {
            'method': 't-test', 'groupby': 'condition',
            'comparisons': 'Treat-vs-Ctrl', 'fc_threshold': 2.0,
            'pval_threshold': 0.05,
        }, lambda *_: None,
    ).run(str(input_path))

    summary = result['summary']
    comparison = summary['per_comparison']['Treat-vs-Ctrl']
    diagnostics = comparison['filter_diagnostics']
    assert summary['effective_log2fc_threshold'] == pytest.approx(1.0)
    assert summary['significance_metric'] == 'padj'
    assert summary['multiple_testing_scope'] == 'per_comparison'
    assert diagnostics['n_tested_genes'] == 16
    assert diagnostics['n_both_padj_and_fc'] == comparison['n_up'] + comparison['n_down']

    complete_csv = next(
        item['file_path'] for item in result['result_files']
        if item['label'].startswith('差异表达基因列表')
    )
    columns = pd.read_csv(complete_csv).columns
    assert {'passes_padj', 'passes_fc', 'significant'} <= set(columns)
