"""Regression tests for direct Bulk module scale handling."""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def _write_fpkm_table(tmp_path):
    path = tmp_path / 'example.fpkm.tsv'
    path.write_text(
        'gene\tctrl_1\tctrl_2\ttreat_1\ttreat_2\n'
        'G1\t1.2\t1.5\t8.1\t7.9\n'
        'G2\t2.0\t2.1\t1.8\t1.9\n'
        'G3\t0.1\t0.2\t0.3\t0.1\n',
        encoding='utf-8',
    )
    return str(path)


def test_bulk_deg_rejects_fpkm_for_count_model(tmp_path):
    from modules.bulk_deg import BulkDEGAnalysis
    path = _write_fpkm_table(tmp_path)
    analysis = BulkDEGAnalysis(str(tmp_path), {'method': 'deseq2'}, lambda *_: None)
    with pytest.raises(ValueError, match='原始整数 counts'):
        analysis.run(path)


def test_bulk_qc_keeps_continuous_samples_out_of_count_threshold_filter(tmp_path):
    from modules.bulk_qc import BulkQCAnalysis
    path = _write_fpkm_table(tmp_path)
    result = BulkQCAnalysis(str(tmp_path), {}, lambda *_: None).run(path)
    assert result['summary']['input_measurement'] == 'continuous_expression'
    assert result['summary']['samples_after'] == 4


def test_bulk_normalize_rejects_generic_raw_count_table(tmp_path):
    """A legacy CSV/TSV must not bypass count/design validation."""
    from modules.bulk_normalize import BulkNormalizeAnalysis

    path = tmp_path / 'raw_counts.tsv'
    path.write_text(
        'gene\tctrl_1\tctrl_2\ttreat_1\ttreat_2\n'
        'G1\t10\t12\t20\t22\n'
        'G2\t1\t2\t4\t5\n',
        encoding='utf-8',
    )

    with pytest.raises(ValueError, match='通用上传的 raw count 表'):
        BulkNormalizeAnalysis(
            str(tmp_path), {'method': 'deseq2', 'min_expr_samples': 0}, lambda *_: None,
        ).run(str(path))


def test_bulk_qc_rejects_invalid_h5ad_declared_as_raw_counts(tmp_path):
    """A user declaration cannot relabel fractional data as DESeq2 counts."""
    from modules.bulk_qc import BulkQCAnalysis

    path = tmp_path / 'not_counts.h5ad'
    adata = ad.AnnData(
        X=np.array([[1.2, 2.0], [3.0, 4.0]]),
        obs=pd.DataFrame(index=['S1', 'S2']),
        var=pd.DataFrame(index=['G1', 'G2']),
    )
    adata.uns['input_measurement'] = 'raw_counts'
    adata.write_h5ad(path)

    with pytest.raises(ValueError, match='非整数值'):
        BulkQCAnalysis(str(tmp_path), {}, lambda *_: None).run(str(path))
