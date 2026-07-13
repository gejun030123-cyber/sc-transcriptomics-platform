"""Regression tests for direct Bulk module scale handling."""
import os
import sys

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
