import sys
import types

import numpy as np
import pandas as pd


class _FakeAdata:
    def __init__(self):
        self.obs = pd.DataFrame({'batch': ['batch_a', 'batch_a', 'batch_b', 'batch_b']})
        self.obsm = {'X_pca': np.zeros((4, 2))}


def test_bbknn_neighbors_auto_downshifts_to_min_batch(monkeypatch, tmp_path):
    from modules.batch_correct import BatchCorrectAnalysis

    calls = []

    def fake_bbknn(adata, **kwargs):
        calls.append(kwargs)

    monkeypatch.setitem(sys.modules, 'bbknn', types.SimpleNamespace(bbknn=fake_bbknn))

    mod = BatchCorrectAnalysis(
        project_dir=str(tmp_path),
        params={'bbknn_neighbors_within_batch': 3},
        progress_callback=lambda pct, msg: None,
    )
    corrected_key, method_info = mod._run_bbknn(_FakeAdata(), 'batch')

    assert corrected_key == 'X_pca'
    assert calls[0]['neighbors_within_batch'] == 2
    assert method_info['requested_neighbors_within_batch'] == 3
    assert method_info['neighbors_within_batch'] == 2
    assert method_info['min_batch_size'] == 2
    assert method_info['batch_counts'] == {'batch_a': 2, 'batch_b': 2}
    assert 'warning' in method_info
