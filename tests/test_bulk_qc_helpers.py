# tests/test_bulk_qc_helpers.py
import os
import numpy as np
import pytest


def test_gini_uniform():
    """完全均匀分布 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([100.0, 100.0, 100.0, 100.0])
    assert abs(_gini(vals)) < 1e-6


def test_gini_concentrated():
    """一个值占主导 → Gini 接近 1"""
    from modules.bulk_qc import _gini
    vals = np.array([0.0] * 20 + [1000.0])
    assert _gini(vals) > 0.9


def test_gini_all_zeros():
    """全零输入 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([0.0, 0.0, 0.0])
    assert _gini(vals) == 0.0


def test_gini_empty():
    """空数组 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([])
    assert _gini(vals) == 0.0


def test_infer_groups_dash():
    """样本名含 '-' → 取前缀"""
    from modules.bulk_qc import _infer_groups
    names = ['ctrl-1', 'ctrl-2', 'hmc3-1', 'hmc3-2']
    assert _infer_groups(names) == ['ctrl', 'ctrl', 'hmc3', 'hmc3']


def test_infer_groups_underscore():
    """样本名含 '_' → 取前缀"""
    from modules.bulk_qc import _infer_groups
    names = ['ctrl_1', 'rapa_3']
    assert _infer_groups(names) == ['ctrl', 'rapa']


def test_infer_groups_no_sep():
    """无分隔符 → 样本名本身"""
    from modules.bulk_qc import _infer_groups
    names = ['sampleA', 'sampleB']
    assert _infer_groups(names) == ['sampleA', 'sampleB']


def test_infer_multifactor_groups_preserves_combined_group():
    """Ctr_B_1 应默认分到 Ctr_B，不应与 Ctr_En 合并。"""
    from modules.bulk_qc import _infer_groups

    names = ['Ctr_B_1', 'Ctr_B_2', 'Ctr_En_1', 'Ctr_En_2']
    assert _infer_groups(names) == ['Ctr_B', 'Ctr_B', 'Ctr_En', 'Ctr_En']


def test_multifactor_group_candidates_include_combined_and_factors():
    from modules.io_utils import infer_sample_group_candidates

    names = [
        'Ctr_B_1', 'Ctr_B_2', 'Ctr_En_1', 'Ctr_En_2',
        'PEA_B_1', 'PEA_B_2', 'PEA_En_1', 'PEA_En_2',
    ]
    candidates = infer_sample_group_candidates(names)

    assert [c['key'] for c in candidates] == ['combined', 'factor_1', 'factor_2']
    assert candidates[0]['values'] == ['Ctr_B', 'Ctr_En', 'PEA_B', 'PEA_En']
    assert candidates[0]['group_sizes'] == {
        'Ctr_B': 2, 'Ctr_En': 2, 'PEA_B': 2, 'PEA_En': 2,
    }
    assert candidates[1]['values'] == ['Ctr', 'PEA']
    assert candidates[2]['values'] == ['B', 'En']


def test_simple_group_names_do_not_duplicate_candidate():
    from modules.io_utils import infer_sample_group_candidates

    candidates = infer_sample_group_candidates(['Ctrl_1', 'Ctrl_2', 'Drug_1', 'Drug_2'])
    assert len(candidates) == 1
    assert candidates[0]['values'] == ['Ctrl', 'Drug']


def test_unstructured_sample_names_have_no_auto_group():
    from modules.io_utils import infer_sample_group_candidates

    assert infer_sample_group_candidates(['sampleA', 'sampleB']) == []


def test_obs_columns_api_returns_multifactor_candidates(test_project):
    from app import create_app
    from config import Config

    path = os.path.join(Config.uploads_dir(test_project), 'bulk.tsv')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(
            'gene\tCtr_B_1\tCtr_B_2\tCtr_En_1\tCtr_En_2\tPEA_B_1\tPEA_B_2\tPEA_En_1\tPEA_En_2\n'
            'G1\t1\t2\t3\t4\t5\t6\t7\t8\n'
            'G2\t8\t7\t6\t5\t4\t3\t2\t1\n'
        )

    app = create_app()
    app.config['TESTING'] = True
    response = app.test_client().get('/api/obs-columns', query_string={'file_path': path})
    payload = response.get_json()

    assert response.status_code == 200
    assert payload['sample_groups']['auto_group']['values'] == ['Ctr_B', 'Ctr_En', 'PEA_B', 'PEA_En']
    assert len(payload['sample_groups']['auto_group_candidates']) == 3


def test_obs_columns_api_separates_group_batch_and_cluster_suggestions(test_project):
    import anndata
    import pandas as pd
    from app import create_app
    from config import Config

    n_obs = 24
    obs = pd.DataFrame({
        # Deliberately put technical columns first: this reproduces the old
        # physical-order bug in the web form.
        'barcode': [f'cell_{i}' for i in range(n_obs)],
        'batch': pd.Categorical(['B1'] * 12 + ['B2'] * 12),
        'sample': pd.Categorical([f'S{i // 6 + 1}' for i in range(n_obs)]),
        'condition': pd.Categorical(['Control'] * 12 + ['Treatment'] * 12),
        'leiden': pd.Categorical(['0'] * 12 + ['1'] * 12),
    }, index=[f'cell_{i}' for i in range(n_obs)])
    adata = anndata.AnnData(
        np.ones((n_obs, 4)),
        obs=obs,
        var=pd.DataFrame(index=['G1', 'G2', 'G3', 'G4']),
    )
    path = os.path.join(Config.uploads_dir(test_project), 'grouping.h5ad')
    adata.write_h5ad(path)

    app = create_app()
    app.config['TESTING'] = True
    response = app.test_client().get(
        '/api/obs-columns',
        query_string={'file_path': path, 'module_name': 'deg'},
    )
    payload = response.get_json()

    assert response.status_code == 200
    assert payload['suggestions']['groupby'] == 'condition'
    assert payload['suggestions']['color_by'] == 'condition'
    assert payload['suggestions']['batch_key'] == 'batch'
    assert payload['suggestions']['cluster_key'] == 'leiden'
    assert payload['suggestions']['sample_key'] == 'sample'
    assert payload['grouping_candidates'][0] == 'condition'


def test_detect_outliers_no_outlier():
    """正常聚类数据 → 无离群点"""
    from modules.bulk_qc import _detect_outliers_mahal
    np.random.seed(42)
    coords = np.random.randn(20, 3)
    outliers = _detect_outliers_mahal(coords, [f's{i}' for i in range(20)])
    assert len(outliers) == 0


def test_detect_outliers_with_outlier():
    """含明显离群点 → 检测到"""
    from modules.bulk_qc import _detect_outliers_mahal
    np.random.seed(42)
    coords = np.random.randn(20, 3)
    coords = np.vstack([coords, [[50, 50, 50]]])  # 离群点
    samples = [f's{i}' for i in range(21)]
    outliers = _detect_outliers_mahal(coords, samples)
    assert 's20' in outliers


def test_detect_outliers_too_few():
    """样本数 < 4 → 返回空"""
    from modules.bulk_qc import _detect_outliers_mahal
    coords = np.array([[0, 0], [1, 1], [2, 2]])
    outliers = _detect_outliers_mahal(coords, ['a', 'b', 'c'])
    assert outliers == []
