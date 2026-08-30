"""簇 Marker（deg）比较组选择与样本范围（scope）的契约测试。

覆盖本次修复的核心问题：
1. deg 现在可以选择“哪个细胞 vs 哪些细胞”（custom A-vs-B / pairwise）。
2. deg 现在可以限制在哪个样品/条件下比较（scope_key + scope_values）。
3. 旧参数（仅 reference，无 comparison_mode）保持向后兼容。
4. scope 校验错误路径（取值不存在、未选值）。
"""

import os

import numpy as np
import pandas as pd
import pytest


def _make_sc_anndata(n_obs=90, n_vars=200, n_clusters=3, n_samples=3):
    import anndata
    import scanpy as _sc

    np.random.seed(42)
    raw = np.random.poisson(lam=5, size=(n_obs, n_vars)).astype(np.float32)
    raw[raw == 0] = 1
    obs = pd.DataFrame(index=[f'cell_{i}' for i in range(n_obs)])
    obs['leiden'] = pd.Categorical([str(i % n_clusters) for i in range(n_obs)])
    sample_ids = []
    for sample_index in range(n_samples):
        sample_ids.extend([f'S{sample_index + 1}'] * (n_obs // n_samples))
    sample_ids = sample_ids[:n_obs]
    while len(sample_ids) < n_obs:
        sample_ids.append(sample_ids[-1])
    obs['sample_id'] = sample_ids
    var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_vars)])
    adata = anndata.AnnData(X=raw.copy(), obs=obs, var=var)
    adata.layers['counts'] = raw.copy()
    _sc.pp.normalize_total(adata, target_sum=1e4)
    _sc.pp.log1p(adata)
    _sc.pp.pca(adata, n_comps=10, svd_solver='arpack')
    _sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
    _sc.tl.umap(adata)
    return adata


def _run_deg(tmp_path, params):
    from modules.deg import DEGAnalysis

    adata = _make_sc_anndata()
    input_path = str(tmp_path / 'input.h5ad')
    adata.write_h5ad(input_path)
    module = DEGAnalysis(
        project_dir=str(tmp_path), params=params,
        progress_callback=lambda p, m: None,
    )
    return module.run(input_path)


class TestDegComparisonSelection:
    def test_custom_pairs_only_run_listed_comparisons(self, tmp_path):
        result = _run_deg(tmp_path, {
            'groupby': 'leiden', 'method': 'wilcoxon', 'n_genes': 10,
            'comparison_mode': 'custom', 'comparisons': '0-vs-2',
            'min_cells_per_group': 5,
        })
        summary = result['summary']
        assert summary['comparison_mode'] == 'custom'
        assert summary['comparisons'] == ['0 vs 2']
        df = pd.read_csv(os.path.join(str(tmp_path), 'results', 'deg_results.csv'))
        assert set(df['comparison'].astype(str).unique()) == {'0 vs 2'}
        assert len(df) == summary['total_deg_genes']

    def test_pairwise_mode_generates_all_pairs(self, tmp_path):
        result = _run_deg(tmp_path, {
            'groupby': 'leiden', 'method': 'wilcoxon', 'n_genes': 5,
            'comparison_mode': 'pairwise', 'min_cells_per_group': 5,
        })
        assert result['summary']['comparison_mode'] == 'pairwise'
        df = pd.read_csv(os.path.join(str(tmp_path), 'results', 'deg_results.csv'))
        assert sorted(df['comparison'].astype(str).unique()) == ['0 vs 1', '0 vs 2', '1 vs 2']

    def test_legacy_reference_behavior_preserved(self, tmp_path):
        # 不传 comparison_mode：保持旧的 reference/rest 行为。
        result = _run_deg(tmp_path, {
            'groupby': 'leiden', 'method': 'wilcoxon', 'n_genes': 5,
            'reference': 'rest',
        })
        summary = result['summary']
        assert summary['comparison_mode'] == 'reference'
        assert summary['reference'] == 'rest'
        df = pd.read_csv(os.path.join(str(tmp_path), 'results', 'deg_results.csv'))
        assert sorted(df['comparison'].astype(str).unique()) == [
            '0 vs rest', '1 vs rest', '2 vs rest',
        ]

    def test_reference_specific_group(self, tmp_path):
        result = _run_deg(tmp_path, {
            'groupby': 'leiden', 'method': 'wilcoxon', 'n_genes': 5,
            'reference': '1',
        })
        df = pd.read_csv(os.path.join(str(tmp_path), 'results', 'deg_results.csv'))
        assert sorted(df['comparison'].astype(str).unique()) == ['0 vs 1', '2 vs 1']

    def test_invalid_custom_pair_raises(self, tmp_path):
        from modules.deg import DEGAnalysis

        adata = _make_sc_anndata()
        input_path = str(tmp_path / 'input.h5ad')
        adata.write_h5ad(input_path)
        module = DEGAnalysis(
            project_dir=str(tmp_path),
            params={'groupby': 'leiden', 'comparison_mode': 'custom',
                    'comparisons': '0-vs-99'},
            progress_callback=lambda p, m: None,
        )
        with pytest.raises(ValueError, match='比较组不存在'):
            module.run(input_path)


class TestDegScope:
    def test_scope_restricts_cells_and_is_recorded(self, tmp_path):
        result = _run_deg(tmp_path, {
            'groupby': 'leiden', 'method': 'wilcoxon', 'n_genes': 5,
            'comparison_mode': 'custom', 'comparisons': '0-vs-2',
            'scope_key': 'sample_id', 'scope_values': 'S1,S2',
            'min_cells_per_group': 5,
        })
        summary = result['summary']
        assert summary['scope_key'] == 'sample_id'
        assert summary['scope_values'] == ['S1', 'S2']
        output = summary['output_adata'] if 'output_adata' in summary else result['output_adata']
        import anndata
        out = anndata.read_h5ad(output)
        assert set(out.obs['sample_id'].astype(str).unique()) == {'S1', 'S2'}

    def test_scope_missing_values_raises(self, tmp_path):
        from modules.deg import DEGAnalysis

        adata = _make_sc_anndata()
        input_path = str(tmp_path / 'input.h5ad')
        adata.write_h5ad(input_path)
        module = DEGAnalysis(
            project_dir=str(tmp_path),
            params={'groupby': 'leiden', 'comparison_mode': 'custom',
                    'comparisons': '0-vs-2',
                    'scope_key': 'sample_id', 'scope_values': 'S9,S10'},
            progress_callback=lambda p, m: None,
        )
        with pytest.raises(ValueError, match='不存在'):
            module.run(input_path)

    def test_scope_key_without_values_raises(self, tmp_path):
        from modules.deg import DEGAnalysis

        adata = _make_sc_anndata()
        input_path = str(tmp_path / 'input.h5ad')
        adata.write_h5ad(input_path)
        module = DEGAnalysis(
            project_dir=str(tmp_path),
            params={'groupby': 'leiden', 'scope_key': 'sample_id'},
            progress_callback=lambda p, m: None,
        )
        with pytest.raises(ValueError, match='未选择具体取值'):
            module.run(input_path)


class TestScopeHelper:
    def test_apply_scope_subsets_and_warns(self, tmp_path):
        from modules.base import BaseAnalysis

        class _Dummy(BaseAnalysis):
            MODULE_NAME = 'dummy'

            def run(self, input_path):
                raise NotImplementedError

        adata = _make_sc_anndata()
        module = _Dummy(
            project_dir=str(tmp_path),
            params={'scope_key': 'sample_id', 'scope_values': 'S1,S3'},
            progress_callback=lambda p, m: None,
        )
        out = module.apply_scope(adata)
        assert out.n_obs == 60
        assert set(out.obs['sample_id'].astype(str).unique()) == {'S1', 'S3'}

    def test_apply_scope_ignores_unknown_values_with_warning(self, tmp_path):
        from modules.base import BaseAnalysis

        class _Dummy(BaseAnalysis):
            MODULE_NAME = 'dummy'

            def run(self, input_path):
                raise NotImplementedError

        adata = _make_sc_anndata()
        module = _Dummy(
            project_dir=str(tmp_path),
            params={'scope_key': 'sample_id', 'scope_values': 'S1,NOT_A_SAMPLE'},
            progress_callback=lambda p, m: None,
        )
        out = module.apply_scope(adata)
        assert set(out.obs['sample_id'].astype(str).unique()) == {'S1'}


class TestProportionScopeAndPairs:
    def test_compare_groups_apply_to_sample_level_tests(self, tmp_path):
        import anndata
        import scanpy as _sc

        np.random.seed(7)
        n_obs, n_vars = 120, 100
        raw = np.random.poisson(lam=5, size=(n_obs, n_vars)).astype(np.float32)
        raw[raw == 0] = 1
        obs = pd.DataFrame(index=[f'cell_{i}' for i in range(n_obs)])
        obs['leiden'] = pd.Categorical([str(i % 2) for i in range(n_obs)])
        obs['celltype'] = pd.Categorical(
            ['A'] * (n_obs // 2) + ['B'] * (n_obs - n_obs // 2))
        obs['sample_id'] = ['S1'] * 30 + ['S2'] * 30 + ['S3'] * 30 + ['S4'] * 30
        obs['condition'] = ['C1'] * 60 + ['C2'] * 60
        var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_vars)])
        adata = anndata.AnnData(X=raw.copy(), obs=obs, var=var)
        adata.layers['counts'] = raw.copy()
        _sc.pp.normalize_total(adata, target_sum=1e4)
        _sc.pp.log1p(adata)
        _sc.pp.pca(adata, n_comps=5, svd_solver='arpack')
        _sc.pp.neighbors(adata, n_neighbors=8, n_pcs=5)
        _sc.tl.umap(adata)
        input_path = str(tmp_path / 'input.h5ad')
        adata.write_h5ad(input_path)

        from modules.proportion import ProportionAnalysis

        module = ProportionAnalysis(
            project_dir=str(tmp_path),
            params={
                'groupby': 'celltype',
                'condition_key': 'condition',
                'sample_key': 'sample_id',
                'compare_groups': 'C1-vs-C2',
                'scope_key': 'sample_id',
                'scope_values': 'S1,S2,S3,S4',
                'analysis_unit': 'sample',
                'min_cells_per_sample': 5,
                'min_samples_per_condition': 2,
                'show_proportion_heatmap': False,
            },
            progress_callback=lambda p, m: None,
        )
        result = module.run(input_path)
        sample_tests = pd.read_csv(
            os.path.join(str(tmp_path), 'results', 'sample_level_proportion_tests.csv'))
        assert 'comparison' in sample_tests.columns
        assert set(sample_tests['comparison'].astype(str).unique()) == {'C1 vs C2'}
        assert 'kruskal' not in set(sample_tests['test'].astype(str).unique())
        assert result['summary']['scope_key'] == 'sample_id'
