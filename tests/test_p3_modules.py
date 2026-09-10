"""P3 批量模块单元测试

测试 bulk_timecourse / bulk_deg_integration / bulk_heatmap 中可独立提取的逻辑函数。
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# 辅助工厂
# ---------------------------------------------------------------------------

def _make_timecourse(params=None):
    """创建 BulkTimecourseAnalysis 实例，使用最小参数。"""
    from modules.bulk_timecourse import BulkTimecourseAnalysis
    base = {'time_column': 'time', 'group_column': 'group'}
    if params:
        base.update(params)
    return BulkTimecourseAnalysis(
        project_dir='/tmp',
        params=base,
        progress_callback=lambda p, m: None,
    )


# ===================================================================
# 1. _build_spline_basis
# ===================================================================

class TestBuildSplineBasis:
    """测试 B-spline 基函数矩阵构造。"""

    def test_shape_df3(self):
        """df=3 时输出形状应为 (6, 3)。"""
        tc = _make_timecourse()
        basis = tc._build_spline_basis(np.array([0, 1, 2, 3, 4, 5]), df=3)
        assert basis.shape == (6, 3)

    def test_shape_df4(self):
        """df=4 时输出形状应为 (6, 4)。"""
        tc = _make_timecourse()
        basis = tc._build_spline_basis(np.array([0, 1, 2, 3, 4, 5]), df=4)
        assert basis.shape == (6, 4)

    def test_df_too_small_raises(self):
        """df < 3（degree=3 无截距的最小值）应抛出 ValueError。"""
        tc = _make_timecourse()
        with pytest.raises(ValueError, match='too small'):
            tc._build_spline_basis(np.array([0, 1, 2, 3, 4, 5]), df=2)

    def test_no_nan(self):
        """输出不应包含 NaN。"""
        tc = _make_timecourse()
        basis = tc._build_spline_basis(np.array([0, 1, 2, 3, 4, 5]), df=3)
        assert not np.isnan(basis).any()

    def test_numeric_dtype(self):
        """输出应为数值型数组。"""
        tc = _make_timecourse()
        basis = tc._build_spline_basis(np.array([0, 1, 2, 3, 4, 5]), df=3)
        assert np.issubdtype(basis.dtype, np.floating)

    def test_different_time_values(self):
        """非等间距时间点也应正确生成基矩阵。"""
        tc = _make_timecourse()
        time_vals = np.array([0.0, 0.5, 2.0, 5.0])
        basis = tc._build_spline_basis(time_vals, df=3)
        assert basis.shape == (4, 3)
        assert not np.isnan(basis).any()


# ===================================================================
# 2. _run_temporal_f_test
# ===================================================================

class TestRunTemporalFTest:
    """测试时序 F 检验。"""

    def test_returns_tuple(self):
        """应返回 (F_stats, pvalues) 元组。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        counts = rng.normal(0, 1, size=(20, 10))  # 20 观测 × 10 基因
        time_basis = tc._build_spline_basis(np.arange(20), df=3)
        result = tc._run_temporal_f_test(counts, time_basis, n_spline_cols=3)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_output_length(self):
        """F_stats 和 pvalues 长度应与基因数一致。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        n_genes = 10
        counts = rng.normal(0, 1, size=(20, n_genes))
        time_basis = tc._build_spline_basis(np.arange(20), df=3)
        F_stats, pvalues = tc._run_temporal_f_test(counts, time_basis, n_spline_cols=3)
        assert len(F_stats) == n_genes
        assert len(pvalues) == n_genes

    def test_f_stats_non_negative(self):
        """所有 F 统计量应 >= 0（代码做了 max(F, 0) 截断）。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        counts = rng.normal(0, 1, size=(20, 10))
        time_basis = tc._build_spline_basis(np.arange(20), df=3)
        F_stats, _ = tc._run_temporal_f_test(counts, time_basis, n_spline_cols=3)
        assert np.all(F_stats >= 0)

    def test_pvalues_in_range(self):
        """所有 p 值应在 [0, 1] 之间。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        counts = rng.normal(0, 1, size=(20, 10))
        time_basis = tc._build_spline_basis(np.arange(20), df=3)
        _, pvalues = tc._run_temporal_f_test(counts, time_basis, n_spline_cols=3)
        assert np.all(pvalues >= 0)
        assert np.all(pvalues <= 1)

    def test_time_effect_detected(self):
        """具有明确时间趋势的基因应产生显著 F 统计量。"""
        tc = _make_timecourse()
        n_obs = 30
        t = np.arange(n_obs, dtype=float)
        # 基因 0: 强线性趋势; 基因 1: 纯噪声
        gene_signal = 2.0 * t + np.random.default_rng(0).normal(0, 0.1, n_obs)
        gene_noise = np.random.default_rng(1).normal(0, 1, n_obs)
        counts = np.column_stack([gene_signal, gene_noise])
        time_basis = tc._build_spline_basis(t, df=3)
        F_stats, pvalues = tc._run_temporal_f_test(counts, time_basis, n_spline_cols=3)
        # 信号基因的 F 应远大于噪声基因
        assert F_stats[0] > F_stats[1]
        assert pvalues[0] < pvalues[1]

    def test_insufficient_obs_returns_zeros(self):
        """当观测数不足以拟合模型时应返回全零 F、全一 p。"""
        tc = _make_timecourse()
        # df=3 -> n_spline_cols=3, df_resid = n_obs - 3 - 1, 需要 n_obs > 4
        # n_obs=4 -> df_resid=0, 应走 early return
        counts = np.ones((4, 5))
        time_basis = tc._build_spline_basis(np.array([0, 1, 2, 3]), df=3)
        F_stats, pvalues = tc._run_temporal_f_test(counts, time_basis, n_spline_cols=3)
        assert np.all(F_stats == 0)
        assert np.all(pvalues == 1.0)


# ===================================================================
# 3. _run_interaction_f_test
# ===================================================================

class TestRunInteractionFTest:
    """测试交互效应 F 检验。"""

    def test_returns_tuple(self):
        """应返回 (F_stats, pvalues) 元组。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        n_obs = 20
        n_genes = 8
        n_interaction_cols = 3
        n_full_model_cols = 7  # intercept(1) + spline(3) + group(3)
        # full_basis 不含 intercept，维度 = n_full_model_cols - 1 = 6
        full_basis = rng.normal(0, 1, size=(n_obs, n_full_model_cols - 1))
        counts = rng.normal(0, 1, size=(n_obs, n_genes))
        result = tc._run_interaction_f_test(
            counts, full_basis,
            n_interaction_cols=n_interaction_cols,
            n_full_model_cols=n_full_model_cols,
        )
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_output_length(self):
        """F_stats 和 pvalues 长度应与基因数一致。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        n_obs = 20
        n_genes = 8
        n_full_model_cols = 7
        full_basis = rng.normal(0, 1, size=(n_obs, n_full_model_cols - 1))
        counts = rng.normal(0, 1, size=(n_obs, n_genes))
        F_stats, pvalues = tc._run_interaction_f_test(
            counts, full_basis,
            n_interaction_cols=3,
            n_full_model_cols=n_full_model_cols,
        )
        assert len(F_stats) == n_genes
        assert len(pvalues) == n_genes

    def test_f_stats_non_negative(self):
        """所有 F 统计量应 >= 0。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        n_obs = 20
        n_full_model_cols = 7
        full_basis = rng.normal(0, 1, size=(n_obs, n_full_model_cols - 1))
        counts = rng.normal(0, 1, size=(n_obs, 10))
        F_stats, _ = tc._run_interaction_f_test(
            counts, full_basis,
            n_interaction_cols=3,
            n_full_model_cols=n_full_model_cols,
        )
        assert np.all(F_stats >= 0)

    def test_pvalues_in_range(self):
        """所有 p 值应在 [0, 1] 之间。"""
        tc = _make_timecourse()
        rng = np.random.default_rng(42)
        n_obs = 20
        n_full_model_cols = 7
        full_basis = rng.normal(0, 1, size=(n_obs, n_full_model_cols - 1))
        counts = rng.normal(0, 1, size=(n_obs, 10))
        _, pvalues = tc._run_interaction_f_test(
            counts, full_basis,
            n_interaction_cols=3,
            n_full_model_cols=n_full_model_cols,
        )
        assert np.all(pvalues >= 0)
        assert np.all(pvalues <= 1)

    def test_insufficient_obs_returns_zeros(self):
        """当 df_resid <= 0 时应返回全零 F、全一 p。"""
        tc = _make_timecourse()
        # n_obs=5, n_full_model_cols=7 -> df_resid = 5 - 7 = -2 <= 0
        n_obs = 5
        n_full_model_cols = 7
        full_basis = np.ones((n_obs, n_full_model_cols - 1))
        counts = np.ones((n_obs, 3))
        F_stats, pvalues = tc._run_interaction_f_test(
            counts, full_basis,
            n_interaction_cols=3,
            n_full_model_cols=n_full_model_cols,
        )
        assert np.all(F_stats == 0)
        assert np.all(pvalues == 1.0)


# ===================================================================
# 4. 一致性评分逻辑（内联复制 bulk_deg_integration L125-148）
# ===================================================================

class TestConsistencyScore:
    """测试 DEG 整合模块的一致性评分公式。"""

    @staticmethod
    def _compute_scores(padj_matrix, logfc_matrix, regulation_matrix,
                        comp_names, fc_threshold=2.0, pval_threshold=0.05,
                        min_comparisons=2, exclude_mixed=False):
        """复制 bulk_deg_integration 中的一致性评分逻辑。"""
        log2fc_thresh = np.log2(fc_threshold)
        consistency_scores = []
        for g in padj_matrix.index:
            pvals = padj_matrix.loc[g, comp_names].values
            abs_fcs = np.abs(logfc_matrix.loc[g, comp_names].values)
            sig_mask = (pvals < pval_threshold) & (abs_fcs >= log2fc_thresh)
            sig_count = int(sig_mask.sum())
            if sig_count >= min_comparisons:
                sig_signs = regulation_matrix.loc[g, comp_names].values[sig_mask]
                sig_neg_log_p = -np.log10(pvals[sig_mask] + 1e-300)
                if np.all(sig_signs == sig_signs[0]):
                    score = float(sig_signs[0] * np.mean(sig_neg_log_p))
                else:
                    score = 0.0
                if exclude_mixed and score == 0:
                    continue
                consistency_scores.append({
                    'gene': g,
                    'consistency_score': round(score, 4),
                    'n_significant': sig_count,
                    'direction': 'Up' if score > 0 else 'Down' if score < 0 else 'Mixed',
                })
        return consistency_scores

    def test_concordant_upregulated_high_score(self):
        """在两个比较中均上调且显著的基因应获得正的一致性评分。"""
        genes = ['GeneA', 'GeneB', 'GeneC']
        comps = ['comp1', 'comp2']
        # GeneA: 两个比较均显著上调
        # GeneB: 仅一个比较显著
        # GeneC: 两个比较方向不一致
        padj = pd.DataFrame(
            [[0.001, 0.002], [0.001, 0.8], [0.001, 0.001]],
            index=genes, columns=comps,
        )
        logfc = pd.DataFrame(
            [[2.0, 3.0], [2.0, 0.1], [2.0, -2.0]],
            index=genes, columns=comps,
        )
        regulation = pd.DataFrame(
            [[1, 1], [1, 0], [1, -1]],
            index=genes, columns=comps,
        )
        scores = self._compute_scores(
            padj, logfc, regulation, comps,
            fc_threshold=2.0, pval_threshold=0.05, min_comparisons=2,
        )
        score_dict = {s['gene']: s for s in scores}
        # GeneA 应为正分 (Up)
        assert 'GeneA' in score_dict
        assert score_dict['GeneA']['consistency_score'] > 0
        assert score_dict['GeneA']['direction'] == 'Up'
        assert score_dict['GeneA']['n_significant'] == 2

    def test_concordant_downregulated_negative_score(self):
        """在两个比较中均下调且显著的基因应获得负的一致性评分。"""
        genes = ['GeneD']
        comps = ['comp1', 'comp2']
        padj = pd.DataFrame([[0.001, 0.002]], index=genes, columns=comps)
        logfc = pd.DataFrame([[-2.0, -3.0]], index=genes, columns=comps)
        regulation = pd.DataFrame([[-1, -1]], index=genes, columns=comps)
        scores = self._compute_scores(
            padj, logfc, regulation, comps,
            fc_threshold=2.0, pval_threshold=0.05, min_comparisons=2,
        )
        assert len(scores) == 1
        assert scores[0]['consistency_score'] < 0
        assert scores[0]['direction'] == 'Down'

    def test_mixed_direction_zero_score(self):
        """方向不一致的基因应得分为 0，direction 为 Mixed。"""
        genes = ['GeneE']
        comps = ['comp1', 'comp2']
        padj = pd.DataFrame([[0.001, 0.001]], index=genes, columns=comps)
        logfc = pd.DataFrame([[2.0, -2.0]], index=genes, columns=comps)
        regulation = pd.DataFrame([[1, -1]], index=genes, columns=comps)
        scores = self._compute_scores(
            padj, logfc, regulation, comps,
            fc_threshold=2.0, pval_threshold=0.05, min_comparisons=2,
        )
        assert len(scores) == 1
        assert scores[0]['consistency_score'] == 0.0
        assert scores[0]['direction'] == 'Mixed'

    def test_exclude_mixed_removes_gene(self):
        """exclude_mixed=True 时方向不一致的基因应被排除。"""
        genes = ['GeneF']
        comps = ['comp1', 'comp2']
        padj = pd.DataFrame([[0.001, 0.001]], index=genes, columns=comps)
        logfc = pd.DataFrame([[2.0, -2.0]], index=genes, columns=comps)
        regulation = pd.DataFrame([[1, -1]], index=genes, columns=comps)
        scores = self._compute_scores(
            padj, logfc, regulation, comps,
            fc_threshold=2.0, pval_threshold=0.05, min_comparisons=2,
            exclude_mixed=True,
        )
        assert len(scores) == 0

    def test_non_significant_gene_excluded(self):
        """不满足显著性阈值的基因应被排除。"""
        genes = ['GeneG']
        comps = ['comp1', 'comp2']
        # padj 值很高，不显著
        padj = pd.DataFrame([[0.5, 0.6]], index=genes, columns=comps)
        logfc = pd.DataFrame([[2.0, 3.0]], index=genes, columns=comps)
        regulation = pd.DataFrame([[1, 1]], index=genes, columns=comps)
        scores = self._compute_scores(
            padj, logfc, regulation, comps,
            fc_threshold=2.0, pval_threshold=0.05, min_comparisons=2,
        )
        assert len(scores) == 0

    def test_min_comparisons_filter(self):
        """min_comparisons 门槛应正确过滤。"""
        genes = ['GeneH']
        comps = ['comp1', 'comp2', 'comp3']
        padj = pd.DataFrame([[0.001, 0.001, 0.5]], index=genes, columns=comps)
        logfc = pd.DataFrame([[2.0, 3.0, 2.0]], index=genes, columns=comps)
        regulation = pd.DataFrame([[1, 1, 1]], index=genes, columns=comps)
        # min_comparisons=3，但只有 2 个显著 -> 应排除
        scores_3 = self._compute_scores(
            padj, logfc, regulation, comps,
            fc_threshold=2.0, pval_threshold=0.05, min_comparisons=3,
        )
        assert len(scores_3) == 0
        # min_comparisons=2 -> 应通过
        scores_2 = self._compute_scores(
            padj, logfc, regulation, comps,
            fc_threshold=2.0, pval_threshold=0.05, min_comparisons=2,
        )
        assert len(scores_2) == 1


# ===================================================================
# 5. clip_range 解析逻辑（内联复制 bulk_heatmap L167-176）
# ===================================================================

class TestClipRangeParsing:
    """测试热图 clip_range 参数解析。"""

    @staticmethod
    def _parse_clip_range(clip_str):
        """复制 bulk_heatmap 中的 clip_range 解析逻辑。
        返回 (min, max) 浮点元组或 None，或抛出 ValueError。
        """
        clip_str = clip_str.strip()
        clip_range = None
        if clip_str:
            parts = [x.strip() for x in clip_str.split(',') if x.strip()]
            if len(parts) != 2:
                raise ValueError(
                    f"clip_range 格式错误: '{clip_str}'，应为 'min,max'，如 '-3,3'"
                )
            try:
                clip_range = (float(parts[0]), float(parts[1]))
            except ValueError:
                raise ValueError(f"clip_range 数值解析失败: '{clip_str}'")
        return clip_range

    def test_positive_range(self):
        """解析 '0,5' 应返回 (0.0, 5.0)。"""
        assert self._parse_clip_range('0,5') == (0.0, 5.0)

    def test_negative_range(self):
        """解析 '-1,1' 应返回 (-1.0, 1.0)。"""
        assert self._parse_clip_range('-1,1') == (-1.0, 1.0)

    def test_default_range(self):
        """解析 '-3,3' 应返回 (-3.0, 3.0)。"""
        assert self._parse_clip_range('-3,3') == (-3.0, 3.0)

    def test_whitespace_handling(self):
        """带空格的输入应正确解析。"""
        assert self._parse_clip_range('  -2 , 4  ') == (-2.0, 4.0)

    def test_empty_string_returns_none(self):
        """空字符串应返回 None。"""
        assert self._parse_clip_range('') is None

    def test_invalid_format_raises(self):
        """非 'min,max' 格式应抛出 ValueError。"""
        with pytest.raises(ValueError, match='clip_range 格式错误'):
            self._parse_clip_range('1,2,3')

    def test_single_value_raises(self):
        """只有一个数值应抛出 ValueError。"""
        with pytest.raises(ValueError, match='clip_range 格式错误'):
            self._parse_clip_range('5')

    def test_non_numeric_raises(self):
        """非数值内容应抛出 ValueError。"""
        with pytest.raises(ValueError, match='clip_range'):
            self._parse_clip_range('abc,def')

    def test_float_values(self):
        """浮点数应正确解析。"""
        assert self._parse_clip_range('-1.5,2.5') == (-1.5, 2.5)


class TestHeatmapSampleSelection:
    """Top-var 选基因与热图展示必须共享同一个样品范围。"""

    @pytest.fixture
    def obs(self):
        return pd.DataFrame(
            {'group': ['ctrl', 'ctrl', 'treat', 'treat', 'dose', 'dose']},
            index=['s1', 's2', 's3', 's4', 's5', 's6'],
        )

    def test_deg_groups_selects_only_the_contrast_samples(self, obs):
        from modules.bulk_heatmap import _select_heatmap_display_samples

        indices, contract = _select_heatmap_display_samples(
            obs, sample_display_mode='deg_groups', groupby='group',
            comparison_label='treat vs ctrl',
        )

        assert indices == [0, 1, 2, 3]
        assert contract['selected_groups'] == ['treat', 'ctrl']
        assert contract['mode'] == 'deg_groups'

    def test_selected_groups_keeps_the_requested_groups_only(self, obs):
        from modules.bulk_heatmap import _select_heatmap_display_samples

        indices, contract = _select_heatmap_display_samples(
            obs, sample_display_mode='selected_groups', groupby='group',
            selected_groups='ctrl,dose',
        )

        assert indices == [0, 1, 4, 5]
        assert contract['n_samples'] == 4

    def test_selected_samples_requires_exact_sample_names(self, obs):
        from modules.bulk_heatmap import _select_heatmap_display_samples

        indices, _ = _select_heatmap_display_samples(
            obs, sample_display_mode='selected_samples', selected_samples='s2,s5',
        )
        assert indices == [1, 4]
        with pytest.raises(ValueError, match='未找到所选样本'):
            _select_heatmap_display_samples(
                obs, sample_display_mode='selected_samples', selected_samples='s2,missing',
            )

    @pytest.mark.parametrize(
        ('requested', 'source', 'comparison', 'expected'),
        [
            ('auto', 'top_var', '', 'all'),
            ('auto', 'top_var', 'treat vs ctrl', 'deg_groups'),
            ('auto', 'deg', 'treat vs ctrl', 'deg_groups'),
            ('all', 'deg', 'treat vs ctrl', 'all'),
        ],
    )
    def test_auto_scope_resolution_preserves_explicit_overrides(
            self, requested, source, comparison, expected):
        from modules.bulk_heatmap import _resolve_heatmap_sample_mode

        requested_mode, effective_mode = _resolve_heatmap_sample_mode(
            requested, source, comparison,
        )

        assert requested_mode == requested
        assert effective_mode == expected

    def test_top_variable_genes_are_ranked_only_in_selected_samples(self):
        from modules.bulk_heatmap import _select_top_variable_gene_indices

        # gene_global only varies in the excluded third group; gene_contrast
        # varies inside the selected ctrl/treat comparison.
        matrix = np.array([
            [0.0, 0.0],
            [0.0, 1.0],
            [0.0, 10.0],
            [0.0, 11.0],
            [100.0, 5.0],
            [-100.0, 5.0],
        ])

        top_idx, scores = _select_top_variable_gene_indices(
            matrix, [0, 1, 2, 3], metric='var', top_n=1,
        )

        assert top_idx == [1]
        assert scores[0] == 0.0
        assert scores[1] > 0.0


def test_deg_integration_deduplicates_gene_symbols():
    """多个 Ensembl ID 映射到同一 symbol 时保留最显著的一行。"""
    import pandas as pd
    from modules.bulk_deg_integration import _deduplicate_gene_results

    df = pd.DataFrame({
        'gene': ['DUP', 'DUP', 'UNIQUE'],
        'log2FC': [1.2, -2.0, 0.5],
        'padj': [0.01, 0.001, 0.2],
        'regulation': ['Up', 'Down', 'NS'],
    })
    result = _deduplicate_gene_results(df)

    assert result['gene'].is_unique
    assert result.loc[result['gene'] == 'DUP', 'log2FC'].item() == -2.0


def test_deg_integration_loads_standalone_comparison_labels(tmp_path):
    """独立 runner 无数据库记录时也应使用真实比较名。"""
    import json
    from modules.bulk_deg_integration import _load_comparison_labels

    results = tmp_path / 'results'
    results.mkdir()
    (results / 'bulk_deg_results_0.csv').write_text('gene,log2FC,padj\nA,1,0.01\n')
    (results / 'bulk_deg_comparison_labels.json').write_text(json.dumps({
        'bulk_deg_results_0.csv': 'Treatment-vs-Control',
    }))

    assert _load_comparison_labels(str(tmp_path)) == {
        'bulk_deg_results_0.csv': 'Treatment-vs-Control',
    }


def test_gsea_export_preserves_pathway_names_from_index():
    import pandas as pd
    from modules.bulk_enrichment import _ensure_gsea_term_column

    df = pd.DataFrame({'nes': [2.1], 'fdr': [0.01]}, index=['interferon signaling'])
    result = _ensure_gsea_term_column(df)

    assert result.loc[0, 'Term'] == 'interferon signaling'


def test_enrichment_resolves_missing_kegg_library_once(monkeypatch, tmp_path):
    """KEGG is fetched on demand rather than relying on OmicVerse's GO-only downloader."""
    import modules.bulk_enrichment as enrichment

    monkeypatch.chdir(tmp_path)
    calls = []

    def fake_download(library_name, destination):
        calls.append(library_name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text('KEGG pathway\t\tGENE1\tGENE2\n')

    monkeypatch.setattr(enrichment, '_download_enrichr_geneset', fake_download)
    first = enrichment._resolve_geneset_path('KEGG_2021_Human', 'Human')
    second = enrichment._resolve_geneset_path('KEGG_2021_Human', 'Human')

    assert first == 'genesets/KEGG_2021_Human.txt'
    assert second == first
    assert calls == ['KEGG_2021_Human']


def test_enrichment_uses_requested_database_when_gene_set_is_internal_label():
    from modules.bulk_enrichment import _enrichment_ontology

    assert _enrichment_ontology('gs_ind', 'KEGG') == 'KEGG'


def test_enrichment_integration_retains_source_metadata(tmp_path):
    import pandas as pd
    from modules.bulk_enrichment import _write_enrichment_integration

    results_dir = tmp_path / 'results'
    results_dir.mkdir()
    pd.DataFrame({
        'Database': ['GO_BP'], 'Method': ['ORA'], 'Direction': ['All'],
        'Term': ['response to virus'], 'Adjusted P-value': [0.02],
    }).to_csv(results_dir / 'enrichment_ora_go_bp_results.csv', index=False)
    pd.DataFrame({
        'Database': ['KEGG'], 'Method': ['ORA'], 'Direction': ['Up'],
        'Term': ['MAPK signaling pathway'], 'Adjusted P-value': [0.001],
    }).to_csv(results_dir / 'enrichment_ora_kegg_directional_results.csv', index=False)

    files = _write_enrichment_integration(results_dir)
    integrated = pd.read_csv(results_dir / 'enrichment_integrated_results.csv')

    assert len(files) >= 1
    assert integrated['Database'].tolist() == ['KEGG', 'GO_BP']
    assert {'Comparison', 'Method', 'Direction', 'Term', 'Source result'}.issubset(integrated.columns)
    assert set(integrated['Comparison']) == {'Unspecified'}


def test_enrichment_integration_builds_one_overview_for_same_comparison(tmp_path):
    import os
    import pandas as pd
    from modules.bulk_enrichment import _write_enrichment_integration

    results_dir = tmp_path / 'results'
    results_dir.mkdir()
    common = {'Comparison': ['Ctrl vs Treat', 'Ctrl vs Treat'], 'Method': ['ORA', 'ORA'],
              'Direction': ['All', 'All']}
    pd.DataFrame({
        **common, 'Database': ['GO_BP', 'GO_BP'],
        'Term': ['response to lipid', 'inflammatory response'],
        'Adjusted P-value': [0.001, 0.008], 'Overlap': ['5/100', '3/100'],
    }).to_csv(results_dir / 'enrichment_ora_go_bp_ctrl_vs_treat_results.csv', index=False)
    pd.DataFrame({
        **common, 'Database': ['KEGG', 'KEGG'],
        'Term': ['PPAR signaling pathway', 'Fatty acid metabolism'],
        'Adjusted P-value': [0.002, 0.02], 'Overlap': ['4/90', '2/75'],
    }).to_csv(results_dir / 'enrichment_ora_kegg_ctrl_vs_treat_results.csv', index=False)

    files = _write_enrichment_integration(results_dir)
    overview_files = [item for item in files if item['category'] == 'enrichment_overview']
    integrated = pd.read_csv(results_dir / 'enrichment_integrated_results.csv')

    assert len(integrated) == 4
    assert set(integrated['Comparison']) == {'Ctrl vs Treat'}
    assert {item['file_type'] for item in overview_files} == {'png', 'svg', 'pdf'}
    assert all(os.path.isfile(item['file_path']) for item in overview_files)


def test_bulk_go_priority_triptych_exports_selection_audit_and_source_table(tmp_path):
    import json
    import os
    import pandas as pd
    from modules.bulk_enrichment import (
        _go_priority_config_from_params,
        _write_enrichment_integration,
    )

    results_dir = tmp_path / 'results'
    plots_dir = tmp_path / 'plots'
    results_dir.mkdir()
    common = {
        'Comparison': ['Ctrl vs Treat'], 'Method': ['ORA'], 'Direction': ['Up'],
        'Significant': [True], 'Overlap': ['4/100'], 'Genes': ['IL6;TNF'],
    }
    for database, term, fdr in (
        ('GO_BP', 'inflammatory response', .001),
        ('GO_CC', 'lipid particle', .004),
        ('GO_MF', 'lipid binding', .009),
    ):
        pd.DataFrame({
            **common, 'Database': [database], 'Term': [term],
            'Adjusted P-value': [fdr],
        }).to_csv(results_dir / f'enrichment_ora_{database.lower()}_results.csv', index=False)

    source_paths = sorted(results_dir.glob('enrichment_ora_*_results.csv'))
    config = _go_priority_config_from_params(
        {
            'go_priority_enabled': True, 'method': 'ORA',
            'focus_terms': 'lipid binding', 'focus_label': '炎症与脂代谢',
        },
        databases=['GO_BP', 'GO_CC', 'GO_MF'], top_n=8,
    )
    files = _write_enrichment_integration(
        results_dir, source_paths=source_paths, plots_dir=plots_dir,
        output_prefix='theme_go', go_priority_config=config,
    )

    categories = [item['category'] for item in files]
    assert categories.count('go_priority_dotplot') == 3
    assert categories.count('go_priority_barplot') == 3
    source_path = results_dir / 'theme_go_go_priority_source.csv'
    audit_path = results_dir / 'theme_go_go_priority_selection_audit.json'
    selected = pd.read_csv(source_path)
    audit = json.loads(audit_path.read_text(encoding='utf-8'))
    assert set(selected['Database']) == {'GO_BP', 'GO_CC', 'GO_MF'}
    assert set(selected['selection_reason']) >= {
        'inflammation_priority', 'lipid_metabolism_priority',
    }
    assert audit['units'][0]['status'] == 'completed'
    assert all(os.path.isfile(item['file_path']) for item in files)


def test_bulk_human_pathway_triptych_exports_fdr_ranked_source_and_audit(tmp_path):
    import json
    import os
    import pandas as pd
    from modules.bulk_enrichment import (
        _pathway_triptych_config_for_batch,
        _write_enrichment_integration,
    )

    results_dir = tmp_path / 'results'
    plots_dir = tmp_path / 'plots'
    results_dir.mkdir()
    common = {
        'Comparison': ['Ctrl vs Treat', 'Ctrl vs Treat'],
        'Method': ['ORA', 'ORA'], 'Direction': ['Up', 'Up'],
        'Significant': [True, True], 'Overlap': ['4/100', '3/100'],
        'Genes': ['IL6;TNF', 'APOA1;APOC3'],
    }
    for database, terms, fdrs in (
        ('KEGG', ['TNF signaling pathway', 'Fatty acid metabolism'], [.001, .008]),
        ('Reactome', ['Immune System', 'Metabolism of lipids'], [.002, .009]),
        ('WikiPathways', ['Inflammation', 'Lipid metabolism'], [.003, .01]),
    ):
        pd.DataFrame({
            **common, 'Database': [database, database], 'Term': terms,
            'Adjusted P-value': fdrs,
        }).to_csv(results_dir / f'enrichment_ora_{database.lower()}_results.csv', index=False)

    config = _pathway_triptych_config_for_batch(
        {'method': 'ORA'}, databases=['KEGG', 'Reactome', 'WikiPathways'], top_n=8,
    )
    files = _write_enrichment_integration(
        results_dir, source_paths=sorted(results_dir.glob('enrichment_ora_*_results.csv')),
        plots_dir=plots_dir, output_prefix='human_pathways',
        pathway_triptych_config=config,
    )

    categories = [item['category'] for item in files]
    assert categories.count('pathway_triptych_dotplot') == 3
    assert categories.count('pathway_triptych_barplot') == 3
    selected = pd.read_csv(results_dir / 'human_pathways_pathway_triptych_source.csv')
    audit = json.loads((results_dir / 'human_pathways_pathway_triptych_selection_audit.json').read_text())
    assert set(selected['Database']) == {'KEGG', 'Reactome', 'WikiPathways'}
    assert set(selected['selection_reason']) == {'fdr_top_pathway'}
    assert audit['units'][0]['status'] == 'completed'
    assert all(os.path.isfile(item['file_path']) for item in files)


def test_single_cell_human_pathway_triptych_selection_keeps_database_fdr_separate():
    import pandas as pd
    from modules.sc_cell_go import _select_pathway_triptych_rows

    rows = []
    for library, terms in {
        'KEGG_2021_Human': [('KEGG top', .001), ('KEGG second', .02)],
        'Reactome_2022': [('Reactome top', .002), ('Reactome second', .03)],
        'WikiPathway_2021_Human': [('Wiki top', .003), ('Wiki second', .04)],
    }.items():
        for term, fdr in terms:
            rows.append({
                'status': 'completed', 'Significant': True, 'gene_set': library,
                'Term': term, 'method': 'ORA', 'direction': 'Up',
                'Adjusted P-value': fdr, 'Overlap': '4/100',
            })
    selected, audit = _select_pathway_triptych_rows(pd.DataFrame(rows), top_n=1)

    assert selected['Term'].tolist() == ['KEGG top', 'Reactome top', 'Wiki top']
    assert selected['selection_reason'].eq('fdr_top_pathway').all()
    assert set(audit['databases']) == {'KEGG', 'Reactome', 'WikiPathways'}


def test_enrichment_text_audit_detects_overlapping_labels():
    import matplotlib.pyplot as plt
    from modules.bulk_enrichment import _audit_figure_text_overlap

    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, 'first label')
    ax.text(0.5, 0.5, 'second label')
    try:
        audit = _audit_figure_text_overlap(fig)
        assert audit['n_overlap_pairs'] >= 1
        assert audit['status'] == 'warning'
    finally:
        plt.close(fig)


def test_enrichment_figure_uses_ontology_colors_and_overlap_counts():
    import pandas as pd
    from modules.bulk_enrichment import _enrichment_figure

    df = pd.DataFrame({
        'Term': ['response to virus', 'mRNA binding', 'Measles'],
        'P-value': [1e-6, 1e-3, 1e-4],
        'Overlap': ['6/100', '2/100', '4/100'],
        'Ontology': ['BP', 'MF', 'KEGG'],
    })
    fig = _enrichment_figure(df, 'enrich result', database='GO_BP')
    try:
        assert fig is not None
        assert len(fig.axes) == 1
        assert len(fig.axes[0].patches) == 3
        assert any(text.get_text() == 'n=6' for text in fig.axes[0].texts)
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)
