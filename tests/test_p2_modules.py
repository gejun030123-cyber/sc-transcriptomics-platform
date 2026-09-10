# tests/test_p2_modules.py
"""Tests for P2 single-cell modules — 常量、标记基因、降维、差异表达、质控、HVG、输入校验。"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import pytest
import anndata


# ============================================================
# 1. constants.py — S_GENES, G2M_GENES
# ============================================================

class TestCellCycleConstants:
    """测试细胞周期基因列表常量。"""

    def test_s_genes_not_empty(self):
        """S_GENES 是非空列表。"""
        from modules.constants import S_GENES
        assert isinstance(S_GENES, list)
        assert len(S_GENES) > 0

    def test_g2m_genes_not_empty(self):
        """G2M_GENES 是非空列表。"""
        from modules.constants import G2M_GENES
        assert isinstance(G2M_GENES, list)
        assert len(G2M_GENES) > 0

    def test_no_overlap(self):
        """S_GENES 和 G2M_GENES 无交集。"""
        from modules.constants import S_GENES, G2M_GENES
        overlap = set(S_GENES) & set(G2M_GENES)
        assert overlap == set(), f"重叠基因: {overlap}"

    def test_known_genes_present(self):
        """已知标记基因存在：PCNA 在 S_GENES，MKI67 在 G2M_GENES。"""
        from modules.constants import S_GENES, G2M_GENES
        assert 'PCNA' in S_GENES
        assert 'MKI67' in G2M_GENES


# ============================================================
# 2. annotation.py — MARKER_SETS, DEFAULT_TME_MARKERS
# ============================================================

class TestAnnotationMarkers:
    """测试注释模块的标记基因集。"""

    def test_marker_sets_keys(self):
        """MARKER_SETS 包含内置注释基因集。"""
        from modules.annotation import MARKER_SETS
        assert set(MARKER_SETS.keys()) == {
            'Universal', 'Colorectal', 'Colorectal_refined', 'TME', 'Immune', 'Blood', 'PBMC',
        }

    def test_organoid_marker_sets_cover_common_tissues(self):
        """类器官面板按组织拆分，且每个类型至少有多个 marker。"""
        from modules.annotation import ORGANOID_MARKER_SETS, get_marker_set
        expected = {'intestinal', 'cerebral', 'kidney', 'liver', 'lung', 'pancreatic', 'cardiac'}
        assert set(ORGANOID_MARKER_SETS) == expected
        for organoid_type in expected:
            panel = get_marker_set('Organoid', organoid_type)
            assert panel
            assert all(len(genes) >= 4 for genes in panel.values())
        assert 'LGR5' in get_marker_set('Organoid', 'intestinal')['Intestinal stem cells']
        assert 'NPHS2' in get_marker_set('Organoid', 'kidney')['Podocytes']

    def test_tme_markers_not_empty(self):
        """DEFAULT_TME_MARKERS 中每个细胞类型的基因列表非空。"""
        from modules.annotation import DEFAULT_TME_MARKERS
        for cell_type, genes in DEFAULT_TME_MARKERS.items():
            assert len(genes) > 0, f"{cell_type} 的标记基因列表为空"

    def test_marker_genes_are_strings(self):
        """所有标记基因名称都是字符串类型。"""
        from modules.annotation import MARKER_SETS
        for set_name, marker_dict in MARKER_SETS.items():
            for cell_type, genes in marker_dict.items():
                for gene in genes:
                    assert isinstance(gene, str), (
                        f"{set_name}/{cell_type} 中存在非字符串基因: {gene!r}"
                    )

    def test_marker_resolution_is_case_insensitive_and_applies_coverage_gate(self):
        """鼠/人符号大小写差异应匹配，低覆盖类型不可参与打分。"""
        from modules.annotation import resolve_marker_genes
        usable, coverage = resolve_marker_genes(
            {'T': ['CD3D', 'TRAC'], 'B': ['MS4A1', 'CD79A']},
            ['Cd3d', 'Trac', 'Ms4a1'], min_markers_per_type=2,
        )
        assert usable == {'T': ['Cd3d', 'Trac']}
        assert coverage['B'] == {'matched': 1, 'total': 2}

    def test_cluster_marker_selector_is_data_driven_and_excludes_housekeeping(self):
        """Dotplot genes should be cluster-specific, filtered, and globally unique."""
        import anndata
        pytest.importorskip('scanpy')
        from modules.annotation import select_cluster_marker_genes

        rng = np.random.default_rng(4)
        genes = ['A_MARK', 'B_MARK', 'C_MARK', 'ACTB', 'RPL13', 'MT-GENE'] + [f'G{i}' for i in range(12)]
        matrix = rng.poisson(0.05, size=(36, len(genes))).astype(float)
        for start, gene_index in ((0, 0), (12, 1), (24, 2)):
            matrix[start:start + 12, gene_index] += 8
        adata = anndata.AnnData(
            matrix,
            obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 12 + ['1'] * 12 + ['2'] * 12)}),
            var=pd.DataFrame(index=genes),
        )

        selected = select_cluster_marker_genes(adata, 'leiden', classic_markers={})
        by_cluster = selected['cluster_markers']
        assert {'0', '1', '2'} <= set(by_cluster)
        assert {'A_MARK', 'B_MARK', 'C_MARK'} <= set(selected['genes'])
        assert set(selected['genes']).isdisjoint({'ACTB', 'RPL13', 'MT-GENE'})
        assert len(selected['genes']) == len(set(selected['genes']))
        assert selected['n_data_driven'] >= 3
        assert all(
            item['source'] == 'data_driven'
            for details in selected['marker_details'].values()
            for item in details
        )


# ============================================================
# 3. dimred.py — 肘部法 PC 选择逻辑
# ============================================================

class TestElbowLogic:
    """测试肘部法自动选择主成分数量的逻辑。"""

    @staticmethod
    def _elbow_select(variance_ratio):
        """调用降维模块实际使用的肘部法逻辑。"""
        from modules.dimred import _select_elbow_n_comps

        return _select_elbow_n_comps(variance_ratio, len(variance_ratio))

    def test_elbow_with_plateau(self):
        """方差比在索引 5 附近平台化 → 选择约 5-7 个 PC。"""
        # 前 6 个 PC 有明显递减的方差比，之后趋于平坦
        variance_ratio = np.array([
            0.15, 0.10, 0.08, 0.06, 0.05, 0.04,
            0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
            0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
            0.02, 0.02,
        ])
        result = self._elbow_select(variance_ratio)
        assert 4 <= result <= 8, f"肘部法结果 {result} 不在预期范围 [4, 8]"

    def test_clamping_min(self):
        """结果应 >= 5（下限钳制）。"""
        # 所有方差比相同 → diffs 全零 → argmax 返回 0 → 0 + 2 = 2 → clamp to 5
        variance_ratio = np.ones(10) * 0.1
        result = self._elbow_select(variance_ratio)
        assert result >= 5, f"结果 {result} 低于下限 5"

    def test_clamping_max(self):
        """结果应 <= n_comps（上限钳制）。"""
        # 陡峭递减但 n_comps 很小
        n_comps = 5
        variance_ratio = np.array([0.3, 0.25, 0.20, 0.15, 0.10])
        result = self._elbow_select(variance_ratio)
        assert result <= n_comps, f"结果 {result} 超过 n_comps={n_comps}"

    def test_clear_elbow_at_index_4(self):
        """在索引 4 处有明显肘部弯曲。"""
        # 前 5 个 PC 下降快，之后趋平
        variance_ratio = np.array([
            0.20, 0.15, 0.10, 0.08, 0.06,
            0.02, 0.02, 0.02, 0.02, 0.02,
            0.02, 0.02, 0.02, 0.02, 0.02,
        ])
        result = self._elbow_select(variance_ratio)
        # argmax of diffs2 should be at a small index, result ~ 4-6
        assert 3 <= result <= 7, f"肘部法结果 {result} 不在预期范围 [3, 7]"


# ============================================================
# 4. deg.py — correction_map 和显著性掩码
# ============================================================

class TestDEGLogic:
    """测试差异表达分析的辅助逻辑。"""

    def test_three_colour_volcano_classification(self):
        """火山图应按 padj 和 logFC 分为上调、下调和不显著三类。"""
        from modules.deg import classify_volcano_regulation, VOLCANO_COLOR_MAP

        labels = classify_volcano_regulation(
            logfc=[2.0, -2.0, 0.5, 2.0, -2.0, np.nan],
            padj=[0.01, 0.01, 0.01, 0.10, 0.10, 0.01],
            logfc_cutoff=1.0,
            pval_cutoff=0.05,
        )

        assert labels.tolist() == ['Up', 'Down', 'NS', 'NS', 'NS', 'NS']
        assert set(labels) <= set(VOLCANO_COLOR_MAP)
        assert VOLCANO_COLOR_MAP['Up'] == '#B64342'
        assert VOLCANO_COLOR_MAP['Down'] == '#0F4D92'
        assert VOLCANO_COLOR_MAP['NS'] == '#98A2B3'

    def test_detection_fraction_is_joined_by_gene_not_rank_position(self):
        """Scanpy pts uses var order and must not be zipped to ranked names."""
        from modules.deg import max_detection_fraction_by_gene

        rank_result = {
            'names': np.array(
                [('A_MARKER', 'B_MARKER'), ('UNUSED', 'UNUSED'),
                 ('B_MARKER', 'A_MARKER')],
                dtype=[('A', object), ('B', object)],
            ),
            'pts': pd.DataFrame(
                {'A': [0.0, 0.0, 1.0], 'B': [1.0, 0.0, 0.0]},
                index=['B_MARKER', 'UNUSED', 'A_MARKER'],
            ),
        }

        detected = max_detection_fraction_by_gene(rank_result)

        assert detected['A_MARKER'] == 1.0
        assert detected['B_MARKER'] == 1.0
        assert detected['UNUSED'] == 0.0

    def test_correction_map(self):
        """校正方法映射正确。"""
        correction_map = {
            'benjamini_hochberg': 'fdr_bh',
            'bonferroni': 'bonferroni',
            'BY': 'fdr_by',
        }
        assert correction_map['benjamini_hochberg'] == 'fdr_bh'
        assert correction_map['bonferroni'] == 'bonferroni'
        assert correction_map['BY'] == 'fdr_by'

    def test_significance_mask(self):
        """显著性掩码：pval_adj < cutoff 且 |logfc| > cutoff。"""
        df = pd.DataFrame({
            'gene': ['A', 'B', 'C', 'D'],
            'pval_adj': [0.01, 0.10, 0.03, 0.08],
            'logfc': [2.0, 1.5, 0.5, 3.0],
        })
        pval_cutoff = 0.05
        logfc_cutoff = 1.0
        sig = (df['pval_adj'] < pval_cutoff) & (df['logfc'].abs() > logfc_cutoff)

        # A: pval=0.01<0.05, |logfc|=2.0>1.0 → True
        # B: pval=0.10>=0.05 → False
        # C: pval=0.03<0.05, |logfc|=0.5<=1.0 → False
        # D: pval=0.08>=0.05 → False
        expected = pd.Series([True, False, False, False])
        pd.testing.assert_series_equal(sig, expected)

    def test_significance_mask_all_significant(self):
        """所有基因均显著。"""
        df = pd.DataFrame({
            'pval_adj': [0.001, 0.002, 0.003],
            'logfc': [3.0, 2.0, 1.5],
        })
        sig = (df['pval_adj'] < 0.05) & (df['logfc'].abs() > 1.0)
        assert sig.all()

    def test_significance_mask_none_significant(self):
        """无基因显著。"""
        df = pd.DataFrame({
            'pval_adj': [0.10, 0.20, 0.30],
            'logfc': [0.1, 0.2, 0.3],
        })
        sig = (df['pval_adj'] < 0.05) & (df['logfc'].abs() > 1.0)
        assert not sig.any()


# ============================================================
# 5. qc.py — 基因标记正则匹配
# ============================================================

class TestGeneFlagging:
    """测试质控模块的基因标记逻辑（复现 qc.py 中的正则）。"""

    def test_mito_flagging(self):
        """以 'MT-' 开头的基因应被标记为线粒体基因。"""
        genes = pd.Series(['MT-ND1', 'MT-CO1', 'GAPDH', 'MT-ATP6', 'ACTB'])
        result = genes.str.startswith('MT-')
        expected = pd.Series([True, True, False, True, False])
        pd.testing.assert_series_equal(result, expected)

    def test_ribo_flagging(self):
        """只标记核糖体蛋白，不误标 RPS6 激酶。"""
        from modules.qc import _ribosomal_gene_mask

        genes = ['RPS2', 'RPL13A', 'GAPDH', 'RPS27', 'RPLP0', 'RPS6KA1']
        result = _ribosomal_gene_mask(genes)

        assert result == [True, True, False, True, True, False]

    def test_hb_flagging(self):
        """'HBA1' 应匹配 hb 模式，'HPRT1' 不应匹配。"""
        # 复现 qc.py 中的正则: ^HB[^(P)]
        genes = pd.Series(['HBA1', 'HBB', 'HPRT1', 'HBZ', 'GAPDH'])
        result = genes.str.contains('^HB[^(P)]')
        # HBA1: H + B + A (不是 P) → True
        # HBB: H + B + B (不是 P) → True
        # HPRT1: H + P → 不匹配（^HB 要求第二个字符是 B）
        # HBZ: H + B + Z (不是 P) → True
        # GAPDH: 不以 HB 开头 → False
        expected = pd.Series([True, True, False, True, False])
        pd.testing.assert_series_equal(result, expected)


# ============================================================
# 6. hvg.py — 强制基因解析
# ============================================================

class TestForceGeneParsing:
    """测试 HVG 模块中强制包含基因的字符串解析逻辑。"""

    @staticmethod
    def _parse_force_genes(s):
        """复现 hvg.py 第 70 行的解析逻辑。"""
        if not s.strip():
            return []
        return [g.strip() for g in s.replace('\n', ',').split(',') if g.strip()]

    def test_parse_gene_string(self):
        """逗号和换行符分隔的基因字符串应正确解析。"""
        result = self._parse_force_genes('Gene1,Gene2\nGene3')
        assert result == ['Gene1', 'Gene2', 'Gene3']

    def test_parse_empty_string(self):
        """空字符串应返回空列表。"""
        result = self._parse_force_genes('')
        assert result == []

    def test_parse_with_spaces(self):
        """带空格的基因名应被 trim。"""
        result = self._parse_force_genes(' Gene1 , Gene2 ')
        assert result == ['Gene1', 'Gene2']

    def test_parse_whitespace_only(self):
        """纯空白字符串应返回空列表。"""
        result = self._parse_force_genes('   \n  \t  ')
        assert result == []

    def test_parse_multiple_newlines(self):
        """多个换行符和逗号混合。"""
        result = self._parse_force_genes('A\nB,C\n\nD')
        assert result == ['A', 'B', 'C', 'D']

    def test_parse_trailing_comma(self):
        """尾部逗号不应产生空元素。"""
        result = self._parse_force_genes('Gene1,Gene2,')
        assert result == ['Gene1', 'Gene2']


# ============================================================
# 7. validate_input — batch_correct, clustering, trajectory
# ============================================================

def _make_adata_for_validation(n_obs=10, n_vars=5, obs_dict=None, obsm=None, uns=None):
    """创建用于输入校验测试的最小 AnnData。"""
    X = np.random.rand(n_obs, n_vars)
    obs = pd.DataFrame(obs_dict or {}, index=[f'cell_{i}' for i in range(n_obs)])
    adata = anndata.AnnData(X=X, obs=obs)
    adata.var_names = [f'Gene{i}' for i in range(n_vars)]
    if obsm:
        for k, v in obsm.items():
            adata.obsm[k] = v
    if uns:
        for k, v in uns.items():
            adata.uns[k] = v
    return adata


class _StubBaseForValidation:
    """最小化 BaseAnalysis 接口桩，用于实例化真实模块类调用 validate_input。"""

    def __init__(self, project_dir='/tmp/test', params=None, progress_callback=None):
        self.project_dir = project_dir
        self.params = params or {}
        self.progress_callback = progress_callback or (lambda pct, msg: None)


class TestBatchCorrectValidateInput:
    """测试 batch_correct 模块的输入校验。"""

    def _get_module(self):
        from modules.batch_correct import BatchCorrectAnalysis
        mod = BatchCorrectAnalysis(
            project_dir='/tmp/test',
            params={'batch_key': 'batch'},
            progress_callback=lambda pct, msg: None,
        )
        return mod

    def test_missing_pca_returns_error(self):
        """adata 无 X_pca → 返回错误字符串。"""
        mod = self._get_module()
        adata = _make_adata_for_validation()
        result = mod.validate_input(adata)
        assert isinstance(result, str)
        assert 'PCA' in result or 'pca' in result.lower()

    def test_with_pca_returns_none(self):
        """adata 有 X_pca → 返回 None（校验通过）。"""
        mod = self._get_module()
        adata = _make_adata_for_validation(obsm={'X_pca': np.random.rand(10, 10)})
        result = mod.validate_input(adata)
        assert result is None

    def test_corrected_requested_without_corrected_representation_is_rejected(self):
        """请求校正表示时不能静默退回原始 PCA。"""
        from modules.clustering import ClusteringAnalysis
        mod = ClusteringAnalysis(
            project_dir='/tmp/test',
            params={'resolutions': '0.8', 'use_corrected': True},
            progress_callback=lambda pct, msg: None,
        )
        adata = _make_adata_for_validation(obsm={'X_pca': np.random.rand(10, 5)})
        result = mod.validate_input(adata)
        assert isinstance(result, str)
        assert '不会静默回退' in result

    def test_corrected_requested_prefers_harmony_representation(self):
        """Harmony 表示存在时，校正请求解析为 X_pca_harmony。"""
        from modules.clustering import ClusteringAnalysis
        mod = ClusteringAnalysis(
            project_dir='/tmp/test',
            params={'resolutions': '0.8', 'use_corrected': True},
            progress_callback=lambda pct, msg: None,
        )
        adata = _make_adata_for_validation(obsm={
            'X_pca': np.random.rand(10, 5),
            'X_pca_harmony': np.random.rand(10, 5),
        })
        assert mod.validate_input(adata) is None
        assert mod._resolve_representation(adata, True) == ('X_pca_harmony', False)


class TestClusteringValidateInput:
    """测试 clustering 模块的输入校验。"""

    def _get_module(self):
        from modules.clustering import ClusteringAnalysis
        mod = ClusteringAnalysis(
            project_dir='/tmp/test',
            params={'resolutions': '0.8'},
            progress_callback=lambda pct, msg: None,
        )
        return mod

    def test_missing_umap_is_allowed_when_representation_exists(self):
        """聚类会从 PCA/校正表示重算 UMAP，不要求旧 X_umap。"""
        mod = self._get_module()
        adata = _make_adata_for_validation(obsm={
            'X_pca': np.random.rand(10, 5),
        })
        result = mod.validate_input(adata)
        assert result is None

    def test_with_umap_returns_none(self):
        """adata 有 UMAP 和 PCA 表示 → 返回 None（校验通过）。"""
        mod = self._get_module()
        adata = _make_adata_for_validation(obsm={
            'X_umap': np.random.rand(10, 2),
            'X_pca': np.random.rand(10, 5),
        })
        result = mod.validate_input(adata)
        assert result is None

    def test_missing_representation_returns_error(self):
        """仅有 UMAP 不能支撑 KNN/聚类。"""
        mod = self._get_module()
        adata = _make_adata_for_validation(obsm={'X_umap': np.random.rand(10, 2)})
        result = mod.validate_input(adata)
        assert isinstance(result, str)
        assert 'representation' in result.lower()


class TestTrajectoryValidateInput:
    """测试 trajectory 模块的输入校验。"""

    def _get_module(self):
        from modules.trajectory import TrajectoryAnalysis
        mod = TrajectoryAnalysis(
            project_dir='/tmp/test',
            params={},
            progress_callback=lambda pct, msg: None,
        )
        return mod

    def test_missing_neighbors_returns_error(self):
        """adata 无 neighbors → 返回错误字符串。"""
        mod = self._get_module()
        adata = _make_adata_for_validation()
        result = mod.validate_input(adata)
        assert isinstance(result, str)
        assert '邻居' in result or 'neighbor' in result.lower()

    def test_with_neighbors_returns_none(self):
        """adata 有 neighbors 且含 X_umap → 返回 None（校验通过）。"""
        mod = self._get_module()
        adata = _make_adata_for_validation(
            uns={'neighbors': {}},
            obsm={'X_umap': np.random.rand(10, 2)},
        )
        result = mod.validate_input(adata)
        assert result is None

    def test_with_neighbors_but_missing_umap_returns_error(self):
        """adata 有 neighbors 但无 X_umap → 返回错误字符串（防止下游 KeyError）。"""
        mod = self._get_module()
        adata = _make_adata_for_validation(uns={'neighbors': {}})
        result = mod.validate_input(adata)
        assert isinstance(result, str)
        assert 'UMAP' in result or 'umap' in result
