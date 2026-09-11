# tests/test_signature_scoring.py
"""细胞类型签名评分测试"""
import pytest
import os
import numpy as np


class TestCellMarkers:
    """测试 cell_markers 模块."""

    def test_get_builtin_markers(self):
        """测试内置 marker 库查询."""
        from modules.cell_markers import get_markers

        result = get_markers('microglia')
        assert result['source'] == 'builtin'
        assert 'P2RY12' in result['positive']
        assert 'TMEM119' in result['positive']
        assert 'S100A8' in result['negative']

    def test_get_markers_unknown_type(self):
        """测试未知细胞类型返回 warning."""
        from modules.cell_markers import get_markers

        result = get_markers('unknown_cell_type')
        assert result['source'] == 'none'
        assert len(result['warnings']) > 0

    def test_user_markers_override(self):
        """测试用户自定义 marker 优先于内置库."""
        from modules.cell_markers import get_markers

        result = get_markers(
            'microglia',
            user_positive=['MYGENE1', 'MYGENE2'],
            user_negative=['BADGENE1'],
        )
        assert result['source'] == 'user'
        assert 'MYGENE1' in result['positive']
        assert 'BADGENE1' in result['negative']

    def test_mixed_markers(self):
        """测试部分自定义 marker（用户 positive + 内置 negative）."""
        from modules.cell_markers import get_markers

        result = get_markers(
            'microglia',
            user_positive=['MYGENE1'],
        )
        assert result['source'] == 'mixed'

    def test_marker_case_normalization(self):
        """测试 gene symbol 大小写处理."""
        from modules.cell_markers import get_markers

        result = get_markers(
            'microglia',
            user_positive=['p2ry12', 'tmem119'],
            user_negative=['S100a8'],
        )
        assert 'P2RY12' in result['positive']
        assert 'TMEM119' in result['positive']
        assert 'S100A8' in result['negative']

    def test_list_builtin_cell_types(self):
        """测试列出所有内置细胞类型."""
        from modules.cell_markers import list_builtin_cell_types

        types = list_builtin_cell_types()
        assert 'microglia' in types
        assert 't_cell' in types
        assert 'b_cell' in types
        assert 'macrophage' in types
        assert types['microglia']['display'] == '小胶质细胞'


class TestSignatureScoring:
    """测试 sc_cluster 评分器（确定性算法）."""

    @pytest.fixture
    def mini_adata(self, tmp_path):
        """创建最小 AnnData 用于测试."""
        try:
            import scanpy as sc
            import anndata
        except ImportError:
            pytest.skip("scanpy/anndata 未安装")

        np.random.seed(42)
        n_cells = 100
        n_genes = 50

        # 创建表达矩阵
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)

        # gene names: 前5个是 marker 基因
        gene_names = ['P2RY12', 'TMEM119', 'CX3CR1', 'AIF1',
                       'S100A8', 'FCGR3A', 'CD3D', 'MS4A1'] + \
                      [f'GENE_{i}' for i in range(8, n_genes)]

        adata = anndata.AnnData(X)
        adata.var_names = gene_names

        # 创建两个 cluster: cluster 0 有 microglia-like 表达，cluster 1 没有
        adata.obs['leiden'] = ['0'] * 50 + ['1'] * 50

        # cluster 0 中高表达 P2RY12, TMEM119, CX3CR1, AIF1
        for g in ['P2RY12', 'TMEM119', 'CX3CR1', 'AIF1']:
            g_idx = list(adata.var_names).index(g)
            adata.X[0:50, g_idx] = np.random.negative_binomial(20, 0.3, size=50).astype(np.float32)
            adata.X[50:100, g_idx] = np.random.negative_binomial(3, 0.3, size=50).astype(np.float32)

        # cluster 0 中低表达 S100A8
        g_idx = list(adata.var_names).index('S100A8')
        adata.X[0:50, g_idx] = np.random.negative_binomial(2, 0.3, size=50).astype(np.float32)
        adata.X[50:100, g_idx] = np.random.negative_binomial(10, 0.3, size=50).astype(np.float32)

        # 添加 batch 列
        adata.obs['batch'] = ['batch_1'] * 30 + ['batch_2'] * 20 + ['batch_1'] * 25 + ['batch_2'] * 25

        # 添加 UMAP
        adata.obsm['X_umap'] = np.random.randn(n_cells, 2).astype(np.float32)

        # 保存
        path = str(tmp_path / 'test.h5ad')
        adata.write(path)
        return path

    def test_score_cluster_signature_identifies_target(self, mini_adata):
        """测试评分器正确识别 microglia-like cluster."""
        from modules.evaluators.sc_cluster import score_cluster_signature

        result = score_cluster_signature(
            adata_path=mini_adata,
            cluster_key='leiden',
            positive_markers=['P2RY12', 'TMEM119', 'CX3CR1', 'AIF1'],
            negative_markers=['S100A8', 'FCGR3A'],
        )

        assert 'error' not in result
        assert len(result['cluster_scores']) == 2

        # cluster 0 应该得分更高
        best = result['cluster_scores'][0]
        assert best['cluster'] == '0'
        assert best['total_score'] > 0  # 核心断言：评分器能区分

    def test_missing_markers_return_warning(self, mini_adata):
        """测试 marker 不存在时返回 warning 而非崩溃."""
        from modules.evaluators.sc_cluster import score_cluster_signature

        result = score_cluster_signature(
            adata_path=mini_adata,
            cluster_key='leiden',
            positive_markers=['NONEXISTENT_GENE'],
            negative_markers=[],
        )

        assert len(result['warnings']) > 0
        assert 'NONEXISTENT_GENE' in str(result['warnings'])

    def test_invalid_adata_path_rejected(self):
        """测试非法路径被拒绝."""
        from modules.evaluators.sc_cluster import score_cluster_signature

        result = score_cluster_signature(
            adata_path='/nonexistent/path.h5ad',
            cluster_key='leiden',
            positive_markers=['GENE'],
            negative_markers=[],
        )

        assert 'error' in result

    def test_score_range_zero_to_one(self, mini_adata):
        """测试 total_score 始终在 [0, 1] 范围内."""
        from modules.evaluators.sc_cluster import score_cluster_signature

        result = score_cluster_signature(
            adata_path=mini_adata,
            cluster_key='leiden',
            positive_markers=['P2RY12', 'TMEM119'],
            negative_markers=['S100A8'],
        )

        for cs in result['cluster_scores']:
            assert 0.0 <= cs['total_score'] <= 1.0

    def test_confidence_levels(self, mini_adata):
        """测试 score_cell_type_signature 返回合理的置信度."""
        from modules.evaluators.sc_cluster import score_cell_type_signature

        result = score_cell_type_signature(
            adata_path=mini_adata,
            cluster_key='leiden',
            target_cell_type='microglia',
        )

        assert 'best_cluster' in result
        assert 'confidence' in result
        assert result['confidence'] in ('high', 'medium', 'low', 'very_low', 'none')

    def test_deterministic_scoring(self, mini_adata):
        """测试评分结果是确定性的（相同输入 → 相同输出）."""
        from modules.evaluators.sc_cluster import score_cluster_signature

        result1 = score_cluster_signature(
            adata_path=mini_adata,
            cluster_key='leiden',
            positive_markers=['P2RY12', 'TMEM119'],
            negative_markers=['S100A8'],
        )

        result2 = score_cluster_signature(
            adata_path=mini_adata,
            cluster_key='leiden',
            positive_markers=['P2RY12', 'TMEM119'],
            negative_markers=['S100A8'],
        )

        for cs1, cs2 in zip(result1['cluster_scores'], result2['cluster_scores']):
            assert cs1['total_score'] == cs2['total_score']


class TestSignatureScoringExpressionScale:
    """评分必须与表达量尺度无关，并拒绝无法比较的矩阵。"""

    @staticmethod
    def _write(adata, path):
        import anndata
        adata.write_h5ad(path)

    def _two_cluster_adata(self, counts):
        import anndata
        import pandas as pd

        adata = anndata.AnnData(X=counts.astype('float32'))
        adata.var_names = ['P2RY12', 'TMEM119', 'S100A8', 'ACTB']
        adata.obs['leiden'] = pd.Categorical(['0'] * 20 + ['1'] * 20)
        return adata

    def test_raw_counts_are_normalized_before_scoring(self, tmp_path):
        """原始 counts 输入不再因为量纲而把 total_score 顶到饱和值."""
        import numpy as np
        from modules.evaluators.sc_cluster import score_cluster_signature

        counts = np.full((40, 4), 1.0)
        # cluster 0 高表达 microglia marker，cluster 1 高表达 S100A8
        counts[0:20, 0] = 400.0
        counts[0:20, 1] = 300.0
        counts[20:40, 2] = 500.0
        adata = self._two_cluster_adata(counts)
        path = tmp_path / 'counts.h5ad'
        self._write(adata, path)

        result = score_cluster_signature(
            adata_path=str(path),
            cluster_key='leiden',
            positive_markers=['P2RY12', 'TMEM119'],
            negative_markers=['S100A8'],
        )

        assert 'error' not in result
        assert result['expression_scale'] == 'log1p_normalized_from_X_counts'
        best = result['cluster_scores'][0]
        assert best['cluster'] == '0'
        assert best['total_score'] < 1.0  # 旧的 0.35*原始均值 会直接饱和到 1.0
        assert best['positive_evidence'] == pytest.approx(1.0)

    def test_scoring_is_invariant_to_input_expression_scale(self, tmp_path):
        """同一份数据的 counts 与 log1p 输入必须给出同一结论与同一分数."""
        import numpy as np
        import scanpy as sc
        from modules.evaluators.sc_cluster import score_cluster_signature

        counts = np.full((40, 4), 1.0)
        counts[0:20, 0] = 400.0
        counts[0:20, 1] = 300.0
        counts[20:40, 2] = 500.0
        counts_path = tmp_path / 'counts.h5ad'
        self._write(self._two_cluster_adata(counts), counts_path)

        normalized = self._two_cluster_adata(counts)
        sc.pp.normalize_total(normalized, target_sum=1e4)
        sc.pp.log1p(normalized)
        normalized_path = tmp_path / 'normalized.h5ad'
        self._write(normalized, normalized_path)

        scores = {}
        for label, path in (('counts', counts_path), ('log1p', normalized_path)):
            result = score_cluster_signature(
                adata_path=str(path),
                cluster_key='leiden',
                positive_markers=['P2RY12', 'TMEM119'],
                negative_markers=['S100A8'],
            )
            assert 'error' not in result
            scores[label] = result['cluster_scores'][0]

        assert scores['counts']['cluster'] == scores['log1p']['cluster'] == '0'
        assert scores['counts']['total_score'] == pytest.approx(
            scores['log1p']['total_score'], abs=0.05,
        )

    def test_scaled_matrix_is_rejected_instead_of_scored(self, tmp_path):
        """已 scale（含负值）且无 counts 层时明确报错，而不是给出虚假高分."""
        import numpy as np
        from modules.evaluators.sc_cluster import score_cell_type_signature

        rng = np.random.default_rng(0)
        scaled = rng.normal(size=(40, 4)).astype('float32')
        adata = self._two_cluster_adata(scaled)
        path = tmp_path / 'scaled.h5ad'
        self._write(adata, path)

        result = score_cell_type_signature(
            adata_path=str(path),
            cluster_key='leiden',
            target_cell_type='microglia',
        )

        assert 'error' in result
        assert result['confidence'] == 'none'
        assert result['cluster_scores'] == []

    def test_batch_and_size_do_not_change_marker_confidence(self, tmp_path):
        """簇大小/单批次只作为 qc_flags，不再进入 marker 得分."""
        import numpy as np
        import pandas as pd
        from modules.evaluators.sc_cluster import score_cluster_signature

        counts = np.full((40, 4), 1.0)
        counts[0:20, 0] = 400.0
        counts[0:20, 1] = 300.0
        adata = self._two_cluster_adata(counts)
        adata.obs['batch'] = pd.Categorical(['b1'] * 20 + ['b2'] * 20)
        path = tmp_path / 'batch.h5ad'
        self._write(adata, path)

        result = score_cluster_signature(
            adata_path=str(path),
            cluster_key='leiden',
            positive_markers=['P2RY12', 'TMEM119'],
            negative_markers=[],
        )

        for score in result['cluster_scores']:
            assert 0.0 <= score['total_score'] <= 1.0
            assert isinstance(score['qc_flags'], list)
        # 每个簇都只来自一个批次：罚分只出现在 qc_flags / batch_penalty
        assert result['cluster_scores'][0]['batch_penalty'] == pytest.approx(0.5)
        assert 'single_batch_only' in result['cluster_scores'][0]['qc_flags']
        marker_only = (
            0.40 * result['cluster_scores'][0]['positive_evidence']
            + 0.25 * result['cluster_scores'][0]['detection_fraction']
            + 0.35 * result['cluster_scores'][0]['specificity_score']
            - 0.30 * result['cluster_scores'][0]['negative_ratio']
        )
        assert result['cluster_scores'][0]['total_score'] == pytest.approx(
            max(0.0, min(1.0, marker_only)), abs=5e-4,
        )
