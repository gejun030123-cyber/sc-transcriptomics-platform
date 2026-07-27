# tests/test_io_utils.py
"""Tests for modules/io_utils.py — remap_var_names and read_expression_matrix."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import pytest
import anndata
from scipy import sparse
from scipy.io import mmwrite


class TestRead10xMtxCompat:
    """10x Matrix Market import accepts all supported compression combinations."""

    def test_mixed_compression_preserves_all_three_files(self, tmp_path):
        from modules.io_utils import convert_10x_to_h5ad, read_single_cell_data
        import gzip

        matrix_dir = tmp_path / 'mixed_10x'
        matrix_dir.mkdir()
        mmwrite(matrix_dir / 'matrix.mtx', sparse.coo_matrix([[1, 0], [0, 2], [3, 0]]))
        (matrix_dir / 'barcodes.tsv').write_text('cell_a\ncell_b\n', encoding='utf-8')
        # v2 uses genes.tsv (two columns), which recent Scanpy still expects
        # to receive under the features.tsv.gz filename.
        (matrix_dir / 'genes.tsv').write_text(
            'ENSG000001\tGeneA\n'
            'ENSG000002\tGeneB\n'
            'ENSG000003\tGeneC\n', encoding='utf-8')

        # This is a common mixed set: matrix compressed, metadata files plain.
        with open(matrix_dir / 'matrix.mtx', 'rb') as src, gzip.open(
            matrix_dir / 'matrix.mtx.gz', 'wb'
        ) as dst:
            dst.write(src.read())
        (matrix_dir / 'matrix.mtx').unlink()

        output_path = tmp_path / 'converted.h5ad'
        converted = convert_10x_to_h5ad(str(matrix_dir), str(output_path))
        imported = read_single_cell_data(str(matrix_dir), input_format='10x_mtx')

        assert output_path.exists()
        assert converted.shape == (2, 3)
        assert imported.shape == (2, 3)
        assert 'counts' in imported.layers


# --- remap_var_names tests ---

class TestRemapVarNames:
    """测试 remap_var_names 函数。"""

    def test_already_gene_names(self):
        """var_names 已是基因名（非 Ensembl ID）→ 不做任何修改。"""
        from modules.io_utils import remap_var_names
        adata = anndata.AnnData(
            X=np.random.rand(3, 4),
            var=pd.DataFrame(index=['TP53', 'BRCA1', 'MYC', 'EGFR'])
        )
        result = remap_var_names(adata)
        assert list(result.var_names) == ['TP53', 'BRCA1', 'MYC', 'EGFR']
        assert 'gene_id' not in result.var.columns

    def test_ensembl_with_gene_name_column(self):
        """Ensembl ID + var 中有 gene_name 列 → 映射为基因名。"""
        from modules.io_utils import remap_var_names
        adata = anndata.AnnData(
            X=np.random.rand(3, 4),
            var=pd.DataFrame({
                'gene_name': ['TP53', 'BRCA1', 'MYC', 'EGFR']
            }, index=['ENSG00000141510', 'ENSG00000012048', 'ENSG00000136997', 'ENSG00000146648'])
        )
        result = remap_var_names(adata)
        assert list(result.var_names) == ['TP53', 'BRCA1', 'MYC', 'EGFR']
        assert 'gene_id' in result.var.columns
        assert result.var['gene_id'].iloc[0] == 'ENSG00000141510'

    def test_duplicate_gene_names(self):
        """重复基因名 → 第一次保留原名，后续加 _1, _2 后缀。"""
        from modules.io_utils import remap_var_names
        adata = anndata.AnnData(
            X=np.random.rand(3, 4),
            var=pd.DataFrame({
                'gene_name': ['TP53', 'TP53', 'MYC', 'TP53']
            }, index=['ENSG00000141510', 'ENSG00000141511', 'ENSG00000136997', 'ENSG00000141512'])
        )
        result = remap_var_names(adata)
        names = list(result.var_names)
        assert names[0] == 'TP53'
        assert names[1] == 'TP53_1'
        assert names[2] == 'MYC'
        assert names[3] == 'TP53_2'

    def test_no_gene_name_column(self):
        """Ensembl ID 但 var 中无 gene_name 列 → 保持原样。"""
        from modules.io_utils import remap_var_names
        adata = anndata.AnnData(
            X=np.random.rand(2, 3),
            var=pd.DataFrame(index=['ENSG00000141510', 'ENSG00000012048', 'ENSG00000136997'])
        )
        result = remap_var_names(adata)
        assert list(result.var_names) == ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000136997']

    def test_nan_in_gene_names(self):
        """gene_name 列含 nan 值 → 不崩溃，nan 位置保留原始 ID。"""
        from modules.io_utils import remap_var_names
        adata = anndata.AnnData(
            X=np.random.rand(2, 3),
            var=pd.DataFrame({
                'gene_name': ['TP53', np.nan, 'MYC']
            }, index=['ENSG00000141510', 'ENSG00000012048', 'ENSG00000136997'])
        )
        result = remap_var_names(adata)
        names = list(result.var_names)
        assert names[0] == 'TP53'
        assert names[2] == 'MYC'
        # nan 位置应保留原始 Ensembl ID
        assert 'ENSG00000012048' in names[1]

    def test_empty_adata(self):
        """空 AnnData → 不崩溃。"""
        from modules.io_utils import remap_var_names
        adata = anndata.AnnData(X=np.empty((0, 0)))
        result = remap_var_names(adata)
        assert result.n_vars == 0


# --- read_expression_matrix tests ---

class TestReadExpressionMatrix:
    """测试 read_expression_matrix 函数。"""

    def test_csv_file(self, tmp_path):
        """读取标准 CSV 文件 → 返回正确 AnnData。"""
        from modules.io_utils import read_expression_matrix
        data = {'SampleA': [10.0, 20.0, 30.0], 'SampleB': [40.0, 50.0, 60.0]}
        df = pd.DataFrame(data, index=['ENSG00000141510', 'ENSG00000012048', 'ENSG00000136997'])
        csv_path = str(tmp_path / 'test.csv')
        df.to_csv(csv_path)

        adata = read_expression_matrix(csv_path)
        assert adata.n_obs == 2
        assert adata.n_vars == 3

    def test_tsv_file(self, tmp_path):
        """读取 TSV 文件 → 正确解析 tab 分隔数据。"""
        from modules.io_utils import read_expression_matrix
        df = pd.DataFrame({
            'Ctrl1': [1.0, 2.0],
            'Ctrl2': [3.0, 4.0],
        }, index=['GeneA', 'GeneB'])
        tsv_path = str(tmp_path / 'test.tsv')
        df.to_csv(tsv_path, sep='\t')

        adata = read_expression_matrix(tsv_path)
        assert adata.n_obs == 2
        assert adata.n_vars == 2

    def test_count_fpkm_suffix_cleanup(self, tmp_path):
        """样本名含 _count/_FPKM 后缀 → 清理后保留 count 数据。"""
        from modules.io_utils import read_expression_matrix
        data = {'SampleA_count': [10.0, 20.0], 'SampleA_FPKM': [1.5, 2.5],
                'SampleB_count': [30.0, 40.0], 'SampleB_FPKM': [3.5, 4.5]}
        df = pd.DataFrame(data, index=['ENSG1', 'ENSG2'])
        csv_path = str(tmp_path / 'test_mix.csv')
        df.to_csv(csv_path)

        adata = read_expression_matrix(csv_path)
        assert adata.n_obs == 2
        obs_names = list(adata.obs_names)
        assert all('_FPKM' not in n for n in obs_names)

    def test_count_fpkm_suffix_cleanup_tsv(self, tmp_path):
        """TSV 格式：样本名含 _count/_FPKM 后缀 → 清理后保留 count 数据。"""
        from modules.io_utils import read_expression_matrix
        data = {'SampleA_count': [10.0, 20.0], 'SampleA_FPKM': [1.5, 2.5],
                'SampleB_count': [30.0, 40.0], 'SampleB_FPKM': [3.5, 4.5]}
        df = pd.DataFrame(data, index=['ENSG1', 'ENSG2'])
        tsv_path = str(tmp_path / 'test_mix.tsv')
        df.to_csv(tsv_path, sep='\t')

        adata = read_expression_matrix(tsv_path)
        assert adata.n_obs == 2

    def test_numeric_only_filter_csv(self, tmp_path):
        """CSV 格式：文件含非数值列 → 只保留数值列。"""
        from modules.io_utils import read_expression_matrix
        data = {'SampleA': [10.0, 20.0, 30.0], 'SampleB': [40.0, 50.0, 60.0]}
        df = pd.DataFrame(data, index=['ENSG1', 'ENSG2', 'ENSG3'])
        df['GeneSymbol'] = ['TP53', 'BRCA1', 'MYC']
        csv_path = str(tmp_path / 'test_mixed.csv')
        df.to_csv(csv_path)
        adata = read_expression_matrix(csv_path)
        assert adata.n_obs == 2

    def test_numeric_only_filter_tsv(self, tmp_path):
        """TSV 格式：文件含非数值列 → 只保留数值列，基因名存入 var。"""
        from modules.io_utils import read_expression_matrix
        data = {'SampleA': [10.0, 20.0, 30.0], 'SampleB': [40.0, 50.0, 60.0]}
        df = pd.DataFrame(data, index=['ENSG1', 'ENSG2', 'ENSG3'])
        df['GeneSymbol'] = ['TP53', 'BRCA1', 'MYC']
        tsv_path = str(tmp_path / 'test_mixed.tsv')
        df.to_csv(tsv_path, sep='\t')

        adata = read_expression_matrix(tsv_path)
        assert adata.n_obs == 2
        assert adata.n_vars == 3
        if 'gene_name' in adata.var.columns:
            assert list(adata.var['gene_name']) == ['TP53', 'BRCA1', 'MYC']
