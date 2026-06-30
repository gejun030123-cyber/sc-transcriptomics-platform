import os
import numpy as np
from .base import BaseAnalysis
from .io_utils import convert_10x_to_h5ad


class Convert10x(BaseAnalysis):
    MODULE_NAME = "convert_10x"
    DISPLAY_NAME = "10x 数据转换"
    DESCRIPTION = "将 10x Genomics 三文件格式转换为 h5ad"

    def validate_input(self, adata):
        return None  # 不需要输入 adata

    def run(self, input_path):
        mtx_dir = self.params.get('mtx_dir', '')
        if not mtx_dir:
            return {
                'output_adata': input_path,
                'result_files': [],
                'summary': {},
                'error': '缺少 mtx_dir 参数，请提供 10x 数据目录路径',
            }
        species = self.params.get('species')
        genome = self.params.get('genome')

        self.progress(10, '正在读取 10x 数据...')
        output_path = os.path.join(self.project_dir, 'uploads', 'converted_10x.h5ad')

        self.progress(30, '正在解析矩阵文件...')
        adata = convert_10x_to_h5ad(mtx_dir, output_path, species=species, genome=genome)

        self.progress(90, '正在计算统计信息...')

        n_cells = adata.n_obs
        n_genes = adata.n_vars

        # 计算稀疏度
        total_elements = n_cells * n_genes
        if hasattr(adata.X, 'nnz'):
            # 稀疏矩阵
            nonzero = adata.X.nnz
        else:
            nonzero = np.count_nonzero(adata.X)
        sparsity = round((1 - nonzero / total_elements) * 100, 1) if total_elements > 0 else 0.0

        # 文件大小
        file_size_mb = round(os.path.getsize(output_path) / (1024 * 1024), 1)

        self.progress(100, f'转换完成，共 {n_cells} 个细胞，{n_genes} 个基因')

        return {
            'output_adata': output_path,
            'result_files': [],
            'summary': {
                'n_cells': n_cells,
                'n_genes': n_genes,
                'sparsity': sparsity,
                'file_size_mb': file_size_mb,
                'output_file': os.path.basename(output_path),
            }
        }
