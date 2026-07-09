import os
import re
from .base import BaseAnalysis
from .io_utils import (
    convert_10x_to_h5ad,
    infer_sc_data_format,
    summarize_adata_import,
    write_single_cell_h5ad,
)


class Convert10x(BaseAnalysis):
    MODULE_NAME = "convert_10x"
    DISPLAY_NAME = "单细胞数据导入"
    DESCRIPTION = "将 h5ad / 10x mtx / 10x h5 / loom / zarr 转换为标准 h5ad"

    def validate_input(self, adata):
        return None  # 不需要输入 adata

    def run(self, input_path):
        source_path = str(self.params.get('source_path') or '').strip()
        mtx_dir = str(self.params.get('mtx_dir') or '').strip()
        if not source_path and not mtx_dir:
            return {
                'output_adata': input_path,
                'result_files': [],
                'summary': {},
                'error': '缺少 source_path 或 mtx_dir 参数，请提供单细胞数据路径',
            }
        species = self.params.get('species')
        genome = self.params.get('genome')
        input_format = self.params.get('input_format', 'auto') or 'auto'

        uploads_dir = os.path.join(self.project_dir, 'uploads')
        os.makedirs(uploads_dir, exist_ok=True)

        self.progress(10, '正在准备导入单细胞数据...')

        if not source_path:
            detected_format = '10x_mtx'
            output_path = os.path.join(self.project_dir, 'uploads', 'converted_10x.h5ad')

            self.progress(30, '正在解析矩阵文件...')
            adata = convert_10x_to_h5ad(mtx_dir, output_path, species=species, genome=genome)
        else:
            source_path = os.path.abspath(source_path)
            detected_format = infer_sc_data_format(source_path) if input_format == 'auto' else input_format
            stem = os.path.basename(source_path.rstrip(os.sep))
            stem = re.sub(r'\.(h5ad|h5|hdf5|loom|zarr|csv|tsv|txt|xlsx|xls|mtx)(\.gz)?$', '', stem, flags=re.I)
            stem = re.sub(r'[^A-Za-z0-9_.-]+', '_', stem).strip('._-') or 'single_cell'
            output_path = os.path.join(uploads_dir, f'{stem}_imported.h5ad')

            self.progress(35, '正在标准化 AnnData 结构...')
            adata = write_single_cell_h5ad(
                source_path,
                output_path,
                input_format=input_format,
                species=species,
                genome=genome,
            )

        self.progress(90, '正在计算统计信息...')
        summary = summarize_adata_import(adata, detected_format, output_path)

        self.progress(100, f"转换完成，共 {summary['n_cells']} 个细胞，{summary['n_genes']} 个基因")

        return {
            'output_adata': output_path,
            'result_files': [
                {
                    'file_path': output_path,
                    'file_type': 'h5ad',
                    'category': 'data',
                    'label': 'Imported single-cell h5ad',
                }
            ],
            'summary': summary,
        }
