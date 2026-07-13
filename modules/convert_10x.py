import os
import re
import shutil
import zipfile
from .base import BaseAnalysis
from .io_utils import (
    convert_10x_to_h5ad,
    infer_sc_data_format,
    summarize_adata_import,
    write_single_cell_h5ad,
)


def _safe_extract_zip(zip_path, extract_dir):
    """Extract a user ZIP without allowing traversal or symlink entries."""
    root = os.path.realpath(extract_dir)
    os.makedirs(root, exist_ok=True)
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            name = member.filename.replace('\\', '/')
            if not name or name.startswith('/') or '..' in name.split('/'):
                raise ValueError(f'ZIP 包含不安全路径: {member.filename}')
            # Unix symlink bit; symlinks are not needed for 10x matrices.
            if ((member.external_attr >> 16) & 0o170000) == 0o120000:
                raise ValueError(f'ZIP 不允许包含符号链接: {member.filename}')
            target = os.path.realpath(os.path.join(root, name))
            if target != root and not target.startswith(root + os.sep):
                raise ValueError(f'ZIP 包含越界路径: {member.filename}')
        archive.extractall(root)
    return root


def _find_10x_matrix_dir(root):
    """Find one 10x matrix directory inside a ZIP extraction tree."""
    gene_names = {'features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz'}
    candidates = []
    for current, dirs, files in os.walk(root):
        names = set(files)
        if not ({'matrix.mtx', 'matrix.mtx.gz'} & names):
            continue
        if not ({'barcodes.tsv', 'barcodes.tsv.gz'} & names):
            continue
        if not (gene_names & names):
            continue
        candidates.append(current)
    if not candidates:
        raise ValueError('ZIP 内未找到完整的 10x matrix.mtx + barcodes + genes/features 文件')
    preferred = [path for path in candidates if os.path.basename(path).lower() in {
        'filtered_feature_bc_matrix', 'filtered_gene_bc_matrices', 'filtered_matrix'
    }]
    if len(preferred) == 1:
        return preferred[0]
    if len(candidates) == 1:
        return candidates[0]
    raise ValueError('ZIP 内找到多个 10x 矩阵，无法自动判断；请只保留一套 filtered_feature_bc_matrix')


def _safe_batch_name(value, fallback):
    name = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value or '')).strip('._-')
    return name or fallback


class Convert10x(BaseAnalysis):
    MODULE_NAME = "convert_10x"
    DISPLAY_NAME = "单细胞数据导入"
    DESCRIPTION = "将 h5ad / 10x mtx / 10x h5 / loom / zarr 转换为标准 h5ad"

    def validate_input(self, adata):
        return None  # 不需要输入 adata

    def run(self, input_path):
        batch_sources = self.params.get('batch_sources') or []
        if batch_sources:
            return self._run_batch_import(batch_sources)

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

    def _run_batch_import(self, batch_sources):
        """Import two (or more) zipped 10x batches and concatenate with labels."""
        import anndata as ad

        if len(batch_sources) < 2:
            raise ValueError('批次导入至少需要两个 ZIP 文件')
        if len({str(item.get('batch_name', '')).strip() for item in batch_sources}) != len(batch_sources):
            raise ValueError('批次名称必须唯一且不能为空')

        uploads_dir = os.path.join(self.project_dir, 'uploads')
        extraction_root = os.path.join(uploads_dir, 'batch_imports')
        os.makedirs(extraction_root, exist_ok=True)
        adatas = []
        batch_names = []
        species = self.params.get('species')
        genome = self.params.get('genome')

        for index, item in enumerate(batch_sources, start=1):
            zip_path = os.path.abspath(str(item.get('zip_path') or ''))
            project_root = os.path.realpath(self.project_dir)
            if not zip_path or not os.path.isfile(zip_path):
                raise FileNotFoundError(f'批次 ZIP 不存在: {zip_path}')
            if not os.path.realpath(zip_path).startswith(project_root + os.sep):
                raise ValueError('批次 ZIP 必须位于项目目录内')
            batch_name = _safe_batch_name(item.get('batch_name'), f'batch_{index}')
            extract_dir = os.path.join(extraction_root, f'{index}_{batch_name}')
            if os.path.isdir(extract_dir):
                shutil.rmtree(extract_dir)
            self.progress(10 + int((index - 1) * 30 / len(batch_sources)),
                          f'正在解压批次 {batch_name}...')
            _safe_extract_zip(zip_path, extract_dir)
            matrix_dir = _find_10x_matrix_dir(extract_dir)
            self.progress(20 + int((index - 1) * 30 / len(batch_sources)),
                          f'正在读取批次 {batch_name}...')
            adata = convert_10x_to_h5ad(
                matrix_dir,
                os.path.join(extract_dir, 'batch_imported.h5ad'),
                species=species,
                genome=genome,
            )
            adata.obs_names = [f'{batch_name}_{name}' for name in adata.obs_names]
            adata.obs['batch'] = batch_name
            adatas.append(adata)
            batch_names.append(batch_name)

        self.progress(75, '正在合并两个批次...')
        combined = ad.concat(
            adatas,
            axis=0,
            join='outer',
            label='batch_source',
            keys=batch_names,
            index_unique=None,
            fill_value=0,
        )
        combined.obs['batch'] = combined.obs['batch'].astype(str)
        combined.uns['batch_imports'] = batch_names
        combined.uns['input_format'] = '10x_mtx_zip_batches'
        output_path = os.path.join(uploads_dir, 'combined_batches_imported.h5ad')
        combined.write_h5ad(output_path)
        summary = summarize_adata_import(combined, '10x_mtx_zip_batches', output_path)
        summary.update({
            'n_batches': len(batch_names),
            'batch_names': batch_names,
            'batch_column': 'batch',
        })
        self.progress(100, f'批次导入完成，共 {summary["n_cells"]} 个细胞、{summary["n_genes"]} 个基因')
        return {
            'output_adata': output_path,
            'result_files': [{
                'file_path': output_path,
                'file_type': 'h5ad',
                'category': 'data',
                'label': 'Combined multi-batch single-cell h5ad',
            }],
            'summary': summary,
        }
