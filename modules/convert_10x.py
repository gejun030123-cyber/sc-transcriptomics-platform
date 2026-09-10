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
            # mtx_dir 与单文件导入一样必须限制在项目目录内，防止参数注入
            # 读取服务器任意目录。
            project_root = os.path.realpath(self.project_dir)
            real_mtx_dir = os.path.realpath(mtx_dir)
            if not (real_mtx_dir == project_root
                    or real_mtx_dir.startswith(project_root + os.sep)):
                raise ValueError('mtx_dir 必须位于项目目录内: ' + mtx_dir)
            if not os.path.isdir(real_mtx_dir):
                raise FileNotFoundError('10x 矩阵目录不存在: ' + mtx_dir)
            detected_format = '10x_mtx'
            output_path = os.path.join(self.project_dir, 'uploads', 'converted_10x.h5ad')

            self.progress(30, '正在解析矩阵文件...')
            adata = convert_10x_to_h5ad(real_mtx_dir, output_path, species=species, genome=genome)
            if 'sample_id' not in adata.obs.columns:
                sample_id = _safe_batch_name(
                    os.path.basename(real_mtx_dir.rstrip(os.sep)), 'sample_1')
                adata.obs['sample_id'] = sample_id
        else:
            source_path = os.path.abspath(source_path)
            # 单文件导入与批次导入必须一致地限制在项目目录内，防止通过
            # 参数注入读取服务器任意文件。
            project_root = os.path.realpath(self.project_dir)
            real_source = os.path.realpath(source_path)
            if not (real_source == project_root or real_source.startswith(project_root + os.sep)):
                raise ValueError('导入文件必须位于项目目录内: ' + source_path)
            if not os.path.isfile(real_source):
                raise FileNotFoundError('导入文件不存在: ' + source_path)
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
            if 'sample_id' not in adata.obs.columns:
                adata.obs['sample_id'] = _safe_batch_name(stem, 'single_cell')
                # ``write_single_cell_h5ad`` has already written the
                # standardized object, so persist the generated sample ID as
                # well.  Downstream pseudobulk/time-course modules rely on
                # this column being present in the actual h5ad file.
                adata.write_h5ad(output_path)

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
        sample_ids = []
        conditions = []
        batch_records = []
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
            sample_id = _safe_batch_name(item.get('sample_id'), batch_name)
            condition = str(item.get('condition') or '').strip()
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
            adata.obs['sample_id'] = sample_id
            adata.obs['condition'] = condition
            adatas.append(adata)
            batch_names.append(batch_name)
            sample_ids.append(sample_id)
            conditions.append(condition)
            total_counts = adata.obs['total_counts'] if 'total_counts' in adata.obs else None
            batch_records.append({
                'batch': batch_name,
                'sample_id': sample_id,
                'condition': condition,
                'n_cells': int(adata.n_obs),
                'n_genes': int(adata.n_vars),
                'median_total_counts': round(float(total_counts.median()), 3) if total_counts is not None else None,
            })

        if len(set(sample_ids)) != len(sample_ids):
            raise ValueError('样本 ID（sample_id）必须唯一且不能为空')

        self.progress(75, '正在合并批次...')
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
        combined.obs['sample_id'] = combined.obs['sample_id'].astype(str)
        combined.obs['condition'] = combined.obs['condition'].astype(str)
        combined.uns['batch_imports'] = batch_names
        combined.uns['input_format'] = '10x_mtx_zip_batches'
        combined.uns['sc_batch_import'] = {
            'sample_key': 'sample_id',
            'condition_key': 'condition',
            'source': 'upload_marked',
            'input_format': '10x_mtx_zip_batches',
        }
        output_path = os.path.join(uploads_dir, 'combined_batches_imported.h5ad')
        combined.write_h5ad(output_path)
        import pandas as pd
        batch_summary = pd.DataFrame(batch_records)
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        summary_csv = os.path.join(results_dir, 'batch_import_summary.csv')
        batch_summary.to_csv(summary_csv, index=False)
        result_files = [{
            'file_path': output_path,
            'file_type': 'h5ad',
            'category': 'data',
            'label': 'Combined multi-batch single-cell h5ad',
        }, {
            'file_path': summary_csv,
            'file_type': 'csv',
            'category': 'table',
            'label': 'Batch import summary',
        }]
        try:
            import matplotlib.pyplot as plt
            plots_dir = os.path.join(self.project_dir, 'plots')
            os.makedirs(plots_dir, exist_ok=True)
            fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
            palette = ['#3C5488', '#E64B35', '#00A087', '#4DBBD5']
            bars = ax.bar(batch_summary['batch'], batch_summary['n_cells'],
                          color=[palette[i % len(palette)] for i in range(len(batch_summary))],
                          width=0.62)
            ax.set_title('Cells per imported batch', fontsize=15, fontweight='semibold')
            ax.set_ylabel('Cells')
            ax.grid(axis='y', color='#d8dee9', linewidth=0.6, alpha=0.7)
            ax.set_axisbelow(True)
            for spine in ax.spines.values():
                spine.set_visible(False)
            for bar in bars:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f'{int(bar.get_height()):,}', ha='center', va='bottom', fontsize=10)
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir,
                'batch_import_cell_counts.png', 'batch', 'Batch cell counts',
            ))
            plt.close(fig)
        except Exception as exc:
            self.progress(-1, f'批次组成图生成失败（不影响合并结果）：{exc}')
        summary = summarize_adata_import(combined, '10x_mtx_zip_batches', output_path)
        summary.update({
            'n_batches': len(batch_names),
            'batch_names': batch_names,
            'sample_ids': sample_ids,
            'conditions': conditions,
            'batch_column': 'batch',
            'sample_column': 'sample_id',
            'condition_column': 'condition',
        })
        self.progress(100, f'批次导入完成，共 {summary["n_cells"]} 个细胞、{summary["n_genes"]} 个基因')
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
