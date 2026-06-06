---
title: 10x 三文件转换功能修复与增强
date: 2026-06-05
status: approved
---

# 10x Genomics 三文件转 h5ad 功能修复与增强

## 背景

平台已有基于 `scanpy.read_10x_mtx()` 的 10x 三文件（barcodes.tsv、genes.tsv、matrix.mtx）转 h5ad 功能，但存在 `.tsv` 文件无法上传的问题，且转换完成后缺少数据预览信息。

## 改动范围

### 改动 1：修复 .tsv 文件上传

**文件：`routes/upload.py`**

- 在 `ALLOWED_EXT` 集合（第 9 行）中加入 `.tsv`

**文件：`templates/upload.html`**

- 第 17 行 `<input>` 的 `accept` 属性中加入 `.tsv`
- 第 16 行支持格式说明文字中加入 `.tsv`

**影响分析：** `.tsv` 是通用格式，加入后同时支持 10x 文件和其他 TSV 数据文件。`read_expression_matrix()` 已能处理 TSV，无兼容性问题。

### 改动 2：转换结果预览

**文件：`modules/convert_10x.py`**

在 `run()` 方法（第 14-35 行）中，在 `convert_10x_to_h5ad()` 调用之后、return 之前，计算并添加以下字段到 `summary`：

- `sparsity`: 数据稀疏度百分比（零值占比）
  - 对于稀疏矩阵：`(1 - adata.X.nnz / (n_cells * n_genes)) * 100`
  - 对于密集矩阵：`(1 - np.count_nonzero(adata.X) / (n_cells * n_genes)) * 100`
- `file_size_mb`: 输出 h5ad 文件大小，`os.path.getsize(output_path) / (1024 * 1024)`，保留一位小数

**文件：`templates/upload.html`**

修改 `pollConvertTask()` 函数（第 224-250 行）中完成状态的显示逻辑：

- 读取 `data.result` 中的 `sparsity` 和 `file_size_mb`
- 展示格式：
  ```
  转换完成！10,000 个细胞，20,000 个基因
  稀疏度：92.3% | 文件大小：45.2 MB
  ```

## 不改动的部分

- `io_utils.py` 中的 `convert_10x_to_h5ad()` 函数不需要修改
- `routes/upload.py` 中的 `check_10x_files()` 函数已正确处理 `.tsv` 文件名匹配
- 10x 检测和转换的路由逻辑不变
