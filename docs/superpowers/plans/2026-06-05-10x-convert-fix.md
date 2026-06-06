# 10x 三文件转换功能修复与增强 实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复 .tsv 文件无法上传的问题，并在 10x 转换完成后展示数据稀疏度和文件大小预览。

**Architecture:** 两处独立改动：(1) 后端文件校验和前端 accept 属性加入 .tsv 扩展名；(2) convert_10x 模块计算稀疏度和文件大小，前端展示。

**Tech Stack:** Python, Flask, NumPy, Jinja2, Bootstrap 5

**Spec:** `docs/superpowers/specs/2026-06-05-10x-convert-fix-design.md`

---

### Task 1: 修复后端 .tsv 上传校验

**Files:**
- Modify: `routes/upload.py:9`

- [ ] **Step 1: 在 ALLOWED_EXT 中加入 .tsv**

打开 `routes/upload.py`，将第 9 行：
```python
ALLOWED_EXT = {'.h5ad', '.h5', '.csv', '.txt', '.mtx', '.gz', '.xlsx', '.xls'}
```
改为：
```python
ALLOWED_EXT = {'.h5ad', '.h5', '.csv', '.txt', '.mtx', '.gz', '.xlsx', '.xls', '.tsv'}
```

- [ ] **Step 2: 验证改动**

检查文件内容确认 `.tsv` 已加入集合。

- [ ] **Step 3: 提交**

```bash
git add routes/upload.py
git commit -m "fix: allow .tsv file upload for 10x Genomics data"
```

---

### Task 2: 修复前端 .tsv 上传支持

**Files:**
- Modify: `templates/upload.html:16-17`

- [ ] **Step 1: 更新 input accept 属性**

打开 `templates/upload.html`，将第 17 行的 `accept` 属性从：
```html
<input type="file" id="file-input" class="d-none" accept=".h5ad,.h5,.csv,.txt,.gz,.xlsx,.xls,.mtx,.tsv" multiple>
```
注意：此处 `.tsv` 已存在，但需确认。如果不存在则加入。

同时确认第 16 行说明文字包含 `.tsv`：
```html
<p class="text-muted">支持格式：.h5ad, .h5, .mtx/.tsv (10x Genomics), .csv, .txt, .xlsx, .xls</p>
```

- [ ] **Step 2: 验证改动**

在浏览器中打开上传页面，确认说明文字显示正确，文件选择器能选中 .tsv 文件。

- [ ] **Step 3: 提交**

```bash
git add templates/upload.html
git commit -m "fix: add .tsv to frontend file upload accept list"
```

---

### Task 3: convert_10x 模块增加稀疏度和文件大小计算

**Files:**
- Modify: `modules/convert_10x.py:14-35`

- [ ] **Step 1: 在 run() 方法中添加稀疏度和文件大小计算**

打开 `modules/convert_10x.py`，在文件顶部添加 import：
```python
import os
import numpy as np
from .base import BaseAnalysis
from .io_utils import convert_10x_to_h5ad
```

在 `run()` 方法中，将第 25-35 行（`self.progress(90, ...)` 到 `return`）替换为：
```python
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
```

- [ ] **Step 2: 验证改动**

检查文件语法正确：`python -c "import ast; ast.parse(open('modules/convert_10x.py').read()); print('OK')"`

- [ ] **Step 3: 提交**

```bash
git add modules/convert_10x.py
git commit -m "feat: add sparsity and file size to 10x conversion summary"
```

---

### Task 4: 前端展示转换结果预览

**Files:**
- Modify: `templates/upload.html:236-239` (pollConvertTask 完成状态处理)

- [ ] **Step 1: 修改转换完成后的显示逻辑**

打开 `templates/upload.html`，找到 `pollConvertTask()` 函数中 `data.status === 'completed'` 的分支（约第 233-239 行），将：
```javascript
if (data.status === 'completed') {
    clearInterval(interval);
    document.getElementById('tenx-convert-progress').classList.add('d-none');
    document.getElementById('tenx-convert-result').classList.remove('d-none');
    const summary = data.result || {};
    document.getElementById('tenx-result-msg').textContent =
        `转换完成！${summary.n_cells || '?'} 个细胞，${summary.n_genes || '?'} 个基因`;
```
替换为：
```javascript
if (data.status === 'completed') {
    clearInterval(interval);
    document.getElementById('tenx-convert-progress').classList.add('d-none');
    document.getElementById('tenx-convert-result').classList.remove('d-none');
    const summary = data.result || {};
    let msg = `转换完成！${summary.n_cells || '?'} 个细胞，${summary.n_genes || '?'} 个基因`;
    const details = [];
    if (summary.sparsity !== undefined) details.push(`稀疏度：${summary.sparsity}%`);
    if (summary.file_size_mb !== undefined) details.push(`文件大小：${summary.file_size_mb} MB`);
    if (details.length) msg += '\n' + details.join(' | ');
    document.getElementById('tenx-result-msg').textContent = msg;
```

注意：`\n` 在 `textContent` 中不会渲染为换行。如果需要换行显示，改用 `innerHTML` 并将 `\n` 替换为 `<br>`：
```javascript
    document.getElementById('tenx-result-msg').innerHTML = msg.replace('\n', '<br>');
```

- [ ] **Step 2: 验证改动**

在浏览器中完成一次 10x 转换，确认结果区域显示细胞数、基因数、稀疏度和文件大小。

- [ ] **Step 3: 提交**

```bash
git add templates/upload.html
git commit -m "feat: display sparsity and file size after 10x conversion"
```

---

### Task 5: 端到端验证

- [ ] **Step 1: 启动平台并上传 .tsv 文件**

```bash
cd /data/GJ/platform && python app.py
```

在浏览器中打开 `http://localhost:5000`，创建或选择一个项目，上传一组 10x 三文件（barcodes.tsv、genes.tsv、matrix.mtx，未压缩格式）。

- [ ] **Step 2: 验证 10x 检测和转换**

确认：
1. 上传后自动检测到 10x 文件组合
2. 点击"转换为 h5ad"按钮触发转换
3. 转换完成后显示：细胞数、基因数、稀疏度、文件大小
4. 点击"返回项目"能看到转换后的 h5ad 文件

- [ ] **Step 3: 最终提交（如有修复）**

如果验证过程中发现并修复了问题：
```bash
git add -A
git commit -m "fix: end-to-end verification fixes for 10x conversion"
```
