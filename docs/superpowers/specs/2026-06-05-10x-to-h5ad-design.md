# 10x Genomics 三文件格式转换为 h5ad 设计文档

**日期**: 2026-06-05
**状态**: 待实施
**方案**: 方案 A — 在现有上传页面增加 10x 格式识别 + 独立转换按钮

---

## 背景与目标

单细胞转录组数据常以 10x Genomics 三文件格式（barcodes.tsv、genes.tsv/features.tsv、matrix.mtx）分发。当前平台仅接受 `.h5ad` 作为单细胞分析输入，用户需要手动转换。本设计在上传页面增加自动识别和转换功能，消除这一步骤。

**支持的格式版本：**
- v2: `barcodes.tsv`、`genes.tsv`、`matrix.mtx`（可含 `.gz` 后缀）
- v3: `barcodes.tsv.gz`、`features.tsv.gz`、`matrix.mtx.gz`

---

## 设计概览

在现有上传页面（`/projects/<pid>/upload`）增加 10x 文件组自动检测和一键转换功能。上传三个文件后，页面自动识别为 10x 数据组，显示转换表单，用户点击按钮触发异步转换任务，生成的 h5ad 文件自动出现在项目文件列表中。

---

## 第 1 部分：后端转换逻辑

### 新增函数：`modules/io_utils.py`

```python
def convert_10x_to_h5ad(mtx_dir, output_path, species=None, genome=None):
    """
    将 10x Genomics 三文件格式转换为 h5ad。
    自动检测 v2 (genes.tsv) 和 v3 (features.tsv) 格式。

    参数:
        mtx_dir: 包含 barcodes/genes/features/matrix 文件的目录
        output_path: h5ad 输出路径
        species: 可选，物种名（如 "human"、"mouse"）
        genome: 可选，基因组版本（如 "GRCh38"、"mm10"）
    """
    import scanpy as sc

    adata = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()

    # 保留 Ensembl ID（read_10x_mtx 在 var_names='gene_symbols' 时
    # 将原始 ID 存为 adata.var 的 gene_ids 列）
    # 若 gene_ids 列不存在，从 index 保存
    if 'gene_ids' not in adata.var.columns:
        adata.var['gene_ids'] = adata.var.index.tolist()

    # 可选元数据
    if species:
        adata.uns['species'] = species
    if genome:
        adata.uns['genome'] = genome

    # 保存原始计数
    adata.layers['counts'] = adata.X.copy()

    adata.write_h5ad(output_path)
    return adata
```

**信息保留清单：**

| 数据 | 存储位置 | 说明 |
|------|----------|------|
| 基因符号 | `adata.var_names` | 使用 gene_symbols 作为主索引 |
| Ensembl ID | `adata.var['gene_ids']` | 原始基因 ID |
| 特征类型 | `adata.var['feature_types']` | v3 格式：Gene Expression / Antibody Capture 等 |
| 细胞 barcode | `adata.obs.index` | 细胞标识 |
| 原始计数矩阵 | `adata.layers['counts']` | QC 模块需要 |
| 物种 | `adata.uns['species']` | 可选 |
| 基因组版本 | `adata.uns['genome']` | 可选 |

---

## 第 2 部分：API 端点

### 端点 1：检测 10x 文件组

```
GET /projects/<pid>/upload/check-10x
```

**位置**: `routes/upload.py`

**逻辑**: 扫描 `uploads/` 目录，匹配以下文件名模式：
- barcodes: `barcodes.tsv`, `barcodes.tsv.gz`
- genes: `genes.tsv`, `genes.tsv.gz`, `features.tsv`, `features.tsv.gz`
- matrix: `matrix.mtx`, `matrix.mtx.gz`

**返回 JSON**:
```json
{
  "has_10x": true,
  "files": {
    "barcodes": "barcodes.tsv.gz",
    "genes": "features.tsv.gz",
    "matrix": "matrix.mtx.gz"
  },
  "version": "v3"
}
```

`version` 判定规则：存在 `features.tsv` 或 `features.tsv.gz` → v3；存在 `genes.tsv` 或 `genes.tsv.gz` → v2。

### 端点 2：触发转换

```
POST /projects/<pid>/upload/convert-10x
```

**请求参数**（form-data，均为可选）：
- `species`: 物种名
- `genome`: 基因组版本

**执行流程**：
1. 验证三文件存在（调用 check-10x 逻辑）
2. 创建 `AnalysisTask`（module_name=`convert_10x`，status=`pending`）
3. 通过 `worker.submit_task()` 提交异步任务
4. 返回 `{"task_id": "xxx"}`

**转换完成后的文件位置**: `data/projects/<pid>/uploads/converted_10x.h5ad`

前端通过现有轮询接口 `GET /api/tasks/<task_id>/status` 跟踪进度，任务完成后在 `summary` 中返回 `n_cells`、`n_genes`、`output_file`。

---

## 第 3 部分：上传页面改造

### 3a. 10x 文件自动识别

**改造文件**: `templates/upload.html`

在现有上传区下方新增一个条件显示区域（默认隐藏）：

```
┌─────────────────────────────────────────────────┐
│  ✓ 10x Genomics 数据已就绪                        │
│                                                   │
│  barcodes.tsv.gz  ✓                               │
│  features.tsv.gz  ✓                               │
│  matrix.mtx.gz    ✓                               │
│                                                   │
│  版本: v3 (features.tsv)                           │
│                                                   │
│  物种: [人类 ▾]  (可选)                             │
│  基因组: [______] (可选)                            │
│                                                   │
│  [     转换为 h5ad     ]                           │
│                                                   │
│  转换进度: ████████░░ 80%                          │
│  状态: 正在读取矩阵...                             │
└─────────────────────────────────────────────────┘
```

**JS 交互流程**：

```
文件上传成功 → 调用 GET /check-10x
  ↓
has_10x = true → 显示 10x 就绪区域（列出文件 + 版本标识）
  ↓
用户填写物种/基因组（可选）→ 点击 "转换为 h5ad"
  ↓
POST /convert-10x → 获取 task_id → 每 2 秒轮询 /api/tasks/<id>/status
  ↓
进度更新到进度条 → 完成后显示 "转换成功！n_cells 细胞, n_genes 基因"
```

### 3b. 多文件上传支持

- 现有 `file-input` 已有 `multiple` 属性
- 改造 `fileInput.change` 事件：支持多文件顺序上传（逐个发送 XHR 请求）
- 上传完每个文件后自动触发 `check-10x` 检测
- 拖拽上传也支持多文件（从 `e.dataTransfer.files` 取所有文件）

### 3c. 转换完成后

- 显示成功提示："转换完成！n_cells 个细胞，n_genes 个基因"
- 显示生成的文件名
- "返回项目" 按钮跳转到项目详情页
- h5ad 文件自动出现在项目详情页的 "已上传文件" 列表中

---

## 第 4 部分：模块注册与集成

### 新建模块：`modules/convert_10x.py`

继承 `BaseAnalysis`，注册为 `convert_10x` 模块：

```python
class Convert10x(BaseAnalysis):
    name = 'convert_10x'

    def run(self, input_path):
        mtx_dir = self.params['mtx_dir']
        output_path = os.path.join(self.project_dir, 'uploads', 'converted_10x.h5ad')

        self.progress_callback(10, '正在读取 10x 数据...')
        adata = convert_10x_to_h5ad(
            mtx_dir, output_path,
            species=self.params.get('species'),
            genome=self.params.get('genome')
        )

        self.progress_callback(90, '正在保存 h5ad 文件...')

        return {
            'output_adata': output_path,
            'result_files': [],
            'summary': {
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'output_file': os.path.basename(output_path)
            }
        }
```

### `modules/__init__.py` 注册

在 `MODULE_REGISTRY` 中添加 `convert_10x`，**不加入 `PIPELINE_ORDER`**（数据导入不属于分析管线）。

### 与单细胞分析的衔接

- 转换后的 h5ad 保存在 `uploads/` 目录
- `projects.py` 的 `detail()` 已扫描该目录，无需改动
- 用户进入单细胞分析时选择该 h5ad 文件作为输入
- 现有分析模块无需任何修改

### 错误处理

| 场景 | 处理方式 |
|------|----------|
| 文件不完整 | check-10x 返回 has_10x=false，前端不显示转换按钮 |
| 转换失败 | 任务状态设为 failed，前端显示错误信息和 traceback |
| 文件已存在 | 覆盖写入（同一项目可多次转换） |
| matrix.mtx 损坏 | scanpy 抛出异常，捕获后设为 failed |

---

## 改动文件汇总

| 文件 | 改动类型 | 说明 |
|------|----------|------|
| `modules/io_utils.py` | 修改 | 新增 `convert_10x_to_h5ad()` 函数 |
| `modules/convert_10x.py` | 新建 | 转换任务模块，继承 BaseAnalysis |
| `modules/__init__.py` | 修改 | 注册 convert_10x 到 MODULE_REGISTRY |
| `routes/upload.py` | 修改 | 新增 check-10x 和 convert-10x 端点 |
| `templates/upload.html` | 修改 | 新增 10x 检测区域、转换表单、多文件上传、进度显示 |

---

## 依赖

- `scanpy.read_10x_mtx()`：已安装（scanpy 1.11.3），无需额外依赖
