# Bulk 热图增强设计文档 — 第二阶段

## Context

第一阶段已完成：5 个工具函数、28 个参数、基因选择/数据变换/聚类/视觉样式重构。第二阶段在此基础上补充 8 个高优先级功能。

## 改动清单

### 1. DEG 来源选择器（P0 bug fix）

当前 `gene_import_source=deg` 时取 `deg_files[0]`（第一个文件），多比较场景下不正确。

**新增参数：**
- `deg_comparison` — text，指定 DEG 文件 key（如 `results_0`），留空取最新

**实现：** 在 deg 分支中，若 `deg_comparison` 非空则匹配对应文件，否则取最后一个（最新）。

### 2. 上调/下调分开排列（P0）

**新增参数：**
- `up_down_separate` — checkbox，default False

**实现：** 在聚类排序后，将基因分为上调/下调两组，上调在上、下调在下，中间插入一个 NaN 行作为空隙。

### 3. 基因维度注释条（P1）

**新增参数：**
- `gene_annotation_columns` — text，逗号分隔的 `adata.var` 列名（如 `gene_type,pathway`）

**实现：** 类似样本注释条，在热图顶部或左侧生成基因注释条。从 `adata.var` 读取指定列，为每列生成独立的 Plotly figure。

### 4. 样本相关性热图参数化（P1）

**新增参数：**
- `corr_method` — select: pearson/spearman，默认 pearson
- `corr_colorscale` — select: Blues/RdBu_r/viridis，默认 Blues
- `corr_show_values` — checkbox，默认 False
- `corr_annotation` — checkbox，默认 True（是否添加分组注释条）

**实现：** 替换当前硬编码的 Pearson + Blues 为参数驱动。

### 5. 表达式筛选预设模板（P1）

**实现：** 在 `analysis_select.html` 中为 `filter_expression` 的模板按钮栏增加更多预设：
- `ALL:up` / `ALL:down`
- `ANY:up` / `ANY:down`
- `ONLY[comp]:up`（弹出比较选择）

### 6. clip_range 校验增强（P1）

**实现：** 在 `bulk_heatmap.py` 的 clip_range 解析处增加 try/except 和长度校验，格式错误时抛出清晰的 ValueError。

### 7. log_transform 独立开关（P1）

**新增参数：**
- `log_transform` — select: auto/yes/no，默认 auto

**实现：** `auto` 时根据 `adata.uns` 是否有 `normalization` 标记决定；`yes` 强制做 `np.log2(data + pseudocount)`；`no` 跳过。

### 8. 自定义注释条配色（P2）

**新增参数：**
- `annotation_palette` — text，格式如 `ctrl=#1565c0,hmc3=#e53935,rapa=#4caf50`

**实现：** 解析 `key=color` 对，替换 DEFAULT_PALETTE 中对应分组的颜色。

---

## 改动文件清单

| 文件 | 改动 |
|------|------|
| `routes/analysis.py` | 新增 7 个参数到 bulk_heatmap PARAM_SCHEMAS |
| `modules/bulk_heatmap.py` | 实现 8 个功能 |
| `templates/analysis_select.html` | 扩展筛选表达式模板按钮 |
| `tests/test_heatmap_helpers.py` | 新增测试 |

## 验证方式

1. `python -m pytest tests/ -v` — 全部通过
2. 启动 app，进入 bulk_heatmap 页面验证各新参数生效
3. DEG 来源选择器：多比较场景下能正确选择文件
4. 上调/下调分开：热图中上调基因在上、下调在下，中间有间隙
