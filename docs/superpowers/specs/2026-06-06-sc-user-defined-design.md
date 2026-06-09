# 单细胞转录组用户自定义内容增强设计

**日期**: 2026-06-06
**状态**: 待实施
**方案**: 增量参数扩展（各模块独立提交）

---

## 背景与目标

当前单细胞转录组分析管线的 10 个模块均为固定流程，用户无法指定感兴趣的基因、细胞类型或比较组。常见的自定义需求包括：

- 指定基因列表绘制 Dotplot（而非仅展示 Top N DEG）
- 指定基因在 UMAP 上展示表达分布
- 自定义 Marker 基因进行细胞注释（多行输入更易编辑）
- 指定基因沿拟时序轨迹的表达变化
- 指定特定组间比例比较

本设计在现有模块上增量扩展参数，不引入新模块。

---

## 改动范围总览

| 模块 | 新增功能 | 新增参数 |
|------|----------|----------|
| `modules/deg.py` | 自定义 Dotplot 基因列表 | `custom_dotplot_genes` |
| `modules/annotation.py` | Marker 基因 textarea 升级 | `custom_markers` 改为 textarea |
| `modules/trajectory.py` | 基因沿拟时序表达曲线 | `plot_genes` |
| `modules/proportion.py` | 指定组间比例比较 | `compare_groups` |
| `routes/analysis.py` | 更新 4 个模块的 PARAM_SCHEMAS | — |
| `templates/analysis_select.html` | textarea 渲染（bulk 侧已完成时复用） | — |

---

## 第 1 部分：deg.py — 自定义 Dotplot 基因

### 新增参数

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `custom_dotplot_genes` | 自定义 Dotplot 基因（可选） | textarea | `` | 手动输入基因名，逗号或换行分隔。填写后 Dotplot 使用此列表而非自动 Top N DEG。留空则使用默认行为 |

### 逻辑流程

```
现有逻辑：show_dotplot=True 时，取每组 Top 5 DEG 合并去重，最多 30 个基因
  ↓
新增分支：
  若 custom_dotplot_genes 非空：
    解析基因名列表（逗号/换行分割，去空白）
    在 adata.var_names 中验证
    用验证后的基因列表绘制 Dotplot
    未找到的基因在 summary 中报告
  否则：
    保持原有 Top N DEG 逻辑
```

### 向后兼容

`custom_dotplot_genes` 为空时，行为与当前完全一致。

---

## 第 2 部分：annotation.py — Marker 基因 textarea 升级

### 参数变更

| key | 变更 | 原 type | 新 type |
|-----|------|---------|---------|
| `custom_markers` | 升级为 textarea | text | textarea |

### 格式支持

保持原有格式 `CellType1:GENE1,GENE2;CellType2:GENE3,GENE4`，同时新增换行分隔支持：

```
# 两种格式均支持
CellType1:GENE1,GENE2;CellType2:GENE3,GENE4

# 或每行一个（更易编辑）
CellType1:GENE1,GENE2
CellType2:GENE3,GENE4
```

### 逻辑变更

在现有 `custom_markers` 解析逻辑中，增加对换行符的处理（将 `\n` 视为与 `;` 等价的分隔符）。

---

## 第 3 部分：trajectory.py — 基因沿拟时序表达曲线

### 新增参数

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `plot_genes` | 拟时序基因表达（可选） | textarea | `` | 手动输入基因名，逗号或换行分隔。生成这些基因沿拟时序的表达曲线图。最多 10 个基因 |

### 逻辑流程

```
在 pseudotime 计算完成后：
  ↓
解析 plot_genes → 基因名列表（最多 10 个）
  ↓
对每个基因：
  - 取该基因在所有细胞中的表达值
  - 按 pseudotime 排序
  - 计算 rolling mean（窗口 = n_cells / 50）
  - 生成折线图：X=pseudotime, Y=expression
  ↓
合并为一张多线图（每条线一个基因）
  - 输出：trajectory_gene_expression.json
```

---

## 第 4 部分：proportion.py — 指定组间比例比较

### 新增参数

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `compare_groups` | 指定比较组（可选） | text | `` | 格式：`GroupA-vs-GroupB`。仅比较这两个组的细胞比例差异，生成独立的堆叠柱状图和统计检验。留空则比较所有组 |

### 逻辑流程

```
若 compare_groups 非空：
  解析为 (group_a, group_b)
  在 adata.obs[batch_key] 中筛选这两组
  计算这两组的细胞比例 + 卡方检验
  生成独立的堆叠柱状图
  输出：proportion_compare_GroupA_vs_GroupB.json + CSV
否则：
  保持原有全组比较逻辑
```

---

## 第 5 部分：routes/analysis.py 参数更新

### deg 新增参数

在 `PARAM_SCHEMAS['deg']` 末尾追加：

```python
        {'key': 'custom_dotplot_genes', 'label': '自定义 Dotplot 基因（可选）', 'type': 'textarea', 'default': '',
         'help': '手动输入基因名，逗号或换行分隔。填写后 Dotplot 使用此列表而非自动 Top N DEG。'},
```

### annotation 参数变更

将 `custom_markers` 的 `type` 从 `'text'` 改为 `'textarea'`，help 文本更新：

```python
        {'key': 'custom_markers', 'label': '自定义 Marker（可选）', 'type': 'textarea', 'default': '',
         'help': '每行一个细胞类型，格式：CellType:GENE1,GENE2。示例：\nT_cell:CD3D,CD3E,CD2\nB_cell:CD19,MS4A1,CD79A\nMacrophage:CD68,CD163,MSR1'},
```

### trajectory 新增参数

在 `PARAM_SCHEMAS['trajectory']` 末尾追加：

```python
        {'key': 'plot_genes', 'label': '拟时序基因表达（可选）', 'type': 'textarea', 'default': '',
         'help': '手动输入基因名，逗号或换行分隔。生成这些基因沿拟时序的表达曲线图。最多 10 个基因。'},
```

### proportion 新增参数

在 `PARAM_SCHEMAS['proportion']` 末尾追加：

```python
        {'key': 'compare_groups', 'label': '指定比较组（可选）', 'type': 'text', 'default': '',
         'help': '格式：GroupA-vs-GroupB。仅比较这两个组的细胞比例差异。留空则比较所有组。'},
```

---

## 改动文件汇总

| 文件 | 改动类型 | 说明 |
|------|----------|------|
| `modules/deg.py` | 修改 | custom_dotplot_genes 逻辑 |
| `modules/annotation.py` | 修改 | custom_markers 支持换行分隔 |
| `modules/trajectory.py` | 修改 | plot_genes 拟时序基因表达曲线 |
| `modules/proportion.py` | 修改 | compare_groups 指定组比较 |
| `routes/analysis.py` | 修改 | 4 个模块 PARAM_SCHEMAS 更新 |

---

## 实施状态

| 任务 | 模块 | 状态 |
|------|------|------|
| 1 | deg（custom_dotplot_genes） | ✅ 已实现 |
| 2 | annotation（textarea 升级） | ✅ 已实现 |
| 3 | trajectory（plot_genes） | ✅ 已实现 |
| 4 | proportion（compare_groups） | ✅ 已实现 |
| 5 | routes/analysis.py 更新 | ✅ 已实现 |

---

## 实施顺序

| 顺序 | 模块 | 预计时间 | 依赖 |
|------|------|----------|------|
| 1 | deg（custom_dotplot_genes） | 15min | 无 |
| 2 | annotation（textarea 升级） | 10min | 无 |
| 3 | trajectory（plot_genes） | 20min | 无 |
| 4 | proportion（compare_groups） | 15min | 无 |
| 5 | routes/analysis.py 更新 | 10min | 1-4 |

Task 1-4 可并行实现。
