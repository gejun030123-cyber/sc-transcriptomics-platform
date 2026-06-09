# Bulk DEG 差异分析增强 + 多组整合分析 设计

**日期**: 2026-06-08
**状态**: 待实施
**方案**: 增强现有 bulk_deg.py + 新建 bulk_deg_integration.py

---

## 背景与目标

当前 [modules/bulk_deg.py](modules/bulk_deg.py)（424 行）通过 OmicVerse/pyDEG 调用 3 种方法（t-test/mann-whitney/DESeq2），支持多组比较（`comparisons` 参数）和自定义合并组。

用户数据为 4 组处理（ctrl, hmc3, moclel, rapa）× 3 重复。需要：
- 增加 edgeR / limma-voom 统计方法（通过 inmoose 纯 Python 实现）
- 增加 LRT 多组检验
- 自动生成全部配对比较
- 暴露 pyDEG 已有的高级参数
- 新建多组差异整合模块（Upset 图、一致性评分、logFC 矩阵分析）

---

## 第一部分：bulk_deg 核心增强

### 新增参数总览

| 参数 | type | default | 说明 |
|------|------|---------|------|
| `method` | select | `deseq2` | 新增 `edger`/`limma` 选项 |
| `test_type` | select | `pairwise` | `pairwise`/`lrt`（LRT 仅 edger） |
| `auto_comparisons` | select | `manual` | `manual`/`all_pairwise`/`vs_reference` |
| `reference_group` | dynamic_select | `` | 参考组（vs_reference 模式自动使用） |
| `cooks_filter` | checkbox | true | Cook's 距离过滤 |
| `independent_filter` | checkbox | true | 独立过滤（低表达基因） |
| `padj_method` | select | `fdr_bh` | 校正方法：`fdr_bh`/`bonferroni`/`holm`/`fdr_by` |
| `regulation_filter` | select | `both` | `both`/`up`/`down` |

### 第 1 部分：新增统计方法

| 方法 | key | pyDEG 内部名 | 算法 | 适用场景 |
|------|-----|-------------|------|---------|
| t-test (现有) | `t-test` | `ttest` | 参数检验 | 快速预览 |
| Mann-Whitney (现有) | `mann-whitney` | `wilcox` | 非参数检验 | 不假设正态 |
| DESeq2 (现有) | `deseq2` | `DEseq2` | 负二项分布 | RNA-seq 金标准 |
| **edgeR (新增)** | `edger` | `edgepy` | 负二项+经验贝叶斯 | RNA-seq，组成偏差处理更优 |
| **limma-voom (新增)** | `limma` | `limma` | 线性模型+voom权重 | 大样本、多因素设计 |

**实现**：扩展 `method_map` + 错误处理（inmoose 未安装时给出明确提示）。

```python
method_map = {
    't-test': 'ttest', 'mann-whitney': 'wilcox',
    'deseq2': 'DEseq2', 'edger': 'edgepy', 'limma': 'limma'
}
```

### 第 2 部分：LRT 多组检验

`test_type='lrt'` 模式：使用 edgeR 的 `glmLRT`（似然比检验），一次性回答"所有组间是否有显著差异"。

- 仅 `edger` 方法可用，其他方法选 lrt 时提示并回退到 pairwise
- 输出：整体 LRT 结果 CSV（每个基因的 LRT 统计量 + p-value）
- 可作为两两比较的前置筛选（先做 LRT，再对显著基因做两两比较）

### 第 3 部分：自动生成配对比较

| `auto_comparisons` | 行为 |
|-------------------|------|
| `manual` | 使用现有 `comparisons` 参数（向后兼容） |
| `all_pairwise` | 检测 groupby 列中所有组，生成 C(n,2) 配对。4 组 → 6 个比较 |
| `vs_reference` | 所有组 vs `reference_group`。4 组（ctrl 为参考）→ 3 个比较 |

**参考组**：`reference_group` 参数指定参考组。确保 logFC 方向一致（正值 = 该组 > 参考组）。在 `all_pairwise` 模式下也影响比较顺序。

### 第 4 部分：暴露 pyDEG 高级参数

pyDEG 的 `deg_analysis()` 已支持以下参数，当前代码未暴露：

| 参数 | pyDEG 参数 | 说明 |
|------|-----------|------|
| `cooks_filter` | `cooks_filter` | Cook's 距离过滤，剔除异常高表达驱动的离群结果 |
| `independent_filter` | `independent_filter` | 独立过滤，自动去除低表达基因 |
| `padj_method` | `multipletests_method` | 校正方法选择 |
| `base_mean_filter` | 手动过滤 | 最低平均表达量过滤（当前 schema 已有但未接入逻辑） |

**base_mean_filter 接入**：在结果 DataFrame 上过滤 `baseMean >= base_mean_filter` 的基因。

### 第 5 部分：差异方向过滤 + 结果整合

**方向过滤**：`regulation_filter` 参数（`both`/`up`/`down`），在结果 DataFrame 中按 `regulation` 列过滤。

**多比较合并表**：执行多个比较时生成 `bulk_deg_merged_comparisons.csv`，列如：
```
gene | hmc3_vs_ctrl_logFC | hmc3_vs_ctrl_padj | moclel_vs_ctrl_logFC | ...
```

**多比较统计摘要**：
```python
'summary': {
    'per_comparison': {'hmc3_vs_ctrl': {'n_up': 150, 'n_down': 89}, ...},
    'shared_up': 45,  # >=2 个比较中同时上调
    'shared_down': 32,
}
```

**Volcano 并排展示**：多比较时生成一个 HTML，所有 Volcano 图并排排列。

### 第 6 部分：PARAM_SCHEMAS 更新

在 `PARAM_SCHEMAS['bulk_deg']` 中追加新参数，并将现有 `method` 的 options 扩展。

### 向后兼容

| 场景 | 行为 |
|------|------|
| 所有新参数留空/默认 | 行为与当前完全一致 |
| `auto_comparisons=manual` | 使用现有 comparisons 参数 |
| `method=deseq2` | 不变 |
| inmoose 未安装时选 edger/limma | 抛出明确 ImportError |
| `test_type=lrt` + 非 edger | 提示并回退到 pairwise |

---

## 第二部分：bulk_deg_integration 新模块

### 模块元数据

```python
MODULE_NAME = "bulk_deg_integration"
DISPLAY_NAME = "多组差异整合分析"
DESCRIPTION = "多组比较结果整合：Upset 图、一致性评分、logFC 矩阵分析"
INPUT_REQUIRES = ['bulk_deg']  # 依赖 bulk_deg 多比较输出
```

### 参数

| 参数 | type | default | 说明 |
|------|------|---------|------|
| `min_comparisons` | number | 2 | 基因至少在 N 个比较中显著才纳入 |
| `consistency_n` | number | 50 | Top N 一致性基因数 |
| `fc_threshold` | number | 2.0 | 差异基因判定的 FC 阈值 |
| `pval_threshold` | number | 0.05 | 差异基因判定的 padj 阈值 |

### 功能清单

#### 1. 差异基因集合整合

| 功能 | 输出文件 | 说明 |
|------|---------|------|
| Upset 图 | `deg_integration_upset.json` | 多组比较间差异基因的交集模式 |
| 交集/特有基因提取 | `deg_integration_gene_sets.csv` | 每种交集组合的基因列表 |

Upset 图实现：用 plotly 的 bar chart + matrix 展示各交集组合的基因数。

#### 2. 差异方向一致性

| 功能 | 输出文件 | 说明 |
|------|---------|------|
| 方向一致性矩阵 | `deg_integration_direction_heatmap.json` | 基因 × 比较，值为 +1(Up)/0(NS)/-1(Down) |
| 比较间 logFC 相关性 | `deg_integration_corr.json` | 比较间 Pearson 相关性矩阵 |

#### 3. logFC 矩阵分析

| 功能 | 输出文件 | 说明 |
|------|---------|------|
| logFC 矩阵热图 | `deg_integration_logfc_heatmap.json` | 基因 × 比较的 logFC 热图，层次聚类 |
| 比较间相关性热图 | `deg_integration_corr.json` | 比较间 logFC Pearson 相关性 |

#### 4. Consistency Score

公式：`Score = mean(sign(logFC_i) * -log10(padj_i))` for each gene across comparisons

| 功能 | 输出文件 | 说明 |
|------|---------|------|
| 一致性评分表 | `deg_integration_consistency.csv` | 每个基因的跨比较稳健性评分 |
| Top 一致性基因热图 | `deg_integration_top_consistent.json` | Top N 一致性基因的 logFC 热图 |

### 执行流程

```
输入：项目 results 目录中的 bulk_deg 多比较输出
  ↓
1. 扫描 project_dir/results 中的 bulk_deg DEG CSV 文件，解析每个比较的结果
  ↓
2. 构建基因 × 比较的 logFC 矩阵和 padj 矩阵
  ↓
3. 计算一致性评分 (Consistency Score)
  ↓
4. 识别差异基因集合（Up/Down/NS），生成 Upset 图
  ↓
5. 提取交集/特有基因列表
  ↓
6. 生成方向一致性矩阵热图
  ↓
7. 生成 logFC 矩阵热图（层次聚类）
  ↓
8. 生成比较间 logFC 相关性热图
  ↓
9. 生成 Top 一致性基因热图
  ↓
10. 保存输出
```

---

## 改动文件汇总

| 文件 | 操作 | 说明 |
|------|------|------|
| `modules/bulk_deg.py` | 大幅修改 | 新增方法 + LRT + 自动比较 + 高级参数 + 结果整合 |
| `modules/bulk_deg_integration.py` | **新建** | 多组差异整合模块 |
| `routes/analysis.py` | 修改 | bulk_deg + bulk_deg_integration PARAM_SCHEMAS |
| `modules/__init__.py` | 修改 | 注册新模块 |
| `templates/analysis_select.html` | 不改 | 使用现有 textarea/select 渲染 |

---

## 不做的事情（YAGNI）

- 不添加 K-means/Mfuzz 聚类 → 下一轮
- 不添加 PPI 网络 / Hub 基因 → 需要外部数据库
- 不添加 Meta 分析（Fisher/Stouffer）→ 下一轮
- 不添加批量 GSEA 整合 → 需独立 GSEA 模块
- 不添加转录因子活性推断 → 需 DoRothEA 数据库
- 不添加交互式比较矩阵 UI → 前端重构
- 不添加多因素模型设计公式 → 下一轮
- 不添加配对设计 → 下一轮
