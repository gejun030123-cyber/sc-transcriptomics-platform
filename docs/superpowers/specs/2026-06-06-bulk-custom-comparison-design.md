# Bulk RNA-seq 自定义比较功能增强设计

**日期**: 2026-06-06
**状态**: 待实施
**方案**: 方案 C — 分模块增强（各模块独立提交）

---

## 背景与目标

当前 Bulk RNA-seq 管线支持单次两组比较（group1 vs group2），无法满足以下常见分析需求：
- 多组两两配对比较（如 Control vs DrugA, Control vs DrugB, DrugA vs DrugB）
- 自定义合并组比较（如把 Treated_1h + Treated_3h 合并为 "High" 组再比较）
- 上调/下调基因分别做通路富集
- 时序实验中特定两组在各时间点的差异比较
- 手动指定基因列表绘制热图或做富集

本设计在现有 4 个模块上增量扩展参数，不引入新模块。

---

## 改动范围总览

| 模块 | 新增功能 | 新增参数 |
|------|----------|----------|
| `modules/bulk_deg.py` | 多组两两比较、自定义合并组 | `comparisons`, `custom_groups` |
| `modules/bulk_enrichment.py` | 上调/下调分开富集、自定义基因列表 | `split_direction`, `custom_genes` |
| `modules/bulk_timecourse.py` | 分组配对比较 | `pairwise_groups` |
| `modules/bulk_heatmap.py` | 自定义基因列表 | `custom_genes` |
| `routes/analysis.py` | 更新以上 4 个模块的 PARAM_SCHEMAS | — |
| `templates/analysis_select.html` | textarea 类型参数渲染支持 | — |

---

## 第 1 部分：bulk_deg 增强

### 新增参数

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `comparisons` | 多组比较（可选） | text | `` | 多个比较用分号分隔，格式：`A-vs-B;C-vs-D`。填写后 group1/group2 参数被忽略。示例：`DrugA-vs-Control;DrugB-vs-Control;DrugA-vs-DrugB` |
| `custom_groups` | 自定义合并组（可选） | textarea | `` | 格式：每行一个定义，`新组名=原组1+原组2`。示例：`High=Treated_1h+Treated_3h`。定义后可在 comparisons 中使用新组名 |

### 逻辑流程

```
输入 adata + 参数
  ↓
1. 解析 custom_groups（如有）
   - 将 "High=Treated_1h+Treated_3h" 解析为 {'High': ['Treated_1h', 'Treated_3h']}
   - 在 adata.obs 中创建临时列 '_custom_group'
   - 将原始 groupby 列的值映射到新组名（未映射的样本保留原值）
   - groupby 改为 '_custom_group'
  ↓
2. 解析 comparisons（如有）
   - 将 "DrugA-vs-Control;DrugB-vs-Control" 解析为 [('DrugA','Control'), ('DrugB','Control')]
   - 忽略 group1/group2 参数
   - 若 comparisons 为空，使用原有的 group1 vs group2 逻辑（向后兼容）
  ↓
3. 循环执行每组比较
   - 复用现有 DEG 计算逻辑（OmicVerse pyDEG）
   - 每组独立生成：火山图、MA图、箱线图
   - 文件名加比较名后缀：bulk_deg_volcano_DrugA_vs_Control.json
  ↓
4. 合并输出
   - bulk_deg_all_comparisons.csv：所有比较结果合并，含 comparison 列
   - 每个比较的独立 CSV：bulk_deg_DrugA_vs_Control.csv
   - summary 返回 n_comparisons、每个比较的 n_up/n_down
```

### 向后兼容

- `comparisons` 为空时，行为与当前完全一致（group1 vs group2 单次比较）
- `custom_groups` 为空时，groupby 直接使用原始 obs 列

---

## 第 2 部分：bulk_enrichment 增强

### 新增参数

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `split_direction` | 分开分析上调/下调 | checkbox | false | 开启后将 DEG 结果按 regulation 列拆分为 Up 和 Down 两组，分别做 ORA 富集分析，生成独立的气泡图和条形图 |
| `custom_genes` | 自定义基因列表（可选） | textarea | `` | 手动输入基因名，逗号或换行分隔。填写后忽略 input_source，直接用此列表做 ORA。示例：`TP53,BRCA1,EGFR` |

### 逻辑流程

**split_direction 模式：**

```
从 DEG CSV 读取基因列表
  ↓
按 regulation 列拆分：
  - up_genes = deg_df[regulation == 'Up']['gene']
  - down_genes = deg_df[regulation == 'Down']['gene']
  ↓
分别跑 ORA：
  - up 结果 → enrichment_ora_up_bubble.json + enrichment_ora_up_bar.json
  - down 结果 → enrichment_ora_down_bubble.json + enrichment_ora_down_bar.json
  ↓
合并结果：enrichment_split_results.csv（含 direction 列：Up/Down）
  - summary 返回 n_up_enriched, n_down_enriched
```

**custom_genes 模式：**

```
解析 textarea 中的基因名（按逗号或换行分割，去除空白）
  ↓
验证基因名是否在 pathway 数据库中存在
  ↓
直接用该列表做 ORA（与正常 ORA 流程一致）
  - 忽略 input_source 参数
  - 输出 enrichment_ora_custom_bubble.json + CSV
  - summary 注明 "input_source: custom_gene_list"
```

---

## 第 3 部分：bulk_timecourse 增强

### 新增参数

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `pairwise_groups` | 分组配对比较（可选） | text | `` | 格式：`GroupA-vs-GroupB`。指定 group_column 中的两个组名，在每个时间点做 Welch t-test，生成时序差异热图。需同时填写 group_column |

### 逻辑流程

```
前提：group_column 已填写，pairwise_groups = "Treated-vs-Control"
  ↓
解析为 (group_a='Treated', group_b='Control')
  ↓
对每个时间点 t：
  samples_a = obs[group_column==group_a 且 time==t]
  samples_b = obs[group_column==group_b 且 time==t]
  对每个基因做 Welch t-test → pvalue
  ↓
得到 基因×时间点 的 p-value 矩阵（shape: n_sig_genes × n_times）
  ↓
BH FDR 校正 → q-value 矩阵
  ↓
输出：
  - timecourse_pairwise_results.csv（gene, time, tstat, pvalue, qvalue, log2FC）
  - timecourse_pairwise_heatmap.json（X=时间点, Y=基因, 颜色=-log10(qvalue)）
  - timecourse_pairwise_log2fc_heatmap.json（X=时间点, Y=基因, 颜色=log2FC）
  - summary 返回 n_timepoints、n_pairwise_sig_genes
```

---

## 第 4 部分：bulk_heatmap 增强

### 新增参数

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `custom_genes` | 自定义基因列表（可选） | textarea | `` | 手动输入基因名，逗号或换行分隔。填写后忽略 top_var/deg 模式，直接用此列表绘制热图 |

### 逻辑流程

```
解析 custom_genes → 基因名列表
  ↓
在 adata.var_names 中查找匹配的基因
  ↓
用匹配到的基因绘制热图（与 top_var 流程一致）
  - 忽略 top_n 参数
  - 未找到的基因名在 summary 中报告
```

---

## 第 5 部分：routes/analysis.py 参数更新

在 `PARAM_SCHEMAS` 中为以上 4 个模块添加新参数定义（见各部分参数表）。

新增 `textarea` 类型支持：在模板渲染时使用 `<textarea>` 替代 `<input>`。

---

## 第 6 部分：templates/analysis_select.html 增强

在参数渲染循环中，对 `type == 'textarea'` 的参数使用：

```html
{% elif param.type == 'textarea' %}
<textarea class="form-control" id="param-{{ param.key }}" name="{{ param.key }}"
          rows="3" placeholder="{{ param.help }}">{{ param.default }}</textarea>
```

---

## 改动文件汇总

| 文件 | 改动类型 | 说明 |
|------|----------|------|
| `modules/bulk_deg.py` | 修改 | 多组比较 + 自定义合并组逻辑 |
| `modules/bulk_enrichment.py` | 修改 | split_direction + custom_genes |
| `modules/bulk_timecourse.py` | 修改 | pairwise_groups 配对比较 |
| `modules/bulk_heatmap.py` | 修改 | custom_genes 自定义基因列表 |
| `routes/analysis.py` | 修改 | 4 个模块的 PARAM_SCHEMAS 新增参数 |
| `templates/analysis_select.html` | 修改 | textarea 参数类型渲染 |

---

## 不做的事情（YAGNI）

- 不引入新的分析模块
- 不添加 volcano 图的交互式基因筛选
- 不支持 GMT 文件上传（后续迭代）
- 不添加多元方差分析（MANOVA）
- 不重构 OmicVerse pyDEG 的调用方式

---

## 实施顺序

| 顺序 | 模块 | 预计时间 | 依赖 |
|------|------|----------|------|
| 1 | bulk_heatmap（custom_genes） | 15min | 无 |
| 2 | bulk_enrichment（split_direction + custom_genes） | 30min | 无 |
| 3 | bulk_deg（comparisons + custom_groups） | 45min | 无 |
| 4 | bulk_timecourse（pairwise_groups） | 30min | 无 |
| 5 | routes/analysis.py + templates 更新 | 20min | 1-4 |

Task 1-4 可并行实现。
