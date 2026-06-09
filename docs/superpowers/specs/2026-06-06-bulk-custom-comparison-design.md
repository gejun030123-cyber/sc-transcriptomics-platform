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

在 `PARAM_SCHEMAS` 中为 4 个模块添加新参数定义。

### bulk_deg 新增参数

在 `PARAM_SCHEMAS['bulk_deg']` 列表末尾（`plot_genes` 参数之后）追加：

```python
        {'key': 'comparisons', 'label': '多组比较（可选）', 'type': 'text', 'default': '',
         'help': '多个比较用分号分隔，格式：A-vs-B;C-vs-D。填写后实验组/对照组参数被忽略。示例：DrugA-vs-Control;DrugB-vs-Control;DrugA-vs-DrugB'},
        {'key': 'custom_groups', 'label': '自定义合并组（可选）', 'type': 'textarea', 'default': '',
         'help': '每行一个定义，格式：新组名=原组1+原组2。示例：High=Treated_1h+Treated_3h。定义后可在多组比较中使用新组名。'},
```

### bulk_enrichment 新增参数

已在 Task 2 中完成（`split_direction` checkbox + `custom_genes` textarea）。

### bulk_heatmap 新增参数

已在 Task 1 中完成（`custom_genes` textarea）。

### bulk_timecourse 新增参数

在 `PARAM_SCHEMAS['bulk_timecourse']` 列表末尾（`fdr_threshold` 参数之后）追加：

```python
        {'key': 'pairwise_groups', 'label': '分组配对比较（可选）', 'type': 'text', 'default': '',
         'help': '格式：GroupA-vs-GroupB。在每个时间点对两组做 Welch t-test，生成时序差异热图。需同时填写分组列名。'},
```

### textarea 类型处理

前端需要将 `type == 'textarea'` 的参数渲染为 `<textarea>` 元素（见第 6 部分）。`PARAM_SCHEMAS` 中无需额外配置，类型字段 `'textarea'` 即可驱动模板渲染。

### 参数类型汇总

| 参数 | 模块 | type | 说明 |
|------|------|------|------|
| `comparisons` | bulk_deg | `text` | 单行文本，分号分隔 |
| `custom_groups` | bulk_deg | `textarea` | 多行文本，每行一个定义 |
| `custom_genes` | bulk_enrichment | `textarea` | 多行文本，逗号或换行分隔 |
| `custom_genes` | bulk_heatmap | `textarea` | 同上 |
| `pairwise_groups` | bulk_timecourse | `text` | 单行文本，格式 `A-vs-B` |

---

## 第 6 部分：templates/analysis_select.html 增强

### textarea 渲染

在参数渲染循环中，`checkbox` 块之后、`number` 块之前，插入 textarea 分支：

```html
{% elif param.type == 'textarea' %}
<textarea name="{{ param.key }}" class="form-control" id="param-{{ param.key }}"
          rows="3" placeholder="{{ param.help }}">{{ param.default }}</textarea>
```

### 渲染优先级

模板中的类型判断顺序（从上到下）：

1. `dynamic_select` — 依赖其他参数的动态下拉
2. `select` — 静态下拉
3. `checkbox` — 复选框
4. `textarea` — **新增**，多行文本区域
5. `number` — 数字输入
6. 默认 `<input type="text">` — 单行文本

### 表单提交处理

textarea 的值以字符串形式提交，后端模块自行解析（按换行/逗号分割）。前端无需额外 JavaScript 处理。

### 样式

使用 Bootstrap 5 的 `form-control` 类，`rows="3"` 提供合理的默认高度。`placeholder` 直接显示参数的 `help` 文本，提示用户输入格式。

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

## 错误处理策略

| 场景 | 处理方式 |
|------|----------|
| `comparisons` 中的组名在 obs 中不存在 | 跳过该比较，在 summary 中报告 `skipped_comparisons` |
| `custom_groups` 中的原组名不存在 | 保留原值，不创建映射（静默忽略不存在的成员） |
| `custom_genes` 中所有基因均不在 var_names 中 | 抛出 `ValueError`，提示检查基因名 |
| `pairwise_groups` 格式不正确（缺少 `-vs-`） | 跳过 pairwise 分析，在 summary 中提示格式错误 |
| `split_direction` 开启但 regulation 列缺失 | 回退到标准 ORA（不拆分方向） |
| 某个时间点样本数 < 2 | 跳过该时间点的 pairwise 比较 |

---

## 实施状态

| 任务 | 模块 | 状态 | 提交 |
|------|------|------|------|
| 1 | bulk_heatmap（custom_genes） | ✅ 已完成 | `881c818` |
| 2 | bulk_enrichment（split_direction + custom_genes） | ✅ 已完成 | `f23c95d` |
| 3 | bulk_deg（comparisons + custom_groups） | ✅ 已完成 | — |
| 4 | bulk_timecourse（pairwise_groups） | ✅ 已完成 | — |
| 5 | routes/analysis.py + templates 更新 | ✅ 已完成 | — |

> **说明**: Task 5 中，textarea 模板渲染已完成；bulk_deg 和 bulk_timecourse 的 PARAM_SCHEMAS 已在本次实现中完成（与 sc 模块一起更新）。

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

---

## 测试与验收

### 验收标准

**bulk_deg 多组比较：**
- 输入 `comparisons=DrugA-vs-Control;DrugB-vs-Control` 生成 2 组独立的火山图和 CSV
- 合并 CSV `bulk_deg_all_comparisons.csv` 包含 `comparison` 列
- `comparisons` 为空时，group1/group2 单次比较行为不变

**bulk_deg 自定义合并组：**
- 输入 `custom_groups=High=Treated_1h+Treated_3h` 后 obs._custom_group 列正确映射
- comparisons 中可引用新组名 `High`
- 未在 custom_groups 中列出的样本保留原始分组值

**bulk_enrichment split_direction：**
- 开启后生成 `enrichment_ora_up_bubble.json` 和 `enrichment_ora_down_bubble.json`
- 合并 CSV 包含 `direction` 列（Up/Down）
- regulation 列缺失时回退到标准 ORA

**bulk_timecourse pairwise_groups：**
- `pairwise_groups=Treated-vs-Control` 在每个时间点执行 Welch t-test
- 生成 `timecourse_pairwise_results.csv` 和热图 JSON
- BH FDR 校正应用于每个时间点独立

**textarea 渲染：**
- 访问 bulk_deg 分析页面，`comparisons` 和 `custom_groups` 参数正确渲染为 `<textarea>` 和 `<input type="text">`
- 访问 bulk_enrichment 分析页面，`custom_genes` 参数正确渲染为 `<textarea>`

### 验证命令

```bash
# 验证所有模块可导入
python -c "from modules.bulk_deg import BulkDEGAnalysis; from modules.bulk_timecourse import BulkTimecourseAnalysis; print('OK')"

# 验证 PARAM_SCHEMAS 完整性
python -c "
from routes.analysis import PARAM_SCHEMAS
assert 'comparisons' in [p['key'] for p in PARAM_SCHEMAS['bulk_deg']]
assert 'custom_groups' in [p['key'] for p in PARAM_SCHEMAS['bulk_deg']]
assert 'pairwise_groups' in [p['key'] for p in PARAM_SCHEMAS['bulk_timecourse']]
print('All params OK')
"

# 验证辅助函数
python -c "
from modules.bulk_deg import _parse_comparisons, _parse_custom_groups
assert _parse_comparisons('A-vs-B;C-vs-D') == [('A','B'), ('C','D')]
assert _parse_custom_groups('High=Treated_1h+Treated_3h\nLow=Ctrl') == {'High': ['Treated_1h', 'Treated_3h'], 'Low': ['Ctrl']}
print('Helper functions OK')
"
```
