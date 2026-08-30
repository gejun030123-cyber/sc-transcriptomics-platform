# AI 主题驱动分析能力计划（Theme-Driven Analysis）

**日期：** 2026-08-30
**状态：** 已完成（P0）
**定位：** 让 AI 工作台理解用户的生物学意图（如"脂代谢和炎症"），映射到本地通路基因集，执行或复用富集分析，并交付主题聚焦的图表与解读。
**适用范围：** 单细胞 `sc_cell_go`、Bulk `bulk_enrichment` 及后续同类"结果聚焦"场景。
**上位原则：** AGENTS.md——增量复用、最简可靠实现、按"当前必需 / 可选增强 / 明确延期"分期。

---

## 0. 背景与目标场景

用户完成差异表达后，希望在富集结果中聚焦特定生物学主题：

```text
用户：跑完差异表达了，我要关于脂代谢和炎症的富集结果。
平台：先展示命中的通路 term 供确认（如 GO BP 的 lipid metabolic process、
      inflammatory response 等），确认后运行富集并交付：
      全量统计表（不改变）+ 主题子表 + 主题 dotplot/barplot + 摘要解读。
```

### 现状盘点（2026-08-30 审查结论）

已有能力（大量复用，不重造）：

- 本地基因集为可读 GMT/TXT 快照（`GO BP/CC/MF`、`KEGG_2021_Human`、`Reactome_2022`、`WikiPathway_2021_Human`），term 名为人类可读英文（GO BP 含 87 个 lipid、24 个 inflammation 相关 term，关键词匹配可行）；
- `sc_cell_go` 已支持本地多基因集、ORA/GSEA、方向拆分、pseudobulk/细胞级分级、完整审计；
- figure engine 已支持指定通路绘图（`FigureSpec.pathway_terms` + `pathway_selection`），figure studio 已在用。

缺口（本计划要补）：

1. AI 读不到富集/DEG 表内容（`get_task_results` 只返回摘要元数据）；
2. 没有"自然语言主题 → 通路 term"的确定性映射层；
3. 任务层没有"聚焦出图"参数，`pathway_terms` 未接入 AI 与任务流。

---

## 1. 架构决策：两种执行策略

| 策略 | 做法 | 统计性质 | 结论 |
| --- | --- | --- | --- |
| A（默认） | 富集照常全库运行，FDR 在全库计算；主题仅用于展示层筛选——从全表筛出主题 term 出子表与专题图 | 诚实，不改变结论 | **采用** |
| B（延期） | 将主题 term 子集做成 GMT 再跑富集 | FDR 因检验数减少而改变，有 cherry-picking 风险 | 明确延期，若未来引入必须强制审计标注 |

## 2. 同类场景审查（共性模式）

模式：**自然语言意图 → 读取内部结果 → 语义筛选 → 聚焦交付**。

| 环节 | 场景示例 | 现状 | 处理 |
| --- | --- | --- | --- |
| 分群/注释 | 找最接近某细胞类型的 cluster | 已覆盖（signature scoring/sweep/分支/LLM 注释） | 不动 |
| 降维 | PC 该选几个 | `recommend_analysis_config` 部分覆盖 | 不动 |
| **富集** | 本计划场景 | 三项缺口 | **P0 实现** |
| QC/DEG/热图/比例/通讯/轨迹 | 各类"只看…/聚焦…" | 各自缺口 | P1 起接入同一读取框架 |

因此 P0 交付两个**通用**只读工具（`read_task_table`、`search_pathway_terms`），后续场景逐个复用，不做富集专用的一次性实现。

## 3. P0 设计（当前必需）

### 3.1 新建 `modules/theme_lexicon.py`（主题→通路映射）

- 内置中英文主题词典：首批 `lipid_metabolism`（脂代谢）、`inflammation`（炎症）、`immune_response`（免疫）、`hypoxia`（缺氧）、`emt`（上皮间质化）、`apoptosis`（凋亡）、`cell_cycle`（细胞周期）、`oxidative_stress`（氧化应激）；
- 每个 theme = 中文标签 + 中英文关键词组；中文词用于把用户输入归到主题，英文词用于匹配 GMT term 名；
- API：
  - `list_themes()` → `[{key, label, keywords}]`
  - `resolve_theme(text)` → 主题 key（无命中返回 `None`）
  - `load_library_terms(library)` → `{term: [genes]}`（解析平台本地基因集快照，路径解析与 `sc_cell_go` 一致）
  - `search_terms(query, libraries=None, limit=60)` → `[{term, library, n_genes, example_genes}]`，纯本地确定性匹配；
- 不接触任何表达数据；不上传任何内容到外部服务。

### 3.2 新增 AI 工具（`modules/ai_tools.py`）

| 工具 | 级别 | 行为 |
| --- | --- | --- |
| `search_pathway_terms` | 只读自动执行 | 输入主题词/关键词/库名 → 返回匹配 term 列表（term、库、基因数、示例基因）；支持 `list_themes=true` 列出全部主题 |
| `read_task_table` | 只读自动执行 | 读取任务登记的结果 CSV/JSON：列裁剪、行数上限（默认 50）、可选 term 包含/FDR 阈值筛选；路径必须通过 `Config._validate_path` 项目校验 |

### 3.3 `sc_cell_go` 增加 `focus_terms` 参数

- 参数：`focus_terms`（term 名列表，可由 `search_pathway_terms` 产出后由用户确认）；
- 行为：全量统计与现有输出**完全不变**；额外输出 `focus/` 子目录：
  - `focus_enrichment_<...>.csv`：主题筛选子表（含方向、FDR、重叠基因、命中主题词）；
  - 主题 dotplot/barplot（内部走 figure engine `pathway_selection='selected'`），登记 ResultFile；
- 审计：`focus_audit` 记录主题词 → 匹配规则 → 命中 term → 生成时间；focus 图表明确标注"从全量结果筛选展示，FDR 为全库校正结果"；
- 无命中 term 时如实返回审计状态 `no_match`，不伪造结果。

### 3.4 chat 编排（system prompt 增加富集意图指引）

标准对话流：

1. 用户表达主题意图 → AI 调 `search_pathway_terms`（只读）；
2. AI 把命中 term 列表给用户确认（或确认主题与库）；
3. 确认后提交 `run_analysis(sc_cell_go, focus_terms=...)`（走现有 `proposed_tools` 确认流）；
4. 完成后 `read_task_table` 读取主题子表摘要 + 附图返回。

前端复用现有确认流，无新前端组件。

### 3.5 科研诚信与安全边界

- FDR 永远在全库计算；主题图仅展示子集，图注/审计写明"从全量结果筛选展示"；
- 主题无命中 → 如实返回空并给出相近 term；
- LLM 只处理主题词与 term 名，不接触表达矩阵、细胞条码、样本元数据；
- 富集默认 local 基因集快照，不把基因列表发外部服务；
- focus 结果不替代全量结论，交付中两者并列。

## 4. 实施分期

- **P0 当前必需（本次执行）**：theme_lexicon + 2 个只读 AI 工具 + `sc_cell_go focus_terms` + chat 编排 + 测试与回归。
- **P1 可选增强**：`bulk_enrichment` 支持 focus；主题词典扩充 + LLM 归一化；QC/DEG/比例/通讯/轨迹接入 `read_task_table` 框架；figure studio 打通 focus 版本管理。
- **P2 明确延期**：WES 变异→通路主题；跨任务自动汇总报告；文献知识库注释；基因集版本目录管理（当前固定快照已够用）；策略 B（主题子集重跑）。

## 5. 测试与验收

- 单元：主题词典中英文解析、term 检索（含大小写/GO ID）、`focus_terms` 参数校验、审计字段、`read_task_table` 路径校验与行数限制；
- 集成：小数据集 pseudobulk DEG → `sc_cell_go focus`：全量表逐字节不变 + focus 子集 + PNG/SVG 登记；
- Agent 语义：模拟"跑完 DEG，要脂代谢和炎症的富集"完整工具调用序列（search → confirm → run → read）；
- 回归：不带 `focus_terms` 的现有行为与输出完全一致（AGENTS.md 第 7 条）。

---

## 6. 执行状态与续接点（2026-08-30，持续更新）

> 本节是中断续接的唯一权威记录。每完成一步就更新；续接时从这里开始，不要重做已完成项。

### 6.1 已完成

| 项 | 文件 | 状态 |
| --- | --- | --- |
| 计划文档 | `docs/ai-theme-driven-analysis-plan.md` | ✅ 本文档 |
| 主题词典与检索 | `modules/theme_lexicon.py`（新建） | ✅ 8 个主题、`list_themes/resolve_themes/search_terms`；复用 `sc_cell_go._local_gene_set_path/_read_local_gene_sets`，读 `data/go_gene_sets/` 平台快照 |
| `focus_terms` 聚焦输出 | `modules/sc_cell_go.py` | ✅ 新增 `_parse_focus_terms`/`_focus_term_mask`/`_emit_focus_outputs`；run() 解析参数、每个比较×库调用一次；输出主题 CSV + 聚焦 dotplot/barplot（`pathway_selection='selected'`，`_render_enrichment_variants` 复用）；audit 增加 `focus` 节（no_match 如实记录）；summary 增加 `focus_terms/focus_note` |
| schema 参数 | `modules/schemas.py` | ✅ `sc_cell_go` 增加 `focus_terms`（text，默认空） |
| AI 工具 | `modules/ai_tools.py` | ✅ `execute_tool` 分发 + `_search_pathway_terms`（list_themes/query/libraries/limit）+ `_read_task_table`（file_id/columns/contains/fdr_max/limit；项目归属 + `Config._validate_path` 校验；200MB 上限） |
| 工具 schema + 只读注册 | `modules/ai_adapter.py` | ✅ TOOLS_ANTHROPIC 增加两个工具定义（TOOLS_OPENAI 自动派生）；`AUTO_EXEC_TOOLS` 加入 `search_pathway_terms`、`read_task_table` |
| system prompt 编排 | `modules/ai_adapter.py` | ✅ 能力清单第 9 条 + "主题驱动富集流程"五步指引 |
| 审核修正 | `modules/theme_lexicon.py`、`modules/ai_tools.py`、`modules/sc_cell_go.py` | ✅ 对齐 `emt` 主题 key 与 `load_library_terms()` API；`read_task_table` 在列裁剪前应用 FDR 筛选，并拒绝目录符号链接；字符串化列表按安全字面量解析，保留 term 内逗号；focus audit 记录匹配规则和生成时间 |
| 测试与应用冒烟 | `tests/test_theme_focus_tools.py`、`tests/test_schemas.py` | ✅ 新增 24 项 P0 定向测试；`python -m pytest tests/test_theme_focus_tools.py tests/test_agent_tools.py tests/test_sc_de_contract.py tests/test_schemas.py -q`：83 passed；隔离可写运行目录下 `create_app()` 通过 |

### 6.2 已解决的续接问题

1. **`focus_terms` 列表被字符串化**：已在 `_parse_focus_terms` 中用 `ast.literal_eval` 安全还原 Python 列表字面量，回退到文本分隔解析；因此引号、方括号和 term 内逗号均不会污染 term 名。
2. **`ResultFile.get_by_id` 存在性**：已在 `models.py` 确认并由 `read_task_table` 使用；不存在时返回“结果文件不存在”。
3. **只读表格路径与筛选语义**：审核发现 `Config._validate_path` 仅拒绝末级符号链接，已额外拒绝项目内目录符号链接；FDR 筛选已在列裁剪前执行，且缺少可识别 FDR 列或阈值不在 `[0, 1]` 时明确报错。

### 6.3 P0 完成清单

1. ✅ `theme_lexicon`：8 个主题、确定性本地检索、单库异常隔离、库名与符号链接约束测试完成。
2. ✅ AI 只读工具：`search_pathway_terms` 和 `read_task_table` 已注册、自动执行、具备项目归属/路径/文件大小/筛选上限保护。
3. ✅ `sc_cell_go focus_terms`：全量统计不变，输出主题子表、聚焦图、匹配审计；无命中和无显著命中均如实记录。
4. ✅ 回归：不带 `focus_terms` 时不写 focus summary；列表、字符串化列表、显示名匹配和非子串匹配均有测试。
5. ✅ 提交准备：已通过定向回归、编译和应用工厂冒烟；下一步仅为创建提交并推送当前分支。是否快进 `master` 仍由用户决定。

### 6.4 P1 备忘（本次不做）

- `bulk_enrichment` 同款 focus；主题词典扩充；QC/DEG/比例/通讯/轨迹接入 `read_task_table`。

### 6.5 关键设计决定（已定，不要重议）

- **策略 A**：全量统计不变，focus 仅是展示层筛选；FDR 永远全库校正。策略 B（主题子集重跑）明确延期。
- **不做独立"重新出图"写工具**：聚焦图在 `sc_cell_go` 任务内产出，避免新建 ResultFile 写路径。
- 图形侧零开发：`figure_engine` 的 `pathway_terms/pathway_selection` 已支持，dotplot/barplot/GSEA 模板全部走 `_select_terms`/`requested_term_rows`。
- `search_pathway_terms` 是只读自动执行（在 `AUTO_EXEC_TOOLS`）；真正运行的 `run_analysis` 仍走既有确认流。
- term 匹配语义：全名精确（忽略大小写/首尾空白）或"去掉尾部 `(GO:xxxxx)` ID 后"的显示名精确匹配——不是子串匹配，避免误聚焦。
