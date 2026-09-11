# 单细胞分析部分代码审查报告

> **修复状态**：本报告的 P0 已全部修复，部分 P1 已修复；逐条对照与未修复项的原因见
> `docs/single-cell-review-fixes-2026-09-10.md`。本文件保留审查当时的原始结论（含未修复项），
> 作为问题清单与优先级依据。

- 日期：2026-09-10
- 范围：`modules/` 下单细胞（sc）分析流水线及其共享底座
- 方式：**只读**代码审查（未修改任何文件、未执行写盘命令）+ 已安装依赖源码核对 + 相关测试实际执行
- 复核标记：`【已复核】` 表示该结论由审查者本人再次读代码/跑最小实验确认；其余为分模块深读结论（含行号证据）

## 1. 覆盖范围

| 分组 | 文件 |
| --- | --- |
| 共享底座 | `base.py`、`io_utils.py`、`schemas.py`、`__init__.py`、`figure_style.py`、`constants.py`、`pc_strategy.py`、`design_preflight.py`、`worker.py`、`routes/analysis.py`、`config.py`（路径校验部分） |
| 上游 | `convert_10x.py`、`sc_batch.py`、`qc.py`、`normalize.py`、`hvg.py`、`dimred.py` |
| 整合与聚类 | `batch_correct.py`、`clustering.py`、`subcluster.py`、`qc_reassess.py` |
| 注释 | `annotation.py`、`llm_annotation.py`、`cell_markers.py`、`evaluators/sc_cluster.py` |
| 下游分析 | `functional_state.py`、`scenic.py`、`sc_timecourse.py`、`deg.py`、`sc_cell_deg.py`、`sc_batch_export.py`、`sc_de_utils.py`、`sc_cell_go.py`、`enrichment_statistics.py`、`trajectory.py`、`proportion.py`、`neighborhood_da.py`、`cell_communication.py`、`virtual_ko.py`、`celloracle_worker.py` |
| 出图与登记 | `sc_figure_diagnostics.py`、`native_figures.py`、`visualization.py`、`reporting/result_manifest.py`、`reporting/pipeline_report.py`、`reporting/project_report.py`、`reporting/static_rendering.py` |

实际执行的相关测试：`test_qc_semantics / test_normalize_semantics / test_hvg`（50 passed）、`test_sc_de_contract / test_sc_batch / test_subcluster / test_proportion / test_neighborhood_da / test_scenic / test_functional_state / test_signature_scoring / test_grouping_semantics / test_deg_comparison_scope / test_worker_sc_cell_go`（103 passed）—— **153 项全部通过，未见既有回归**。

## 2. 总体结论

单细胞部分的**设计纪律明显高于同规模科研平台的平均水平**：重复单位（细胞 ≠ 生物学重复）在多数模块被显式建模，pseudobulk 先求和再标准化、背景基因集、路径白名单、任务绑定溯源、证据等级写入 manifest 等做得相当扎实。问题集中在三类：

1. **必然失败或必然错标的缺陷**（scVI 分支、失败任务 manifest 写 completed、pickle 子进程继承密钥）；
2. **"文档/界面声称"与"代码实际行为"不一致**（P 值口径、FDR 检验族、过滤参数被忽略、单位语义混用），这类问题最危险，因为读者会相信声明；
3. **静默降级/静默 fallback**（评分回退、背景丢失、分母重归一化、缺测填 0），它们不会报错，只会悄悄改变科学结论。

分级统计（去重后）：

| 级别 | 数量 | 主要影响面 |
| --- | --- | --- |
| P0 阻断/高危 | 12 | 功能必然失败、结果被静默裁剪、审计结论相反、密钥泄露 |
| P1 中危 | 27 | 统计口径、静默降级、资源与可复现性 |
| P2 低危 | 30+ | 结构重复、参数缺界、日志与命名 |

### 建议优先修复（Top 8）

| # | 位置 | 一句话 |
| --- | --- | --- |
| 1 | `batch_correct.py:292` | scVI 分支引用未定义变量 `layer`，该路径 100% 失败 `【已复核】` |
| 2 | `virtual_ko.py:184` + `celloracle_worker.py:109/137` | 用户上传 pickle 在继承全部环境变量的子进程中反序列化，平台口令与 AI 密钥可被读取 |
| 3 | `deg.py:224/722` | pairwise 模式把"仅参与比较的组"写成标准输出 `deg_output.h5ad`，下游静默按裁剪后数据分析 `【已复核】` |
| 4 | `reporting/result_manifest.py:63` | manifest 的 `status` 硬编码 `"completed"`，失败任务留下"已完成"审计记录 `【已复核】` |
| 5 | `qc_reassess.py:83-99` | 同一阈值同时比较"doublet 比例/平均 score/score>0.5 比例"三种量，且列名固定为 `doublet_fraction` `【已复核】` |
| 6 | `functional_state.py:929-949, 2041` | 逐细胞（非独立）Spearman 的 p 值写进 CSV，与 manifest"intentionally omit P"直接矛盾 |
| 7 | `sc_figure_diagnostics.py:757-759` | 超过 `max_labels` 时截断后**重新归一化**分母，该比例表还会进入 summary `【已复核】` |
| 8 | `sc_cluster.py:145-229` | marker 评分不校验表达尺度，raw counts 下 `0.35*均值` 直接饱和到 1.0 → 输出 `high` 置信度 |

## 3. P0：阻断级 / 高危

### P0-1 scVI 分支必然 NameError
- 位置：`modules/batch_correct.py:292`（`_run_scvi`）
- 证据：`return 'X_scVI', {'max_epochs': max_epochs, 'n_hvg_for_scvi': int(work.n_vars), 'layer': layer}`；`_run_scvi` 内只有 `setup_anndata(work, layer='counts', ...)` 关键字参数，全文无 `layer =` 赋值 `【已复核】`
- 影响：`method='scvi'`（schema 正式选项）先训完 80 epochs 再在最后一行崩溃，任务判失败且无任何产物/图件；属 100% 不可用路径
- 建议：改为 `'layer': 'counts'` 或删除该键，并把实际训练层写入 `method_info`

### P0-2 CellOracle 子进程 pickle 反序列化 + 继承全部环境变量
- 位置：`modules/virtual_ko.py:184-191`、`modules/celloracle_worker.py:104-109,137`
- 证据：`env = dict(os.environ)` 后仅补充 HOME/MPLBACKEND 等，未剔除密钥；`co.utility.load_pickled_object(path)` 对 `.pkl/.pickle/.gpickle/.oracle/.celloracle/.links` 直接反序列化；`config.py:31 _load_dotenv()` 把 `.env` 全量写入 `os.environ`（含 `PLATFORM_ACCESS_PASSWORD`、`AI_API_KEY`、`AI_API_TOKEN`），而上传白名单与 base GRN 选择器都允许这些扩展名
- 影响：能上传文件并运行 virtual_ko 的用户可借 pickle 在子进程执行任意代码，读取平台登录口令与 AI 密钥；`links_file` 分支连提示都没有，base GRN 的"风险提示"只是写进 summary（RCE 已发生）
- 建议：子进程环境白名单构造并剔除 `*_KEY/*_TOKEN/*PASSWORD/SECRET*`；禁用 pickle 载入或改专用沙箱用户；未受信资源在载入前直接阻断

### P0-3 pairwise DEG 把比较子集写成标准输出
- 位置：`modules/deg.py:218-226` 与 `:722`
- 证据：`adata = adata[adata.obs[groupby].astype(str).isin(involved)].copy()` 之后 `output_path = self.save_output(adata, 'deg')` `【已复核】`
- 影响：`intermediate/deg_output.h5ad` 只含参与比较的组，summary 无裁剪记录；`routes/analysis.py:254-266` 的下游输入选择器按 `INPUT_REQUIRES`（该文件字段齐全）收录它，trajectory/proportion/cell_communication 会静默按裁剪后数据计算
- 建议：统计用子集、落盘用全量（或在输出上显式记录 `cells_dropped_for_contrast`）

### P0-4 manifest 状态恒为 completed
- 位置：`modules/reporting/result_manifest.py:63`，配合 `worker.py:254-261`
- 证据：`"status": "completed",` 为字面量，`build_task_manifest` 无状态入参；worker 先 `register_task_outputs`（写 manifest）后才 `_result_error` → `mark_failed` `【已复核】`
- 可达路径：`cell_communication.py:59-70`（LIANA 未安装）、`batch_correct.py:971`
- 影响：交付/审计核心文件对失败任务写"已完成"且 `result_files` 为空；`project_report.py` 还会采信其中的 `review_evidence`
- 建议：`build_task_manifest` 增加 `status/error` 入参（或从 `result['error']`/`summary['error']` 推导），失败时写入 error 字段

### P0-5 QC 重评估的 doublet 指标语义混用
- 位置：`modules/qc_reassess.py:83-99`（输出列在 `:111-119`，阈值默认 0.3）
- 证据：三条分支分别得到「`predicted_doublet` 均值（比例）」「簇平均 `doublet_score`（连续打分的均值）」「`doublet_score > 0.5` 的比例」，全部写入同一列 `doublet_fraction` 并与"Doublet 比例阈值"比较 `【已复核】`；`doublet_frac == 0.0` 的精确浮点比较会让同表不同簇采用不同定义
- 影响：低质量簇标记与 `auto_remove` 删细胞建立在语义错位的指标上；报告读者会把平均 score 当双细胞比例
- 建议：拆为 `doublet_fraction`（仅比例）与 `mean_doublet_score` 两列/两阈值，列名与图注写清定义

### P0-6 逐细胞 p 值写入结果表，与 manifest 声明矛盾
- 位置：`modules/functional_state.py:925-949`（`add_row`）、`:1039-1040`、`:2041`，manifest `:2366-2372`
- 证据：`spearmanr(...)` 的 `pvalue` 被写入 `cell_level_descriptive` 层级（`n_units` = 细胞数）并整表落盘 `activity_correlations.csv`；manifest 却写 "descriptive only ... intentionally omit P"
- 影响：CSV 中出现基于上万非独立细胞的 P 值（通常≈0），正是模块注释声明要避免的伪重复；表与图口径不一致，下游极易引用
- 建议：descriptive 层级将 `pvalue` 置 NaN，或改名为 `descriptive_p` 并在 `analysis_note` 标注非推断性

### P0-7 组成比例截断后重新归一化
- 位置：`modules/sc_figure_diagnostics.py:757-759`，调用方 `annotation.py:4051-4060` → `:4796` 进 summary
- 证据：`counts = counts.loc[:, list(top_labels)]` 后 `fractions = counts.div(counts.sum(axis=1)...)`——分母变成 Top-N 之和 `【已复核】`；同文件 `condition_composition_figure:1044-1053` 用显式 `Other` 段保留了 100% 分母
- 影响：注释标签 >30（细注释常见）时，`donor_celltype_composition` 的比例被系统性放大，计数表也被静默截断
- 建议：对齐 `condition_composition_figure` 补 `Other` 段，或返回 `dropped_labels` 并在图/summary 声明"仅展示 Top-30，比例基于全部细胞"

### P0-8 细胞类型签名评分：无表达尺度校验 + 无界项加权 → 置信度失真
- 位置：`modules/evaluators/sc_cluster.py:145-229, 296-306`（`annotation.py`、`ai_tools.py`、`routes/branches.py` 四条调用链共用）
- 证据：`positive_score` 直接取 `adata[mask, g].X` 均值，无归一化/尺度校验；`total_score = 0.35*pos + 0.25*spec - 0.20*neg + 0.10*size - 0.10*batch` 后 `clip(0,1)`；置信度阈值硬编码 0.75/0.50/0.25 `【已复核】`
- 影响：在 raw counts（例如 `qc_output.h5ad`）上 `0.35 × 平均 UMI（常 2–10）` 即饱和到 1.0 → `confidence='high'`；`X` 已 scale 时 `positive_score` 可为负。评分完全依赖输入尺度，且簇大小/单样本占比这些与 marker 证据无关的量参与"置信度"
- 建议：评分前用 `resolve_expression_measurement` 校验尺度并只接受非负 log1p；`size/batch` 拆为独立质控提示；置信度只用 marker 特异性与检出率

### P0-9 评分失败时静默归零
- 位置：`modules/evaluators/sc_cluster.py:147-154, 160-175, 180-188, 205-219`
- 证据：每个基因/每个簇的表达式读取都包在 `try: ... except Exception: pass` 中 `【已复核】`
- 影响：backed 索引不兼容、列不存在、稀疏转换失败等任一原因都会让 `pos_exprs`/`spec_scores` 为空 → 分数 0.0 → 上层返回 `confidence='very_low'` 而**无任何错误**。AI 工具据此告诉用户"没有匹配的簇"，属静默科学结论
- 建议：至少记录异常计数与首条异常文本；全失败时返回 `error` 而不是 0 分

### P0-10 LLM 注释批次与 max_tokens 结构性失配
- 位置：`modules/llm_annotation.py:321`（Anthropic）、`:339`（OpenAI），配合 `annotation.py:3703-3709`
- 证据：两分支 `"max_tokens": 4096` `【已复核】`；提示词要求"每个 cluster 一项"，每项 `rationale`+`review_note` 各 240 字符，默认一批 30 簇；缺任一 cluster 直接抛错，异常被转成整任务 `RuntimeError`，无降批重试
- 影响：中文输出极易截断 → JSON/覆盖校验失败 → 一条长流水线的末尾任务失败，已成功批次全部丢弃
- 建议：按 `per_batch` 动态设 `max_tokens`（或默认批降到 8–10）；对截断/缺项实现一次降批重试；逐批结果先落盘

### P0-11 零证据簇被恢复为确定 celltype 且不标 needs_review
- 位置：`modules/annotation.py:2787-2801`、`:2835-2839`、`:3033-3042`
- 证据：`(not np.isfinite(candidate_score) or candidate_score <= 0) and mean_detection <= 0` 时，只要 `annotate_all and score_margin > 0 and agreement >= 0.6` 就保留候选标签并给 `evidence_tier='provisional'`；`needs_review` 的判定条件在默认参数（`confidence_method='none'`）下全为假
- 影响：候选 marker 一个都没检出的簇仍得到确定细胞类型与 ontology 层级，复核表 `needs_review=False`，与真实有支持的 provisional 调用无法区分
- 建议：该 `decision_reason` 强制 `needs_review=True` 且层级降为 Unresolved；或要求候选 marker 均值 > 0 才允许恢复

### P0-12 固定输出文件名导致历史任务产物被覆盖
- 位置：`modules/convert_10x.py:105`、`:250`、`:256`
- 证据：`uploads/converted_10x.h5ad`、`uploads/combined_batches_imported.h5ad`、`results/batch_import_summary.csv` 为硬编码名 `【已复核】`；同项目的 `sc_batch.py:470` 已正确使用 `_available_path()`
- 影响：二次导入覆盖前次文件，而 `worker.py:261` 已把该路径写进旧任务的 `output_adata_path` → 旧任务"已完成"但文件内容已变，溯源失效（违反 AGENTS.md 版本不得混用）
- 建议：统一复用 `_available_path()`；必须覆盖时记录被覆盖对象的哈希

## 4. P1：科学性与统计正确性

### 4.1 富集与 FDR
1. **cell_level 与 pseudobulk 的 ORA 走两套实现，FDR 族语义不同**（`sc_cell_go.py:1531-1541` vs `:1598-1602`）。`run_ora_full` 的契约是"背景中至少一个成员的通路都进入 BH 族（含 0 命中 P=1）"，gseapy 会丢弃 0 命中通路再做 BH → 分母变小、FDR 系统性偏松，而 cell_level 图标题仍写"全库 FDR"。`【已复核：两条分支确实分别调用 run_ora_full 与 gp.enrich/gp.enrichr】`
2. **`ora_min_size/ora_max_size` 在 cell_level（默认路径）被静默忽略**（`sc_cell_go.py:1210-1211` → `_run_ora` 只传 `background`），而 schema 明确承诺按背景过滤过小通路 → 1–4 基因的超小通路参与检验并计入 FDR。`【已复核】`
3. **`execution_mode='enrichr'` 丢弃 tested-gene background**（`sc_cell_go.py:214-223`），Enrichr 用自身全基因组背景，但结果表/审计仍写 tested universe 的 `n_background_genes`；代码注释自己就指出这会高估富集。`【已复核】`
4. **GSEA 导出表 `P-value` 为空串**，真实 p 只在额外列 `pval`（`sc_cell_go.py:1561-1563`），按列名读取的下游拿到空值。
5. **cell_level 不做标识符规范化/大写化**（`:1436-1443`），pseudobulk 才建 `identifier_map` → 同一 DEG 表在两个层级命中率不可比。
6. **`local_gene_set_dir`/本地库不校验物种**（`:150-166`、`:118-129`）配合 `normalise_human_gene` 的 `upper()`，小鼠 GMT 会被静默当人类库使用。

### 4.2 差异表达
7. **默认 reference 模式无最小组保护**（`deg.py:160, 265-275, 313-321`）：`min_cells` 只对 pairwise/custom 生效，任一 celltype 只剩 1 个细胞时 scanpy 抛 `ValueError`，代码只捕 `KeyError` → 整任务以英文异常崩溃。
8. **marker 契约由列名猜测**（`deg.py:588-597, 762`）：`marker_contract='celltype_marker_vs_rest'` 与 `analysis_grouping` 与真实分组无关；`groupby='condition'` 的条件比较被贴上 celltype marker 溯源标签，并被 `routes/analysis.py:130-133`、`sc_cell_go.py:1353-1361` 消费。
9. **`pct.2` 一律 `fillna(0.0)`**（`sc_cell_deg.py:331-347`）：对照组缺该基因时伪造"对照组不表达"，无告警。
10. **Welch 的 log2FC 与所检验统计量不同源**（`sc_batch_export.py:112-118`）：用 `log2(均值CPM比)+0.5` 伪计数，而检验基于 log2CPM t 统计量；与 `bulk_deg.py:292` 的"log 尺度均值差"定义也不一致 → 同一列名三处语义不同、阈值不可比。
11. **只要任一单元成功即报 `biological_replicate_model` / `completed`**（`sc_batch_export.py:1069-1073, 1130-1135`），部分 DESeq2 拟合失败仍显示完成。
12. **`sc_cell_deg` 的 `export_full_tables` 不在 schema**（`:253, 387, 523-529`）→ GUI 恒走 False 分支，公开 CSV 在 Web 上不可达（死代码）。
13. **`output_adata: input_path`**（`sc_cell_deg.py:533`、`sc_batch_export.py:630,1150`）使 worker 把整份上游 h5ad 复制进 `task_artifacts/`，GB 级重复占用，结果页"输出"实为上游副本。

### 4.3 时序与功能状态
14. **`sc_timecourse` 图上把缺测填 0**（`:613-631`）：`fillna(0)` 让"无样本的时间点"读成"该细胞类型占比 0/消失"；CSV 侧反而诚实。
15. **n=2 仍标正式样本级推断**（`sc_timecourse.py:520, 143-152`）：与 `functional_state._sample_size_evidence_tier` 的 `exploratory_n2` 口径不一致，`inference_status='sample_level_kruskal'` 易被读成"通过检验"。
16. **`direction` 不看证据状态**（`sc_timecourse.py:418-421, 481-484`）：`descriptive_only`/`test_unavailable` 的行照样得到 `increasing/decreasing`。
17. **TF 网络权重缺失 `fillna(1.0)`**（`functional_state.py:335-343`）：等于把缺失调控方向当作"激活 +1"，对有符号 ULM 活性产生系统性正向偏倚且无告警。
18. **基因集重复 QC 硬失败**（`functional_state.py:452-484`）：检出稀疏时两个基因集的 detected 集合可能相同 → 整个分析 `raise`，UI 无豁免开关。
19. **BH 校正族边界未声明**（`functional_state.py:758-800`）：`gene_expression/pathway_score/tf_activity × 所有 celltype` 混为一个校正族。

### 4.4 比例、邻域、通讯、轨迹
20. **proportion / neighborhood_da 默认 `min_samples_per_condition=2`**：双侧 Mann-Whitney 在 n1=n2=2 时最小 p=0.3333，**任何差异都不可能显著**，且不给"样本量不足"提示 → 用户误读为"没有差异"。
21. **邻域重叠未做空间 FDR**（`neighborhood_da.py:274-278`）：高度共享细胞的 KNN 邻域直接 BH，FDR 被低估（Milo 正是为此设计空间 FDR）。
22. **proportion 把"低细胞数组过滤"后的对象写成 intermediate 输出**（`:369-376, 512`），下游复用该文件会掉细胞类型。
23. **`fisher_exact` 非 2×2 时静默降级为卡方**（`proportion.py:16-22`），summary 仍标 `fisher_exact`，`statistic` 与 `chi2` 字段口径不一致。
24. **cell_communication 用未归一化 counts 作为 LIANA 输入**（`:85-97`），与 LIANA 官方 `normalize_total+log1p` 用法不一致，`lr_means/magnitude_rank` 受文库大小影响；且未记录 liana/资源版本与校验和。
25. **"No significant interactions" 表述错误**（`cell_communication.py:139-149`）：空表只是被 `expr_prop` 过滤，`rank_aggregate` 返回全部 LR 对；`top_interactions = head(top_n)` 也未显式按 `magnitude_rank` 排序。
26. **trajectory 邻居表示候选漏掉平台真实键名**（`trajectory.py:43-49`）：候选里有 `X_harmony`/`X_combat`，而平台实际写 `X_pca_harmony`/`X_pca_combat` → Harmony/ComBat 数据上静默回退到**未校正** `X_pca` 重建邻居图，无警告。
27. **根节点选取任意且复用陈旧 `uns['iroot']`**（`trajectory.py:94-109`）：层内取"第一个细胞"使原点依赖细胞顺序；未清理上游遗留 `iroot` 会把指向别的细胞顺序的索引当根。
28. **unresolved 簇的显示映射被硬编码**（`annotation.py:462-476` + `:3912-3920`）：`'Unresolved epithelial programme: CLDN2/NME1 (cluster 3; review)' → 'Inflammatory stress cells'`，任何复用该面板且 cluster 恰为 3 的数据都会把 `unresolved-review` 显示为确定类型并赋 resolved 层级。
29. **旧 `score_*` 列被当成本次证据**（`annotation.py:3506, 3847`）：换面板重跑时上一次运行的 score 列会污染 `lineage_mixture_score`/`annotation_confidence`/`annotation_top1/2`。
30. **谱系混合判定只覆盖 4 个带亚型谱系**（`annotation.py:3540-3559`），summary 与 schema 却称 broad lineage 比较。

## 5. P1：工程健壮性、资源与可复现性

31. **worker 双重加载 h5ad**（`worker.py:234` 与 `:252`）：preflight 的 `adata` 在局部作用域中一直存活到 `module.run(input_path)` 结束（`run` 内部会再读一次）→ 大数据集峰值内存约 2× 完整数据集。建议在验证后 `del adata`。`【已复核：变量确实跨越 run 调用存活】`
32. **`qc.py:1536-1538` 在重复基因名上抛 `InvalidIndexError`**：`original_var_names.get_indexer(adata.var_names)`，而 `Index.get_indexer` 要求索引唯一 `【已复核：最小实验确认会抛 InvalidIndexError】`；`load_adata → remap_var_names` 在 var_names 已是 symbol 时于 `io_utils.py:1011` 提前返回、不做唯一化，`io_utils.py:1049` 重映射后也不再 `var_names_make_unique` → 带重复基因名的 h5ad 会在**全部 QC 计算与出图完成后**才崩溃。建议 QC 入口调用 `var_names_make_unique`。
33. **QC 散点图错位配对**（`qc.py:1427-1433`，同类 `:1451-1455`）：三列各自独立过滤有限值后按 `[:count]` 切片配对，行不再对应同一细胞；若 `pct_counts_mt` 的有限值少于 `count`，`c` 与 `x` 长度不等直接 `ValueError` 让整个 QC 失败。建议一次构表后 `dropna(how='any')`。`【已复核】`
34. **全基因 PCA 分支在主对象上 scale 并持久化稠密 `layers['scaled']`**（`dimred.py:259-261`）：10 万细胞 × 2 万基因 ≈ 8 GB/文件写进每个中间 h5ad，而 schema 帮助文本恰称该选项"减少内存压力"。
35. **陈旧校正表示未清理**（`dimred.py:18-25` vs `clustering.py:91-93, 116-119`）：重跑 dimred 后旧 `X_pca_harmony` 等仍在，clustering `use_corrected=true` 会把聚类建在与当前 PCA 不一致的旧空间。
36. **PC 肘部维度截断潜在空间**（`pc_strategy.py:38-41` + `clustering.py:771-773`、`subcluster.py:165-167`）：对 scVI/SysVI/Scanorama 也套用 PCA 肘部维度，可能只保留前几维潜变量。
37. **marker 预览固定用首个分辨率**（`clustering.py:878-881` vs `:919-926`）：主分辨率为 0.8/1.0 或 auto_select 选中非首分辨率时，dotplot/heatmap/violin 与 `uns['marker_selection']` 仍按 0.6 的簇编号，同页 cluster 编号互不对应。
38. **subcluster 覆盖 `X_umap` 且写标准输出**（`subcluster.py:144/181/296`）：作为链式输入时，其后所有模块只分析目标簇细胞，且父级全局坐标丢失。
39. **进程级随机种子污染并发任务**（`batch_correct.py:200`）：`scvi.settings.seed = ...` 会改写 torch/numpy/python 全局 RNG，而 worker 是 2 线程进程池；`sysvi_seed`/`sysvi_epochs` 又不在 schema 中，UI 无法设置。
40. **`plt.close('all')` 在 worker 线程内执行**（`annotation.py:4304-4306` `【已复核】`）：会关掉另一并发任务正在构建的 figure，表现为偶发图件缺失/`savefig` 报 figure 已关闭。应改为 `plt.close(fig_dotplot)`。
41. **logreg 注释模式每簇全量 copy**（`annotation.py:2126-2143`）：最多 200 簇 × 整个 AnnData（含 layers/raw/obsm），大图谱直接 OOM；`wilcoxon` 路径同一 run 也可能全量复制 3 次（`:3388/:3679/:3990`）。
42. **`task_artifacts` 无保留策略**：`worker.py:35-77` 把每个任务的完整 `output_adata.h5ad` 复制一份，实测单项目 `results/task_artifacts` 已达 **122 GB**（58 份 2.2 GB 的 h5ad），且同 task 重跑会覆盖路径而旧 `ResultFile` 行仍指向它。
43. **失败/不可用任务仍登记上游 h5ad 为输出**（`cell_communication.py:61, 120, 143`）：返回 `'output_adata': input_path`，worker 在 `_result_error` 之前就完成快照复制 → 失败任务仍占一份完整副本，且用户/AI 可能误读为"产出了数据"。
44. **virtual_ko 超时只 `proc.kill()` 直接子进程**（`virtual_ko.py:198-201, 235, 267`）：无 `start_new_session`/进程组、不记录 PID，CellOracle 的 joblib/loky 子进程会成为孤儿继续占 CPU；默认超时 6 小时，且单细胞侧**没有取消入口**（全仓 cancel 只在 WES/workflow）。与 AGENTS.md"取消必须真正终止子进程"不符。
45. **celloracle_worker 稠密化 Δ 矩阵 + 逐基因整层写入**（`:443-446, 528-530`）：`delta_x.toarray()` 在 3 万细胞 × 2 万基因时是数 GB（列均值本可稀疏计算）；`simulated_count_<label>` 每基因增加一整层；HVG 收敛仅在 `n_vars > 4000` 且已有 `highly_variable` 列时生效。
46. **`n_propagation` 无上下限**（`schemas.py:531` + `virtual_ko.py:134`）：CellOracle 仅接受 1–5，越界会在跑完 GRN 推断后才 `ValueError`；非 base GRN TF 的扰动基因只发警告但 CellOracle 必然硬失败（`celloracle_worker.py:372-375`）。
47. **可复现性口径不一**：`functional_state` manifest 记录 `input_sha256`，而 `scenic.py`/`sc_timecourse.py` 不记录；`result_manifest` 全链路无校验和，也不记录 scanpy/omicverse/anndata 版本（87 份真实 manifest 实测 `has checksum: False`）。
48. **静默降级无记录**：`functional_state.py:602-614`（评分回退到无对照均值）、`scenic.py:84-93`（用户指定的 `auc_obsm_key` 不存在时静默换默认键）、`dimred.py:337-338`（kneed 缺失时回退但 summary 仍写 `kneedle`）、`proportion.py:16-22`（fisher→卡方）。
49. **收尾阶段无整体异常保护**（`sc_cell_go.py:1863-1880`）：统计表已算完落盘，但收尾 import/`to_csv`/审计写入抛错会让任务判失败且结果不登记。
50. **`neighborhood_da` 全 NaN 绘图**（`:407-413`）：`max(np.nanmax(...), 0.25)` 在全 NaN 时返回 nan → `TwoSlopeNorm(vmin=nan)` 抛错。
51. **`min_prop` 等关键参数无界**（`schemas.py:516-517, 531, 540-542`）：`min_prop > 1` 时 LIANA 返回空表，界面显示"无显著相互作用"而非参数错误。

## 6. P1：安全与隐私

52. **绝对路径随 `result_json` 外发**（`sc_batch_export.py:1127, 1137, 1139-1140` 等 → `ai_tools.py:493-500`）：`csv_package_dir`/`pseudobulk_qc_files`/`deg_source_file` 为绝对路径，worker 只重写被快照的 `result_files`，这些字段保持原样并被 `get_task_results` 整个回传给外部 LLM；`:595, 614-616` 还把绝对路径写进可下载的 `*_export_metadata.json`。违反 AGENTS.md 底线 #1。
53. **`manual_map_csv` 读取任意绝对路径**（`annotation.py:3925-3928` → `:708` `pd.read_csv`）：只做 `os.path.isfile`，未限制在项目目录内，任意可读 CSV 的内容会并入 `celltype`、UMAP 图例与结果表（读取-外带通道）。`【已复核】`
54. **LLM 外发消毒仅拦邮箱与个别路径前缀**（`llm_annotation.py:61-70` `【已复核】`）：`/mnt/...`、UNC、`data/projects/P001`、"donor P-2023-014" 都能通过；唯一防线是 UI 文案提醒。
55. **`dimred` 用未净化的 obs 列名拼文件名**（`dimred.py:389, 425`）：`batch_key` 是用户可控的 obs 列名，`base.save_matplotlib_figure` 只做 `splitext` 不做字符清洗 → 形如 `x/../../..` 的列名可让图片与 `*_nature_readiness.json` 写到 `plots/` 之外。（`clustering.py:198-201`、`subcluster.py:111` 已做了清洗，dimred 是唯一遗漏点。）
56. **ZIP 解压无体积/条目上限，全平台无 `MAX_CONTENT_LENGTH`**（`convert_10x.py:19-29`）：高压缩比 ZIP 可撑满数据盘。
57. **virtual_ko 路径白名单包含整个 `Config.DATA_DIR`**（`virtual_ko.py:56-70`）：`base_grn_file`/`links_file` 为自由文本，可跨项目引用他人上传的 pickle，绕过项目边界。
58. **批次名唯一性在清洗前校验**（`convert_10x.py:169` vs `:190`）：中文批次名清洗后重名 → `obs_names` 前缀与 `obs['batch']` 标签冲突，无告警。
59. **绝对服务器路径进入 obs/summary**（`sc_batch.py:431, 439, 455`）：`source_matrix_dir` 等列不在 `TECHNICAL_OBS_COLUMNS` 中，而 AI 工具会把合格 obs 列取值纳入模型输入。
60. **外部服务外传无审计**（`sc_cell_go.py:220-223`）：`enrichr` 会把显著基因符号 POST 到 maayanlab.cloud（仅基因符号），但审计只记 `execution_mode`，未记录"已外传基因列表"。

## 7. 跨模块 / 平台级一致性问题

61. **设计预检在两个入口不一致**：Web 表单只对 `{bulk_deg, sc_timecourse, batch_correct, bulk_normalize}` 执行 `preflight_blockers`（`routes/analysis.py:502`），AI 路径却额外覆盖 `{functional_state, sc_pseudobulk_deg, proportion, neighborhood_da}`（`ai_tools.py:280-289`）。即**人通过网页提交时绕过了 AI 被强制的阻断**，而 `sc_pseudobulk_deg`（正式样本级统计）正是最需要预检的模块；`design_preflight.py:561-640` 为 `sc_cell_deg` 写的检查在两条路径都不会被执行。`【已复核】`
62. **`INPUT_REQUIRES` 缺失导致输入选择器过度宽松**：`sc_cell_deg`、`sc_pseudobulk_deg`、`sc_cell_go`、`sc_batch_import`、`sc_csv_export` 未声明，`functional_state` 声明为空 `[]`；而 `build_input_options`（`routes/analysis.py:263`）用它过滤输入 → 用户可选中不含 `celltype`/样本列的上游输出，直到运行期才失败。`【已复核】`
63. **`PIPELINE_DEPS` 未覆盖 `sc_cell_go`**：`modules/__init__.py:80-110` 中 `sc_cell_go` 无依赖声明，`validate_pipeline_order` 不会校验其位于 DEG 之后。
64. **两套 ORA 实现长期并存**（自研 `run_ora_full` vs gseapy）：统计契约不同（见 4.1），需要收敛为一条路径或至少在审计中登记被检验/被排除通路数。
65. **证据分级口径不统一**：`functional_state` 有 `exploratory_n2` 分级，`sc_timecourse`/`neighborhood_da`/`proportion` 没有；建议抽出共享的证据分级 helper。
66. **文件类型白名单未被强制**：`virtual_ko.py:360-362`/`celloracle_worker.py:552` 登记 `file_type='oracle'`，不在 `base.VALID_RESULT_FILE_TYPES`；全仓只有 `base.py:15` 定义与一处测试引用，`worker.py:146/165` 直接落库。当前 DB 尚未出现该类型，属"可达但未发生"。
67. **散点/直方图在 SVG 外的栅格化整体到位**（`base.py:359-365`、`static_rendering.py:216/224`），但 `visualization.py:274-286` 自带一份未净化的 `save_plotly_json`，与 `base.py:526-534` 重复。
68. **CJK 字体缺失静默回退**（`figure_style.py:156-163, 189-212`）：中文标签会渲染成方框且 manifest/结果页无提示。

## 8. P2：低危与可维护性（精选）

- **重复代码**：`annotation.py:3388/3679/3990` 三次逐字重复的 `select_cluster_marker_genes` 参数块；`deg.py` 两张 CSV 的 regulation 覆盖不同基因集。
- **死参数/未注册参数**：`functional_state.py:2150 show_cell_level_correlation`、`annotation.py:3220/3566/3978 merge_similar_threshold/doublet_score_threshold/sample_key`、`sc_cell_go.py:1216 gsea_weight`、`virtual_ko.py` 的 `allow_non_count_fallback/base_grn_version/random_seed`。
- **审计字段与事实相反**：`deg.py:167-169, 758 full_table_exported` 恒为 False 而完整表每次无条件落盘。
- **参数无界**：`schemas.py` 中 `n_comps`、`max_genes`（sc_timecourse）、`cluster_agreement_threshold`、`min_prop`、`n_propagation` 等；`annotation.py:3332/3409-3418/3564-3571` 的 `int()/float()` 无 try。
- **裸 except 静默吞错**：`annotation.py:660-661/686-687`（邻接证据）、`io_utils.py:613-646/683/721`、`visualization.py:23-26`、`static_rendering.py:110-111`、`sc_figure_diagnostics.py:512-527`。
- **绘图细节**：`sc_figure_diagnostics.py:61` 的 `grid=False` 因传入 kwargs 反而打开网格（matplotlib 语义）；`:608/620/724` 静默截断（前 50 PC / Top-20 组）且图上无提示；`figure_style.py:133-143` 只有 10 色而 composition 图最多 30 类 → 取模重复、图例无法唯一映射；`sc_figure_diagnostics.py:252` 先建 figure 再 `return None`（泄漏）；`:77/373-377` donor 数线性放大画布，50 donor × 300 dpi 可达 100–240 MB/张。
- **可视化误导**：`visualization.py:42-46` 抽样 5 万细胞但图内不声明；`static_rendering.py:318-320` 的"暂无可用静态数据"占位图与真实图件无法区分。
- **manifest/报告**：`result_manifest.py:23-35` 只记 path/exists/size；`pipeline_report.py:30-40` 的 JSON 嵌绝对路径（md 报告反而相对化）；`project_report.py:57,192-196` 只收 png、把整图 base64 内联、summary 原样打印绝对路径。
- **测试盲区**：`tests/test_semantic_full.py:279,300` 的静态正则只匹配 `'file_type': '字面量'`，因此 `virtual_ko` 的变量拼装漏检；`tests/test_proportion.py::test_fisher_exact_fallback_to_chi2` 固定了降级行为但未覆盖 summary 误标。
- **`_available_path` 竞态**（`sc_batch.py:373-380`）：exists→创建非原子，当前 1–2 并发概率低。
- **`_safe_name` 实现不一致**（`functional_state.py:162-164` 仅 ASCII vs `scenic.py:30-35` 保留 Unicode）。
- **`expression_parser.py:101-111`**：`padj<0.01` 静默忽略回退默认阈值，`padj<<abc` 的 `ValueError` 逃逸出 `ParseError` 捕获（调用方多处只捕 ParseError）。

## 9. 测试与验证现状

- 本次实际执行：**153 项相关测试全部通过**（50 + 103），无既有回归。
- 覆盖较好的方向：QC/标准化/HVG 语义（`test_qc_semantics`、`test_normalize_semantics`、`test_hvg`）、DEG 契约与分组语义（`test_sc_de_contract`、`test_deg_comparison_scope`、`test_grouping_semantics`）、pseudobulk（`test_sc_batch`）、富集来源绑定（`test_worker_sc_cell_go`）、注释 Universal 面板、SCENIC、functional_state、signature scoring。
- 明显缺口（建议补测试）：
  1. 没有任何测试覆盖 `batch_correct` 的 scvi 分支（P0-1 因此长期未被发现）；
  2. `qc_reassess` 的 doublet 判定语义无断言；
  3. 失败任务的 manifest `status`（P0-4）无测试；
  4. `convert_10x` 二次导入覆盖（P0-12）无测试；
  5. `sc_cluster.score_cell_type_signature` 在不同表达尺度下的行为无测试；
  6. `deg` 的 pairwise 输出是否保留全量细胞无测试；
  7. 重复基因名 h5ad 走完整 QC 无端到端测试；
  8. `design_preflight` 的"表单 vs AI"两条路径一致性无测试。

## 10. 值得肯定的设计

1. **重复单位纪律**：`functional_state._sample_feature_table`（`:724-743`）与 `_feature_statistics`（`:758-800`）全部聚合到 `sample_id × condition × celltype`，n 不足时不给 p；`sc_timecourse._sample_column`（`:304`）主动拒绝"一细胞一 ID"、技术 batch、同一样本跨时间点/跨条件；`proportion._sample_level_composition`、`neighborhood_da._validate_design` 同样拦住把细胞当重复的用法，并把 `inference_unit`、`cell_level_association_p_value` 分列。
2. **pseudobulk 实现正确**：`sc_timecourse._matrix_row_expression`（`:209-216`）先按单元求和再 CPM+log1p；`sc_batch.py:307-370` 只用原始整数 counts 求和并做整数校验；`sc_de_utils.py:104-178` 强制从 `layers['counts']` 重建 log1p 并设 `uns['log1p']`，保证 scanpy 走 `expm1` 得到正确 logFC 尺度。
3. **拒绝静默回退**：`batch_correct` 明确拒绝 DESeq2→Welch 回退（`sc_batch_export.py:856-868`，有回归测试）；`qc.py:646-664/914-931` 拒绝未授权的双细胞 caller 回退；`hvg.py`/`normalize.py` 拒绝在错误表达尺度上运行；`clustering.py:791-800` 回读 `uns['neighbors']['params']['use_rep']` 校验表示一致。
4. **溯源与路径安全**：`config.py:282-294` 的 `_validate_path`（realpath 包含性 + 拒绝符号链接）、`config.py:297-347` 的 `SC_BATCH_SOURCE_ROOTS` 白名单（默认空=关闭）、ZIP 成员名与符号链接位预校验（`convert_10x.py:19-29`）、`gene_set_registry` 的 `schema_version + sha256`、`functional_state` 的 CollecTRI metadata 校验。
5. **任务绑定而非"扫描最新文件"**：DEG 内部表按 task_id 命名，路由 `_sc_deg_enrichment_sources` 按"已完成任务 + 项目内 + 非符号链接"解析，`sc_cell_go` 再交叉校验；worker 用 `task_artifacts/{module}/{task_id}/` 快照，重跑不覆盖旧任务字节。
6. **证据等级写进产物**：`functional_state` 的 n=2 探索层、TF 覆盖度四级策略、ULM/JSD 释义、`pseudobulk_validation` 均固化在 manifest；`sc_figure_diagnostics` 把"donor 不作推断重复""flags 不对应移除数"等边界直接画进图注。
7. **通用边界守卫**：`base.py:1144-1152` 空数据守卫、`subcluster.py:152/162` 与 `clustering.py:156-161` 的维度/细胞数钳制、`qc_reassess.py:153-156` 防止 `auto_remove` 删空全部细胞。

## 11. 建议的修复顺序

- **第一批（阻断，1–2 天）**：P0-1（scVI）、P0-4（manifest 状态）、P0-12（覆盖历史产物）、P0-2（pickle/环境变量）、P0-3（DEG 子集落盘）。
- **第二批（科学正确性，需要领域确认）**：P0-5、P0-6、P0-7、P0-8/P0-9、P0-10、P0-11；随后 4.1（富集 FDR 与 size 过滤）、4.2（DEG 最小样本与 log2FC 口径）、4.3（fillna(0)、n=2 分级、TF 权重）。
- **第三批（工程与安全）**：第 31–43 条（内存、陈旧表示、并发 figure、task_artifacts 保留策略）、第 52–60 条（绝对路径外发、manual_map_csv、dimred 文件名、ZIP 限额）。
- **第四批（一致性收口）**：第 61–68 条（预检入口统一、INPUT_REQUIRES 补齐、两套 ORA 收敛、证据分级共享 helper、白名单强制）。
- **配套测试**：第 9 节列出的 8 个缺口应随对应修复一起补，尤其是 scvi 分支与 manifest 状态。

---

### 附：本次审查的局限

- 全部结论基于静态阅读 + 依赖源码核对 + 少量最小复现（`Index.get_indexer` 行为、scipy 小样本 p 值、gseapy/omicverse 路径核对）；未在真实数据集上执行各分析模块。
- `liana`、`celloracle`、`celltypist` 不在主环境，相关运行期行为依据其安装源码与文档推断（已在正文标注不确定处）。
- 未做性能基准测试；内存/耗时判断基于代码路径与参数上限的静态推断（`task_artifacts` 的 122 GB 为实测）。
- 复核标记仅覆盖本次亲自确认的条目；其余条目保留了子代理给出的行号证据，建议修复前先按行号确认。
