# 项目整体分析与分步增强路线图

日期：2026-07-08
状态：阶段 1-6 已完成

## 1. 总体判断

当前项目已经具备很宽的能力面：Flask Web 工作台、单细胞与 Bulk RNA-seq 模块、结构化参数 schema、异步任务、Plotly 结果展示、AI 对话工具、PipelineRun、Agent 会话、候选分支和评分模型都已经存在。

下一阶段最值得投入的方向不是继续堆新模块，而是把已有能力收束成一个可信的分析闭环：

```text
用户目标
  -> 当前状态检查
  -> 参数候选/分支实验
  -> 自动评分 + 可视化证据
  -> 人工审批
  -> 采纳为后续分析基线
  -> 导出可复现结果包
```

换句话说，项目最需要补的是“结果是否可信、为什么可信、如何继续调参、如何可复现”，而不是单纯多加一个算法入口。

## 2. 最希望修改和增加的内容

### P0：结果解释与结果包完善

这是最优先内容。当前模块已经能产生大量 h5ad、CSV 和 Plotly JSON，但结果层仍偏“文件列表”，缺少面向生信判断的统一摘要。

建议增加：

- 每次任务生成一个统一的 `analysis_manifest.json`，记录输入、输出、参数、软件版本、关键指标、警告、可继续分析的下游模块。
- 每个项目生成一个 `project_report.md` 或 `project_report.html`，自动汇总已完成步骤、关键图、关键表、异常告警、建议下一步。
- Plotly JSON 之外，提供可直接浏览的离线 HTML 图库，避免用户只能看到原始 JSON。
- 为单细胞关键阶段生成“审批证据卡”：QC 是否过严、UMAP 是否受 batch 支配、cluster 是否过碎、annotation 是否低置信、DEG 是否足够支持 marker。
- 为 Bulk 关键阶段生成“审批证据卡”：样本是否离群、分组是否分开、差异基因是否方向一致、富集结果是否由少数基因驱动。

涉及文件：

- `worker.py`
- `modules/base.py`
- `modules/visualization.py`
- `routes/results.py`
- `templates/analysis_result.html`
- `templates/results_gallery.html`
- 新增 `modules/reporting/result_manifest.py`
- 新增 `modules/reporting.py` 或 `modules/reporting/` 内的报告生成工具

### P1：目标驱动 Agent 闭环补强

项目已有 `modules/agent_orchestrator.py`、`routes/branches.py`、`AgentSession`、`AnalysisBranch`、`CandidateScore` 等基础，但还需要把它从“能运行候选”升级为“能解释候选、比较候选、收敛目标”。

建议增加：

- 候选分支对比视图：参数 diff、评分 diff、关键图并排、警告并排。
- Agent 每一步输出“可审批摘要”，包括本轮目标、使用证据、推荐动作、风险、是否需要确认。
- `continue_goal_agent` 支持更明确的 refine 语义：更细分、更保守、去 batch、保留罕见群、提高 marker 特异性、降低 doublet 风险。
- 评分模型从单一 marker score 扩展为多证据评分：marker、cluster size、batch composition、QC、annotation confidence、DEG support。
- 采纳候选后自动生成“采纳记录”：为什么采纳、替代候选为什么未采纳、后续从哪个 h5ad 继续。

涉及文件：

- `modules/agent_orchestrator.py`
- `modules/ai_tools.py`
- `modules/evaluators/sc_cluster.py`
- `modules/cell_markers.py`
- `routes/branches.py`
- `templates/project_detail.html`
- `tests/test_agent_orchestrator.py`
- `tests/test_branch_routes.py`
- `tests/test_signature_scoring.py`

### P2：PipelineRun 与单任务结果打通

当前 `PipelineRun` 能串联多个模块，但结果查看、失败恢复、局部重跑和 report 聚合仍应加强。否则用户一旦跑完整流程，仍需要逐个任务翻结果。

建议增加：

- PipelineRun 详情页：显示每个模块状态、耗时、输出、失败点、可从失败模块继续。
- PipelineRun 聚合 report：一次完整 sc 或 bulk 运行后自动汇总结果。
- 支持从任意 completed 子任务继续跑后续模块。
- 支持“同一输入 + 不同参数”的 pipeline 级候选运行，服务于 Agent 分支搜索。

涉及文件：

- `worker.py`
- `routes/api.py`
- `routes/results.py`
- `models.py`
- `templates/project_detail.html`
- 新增 `templates/pipeline_run_detail.html`
- `tests/test_pipeline.py`

### P3：参数 schema 与动态推荐增强

参数 schema 已经很完整，但还可以更像“分析助手”而不只是表单。

建议增加：

- 基于数据规模、物种、batch 列、样本数自动推荐默认参数。
- schema 支持参数风险提示：例如 QC 阈值过严、cluster resolution 过高、DEG 组别样本太少。
- schema 支持结果驱动建议：例如 UMAP batch 混合差时建议 batch_correct，annotation 低置信时建议调整 marker 或重分群。
- 增加“参数预设对比”：保存、加载、比较项目级 preset。

涉及文件：

- `modules/schemas.py`
- `routes/api.py`
- `templates/analysis_select.html`
- `templates/sc_analysis.html`
- `templates/bulk_analysis.html`
- `tests/test_schemas.py`

### P4：缺失依赖、失败恢复与可观测性

项目依赖较多，且 README 已说明没有锁定的依赖文件。为了真正可交付，需要降低环境不一致带来的失败。

建议增加：

- `requirements.txt` 或 `pyproject.toml`，区分 core、sc、bulk、ai、optional extras。
- 启动时 dependency health check，展示哪些模块可用、哪些因可选依赖缺失不可用。
- 任务失败分类：输入错误、依赖缺失、参数错误、算法运行错误、系统资源不足。
- 日志中记录内存峰值、运行耗时、输入规模，便于判断大数据运行风险。

涉及文件：

- 新增 `requirements.txt` 或 `pyproject.toml`
- `app.py`
- `routes/api.py`
- `worker.py`
- `modules/base.py`
- `templates/index.html`
- `tests/test_integration.py`

## 3. 分步执行计划

### 阶段 1：项目级结果 manifest 与报告骨架

目标：让每次分析的产物更可审计、更容易汇总。

任务：

- 新增 `modules/reporting/result_manifest.py`，定义统一 manifest 数据结构。
- 在 `worker.py` 的单任务和 PipelineRun 执行结束时写入 manifest。
- manifest 至少包含：project_id、task_id、module_name、params、input_path、output_adata、result_files、summary、warnings、created_at。
- 在 `routes/results.py` 或 API 中暴露 manifest。
- 增加测试，确认每个任务完成后 manifest 可生成且路径安全。

验收：

- 跑任意模块后，项目 `results/` 或任务输出目录中可看到 manifest。
- manifest 能指向真实存在的 h5ad、CSV、Plotly JSON。
- 不改变现有 ResultFile 表和旧结果页行为。

### 阶段 2：离线 HTML 图库与项目 report

目标：让用户直接看到结果，不只看到文件名或 JSON。

任务：

- 新增 `modules/reporting.py` 或 `modules/reporting/html_report.py`。
- 将 Plotly JSON 转换为离线 HTML 图。
- 生成项目级 `project_report.md`，汇总任务顺序、关键指标、图表链接和下一步建议。
- 在结果页增加“查看报告 / 下载报告”入口。

验收：

- 对 PBMC3k 或已有 Bulk reference 输出，生成可直接打开的 HTML 图库。
- report 中的文件链接均存在。
- 大图表数量较多时页面不会明显卡死。

### 阶段 3：单细胞审批证据卡

目标：把 QC、clustering、annotation、DEG 的结果变成可判断的证据。

任务：

- 为 `qc` 输出过滤强度、细胞保留率、MT/ribo/hb 风险、doublet 风险。
- 为 `clustering` 输出 cluster 数、最小簇比例、batch 支配簇、分辨率拆分证据。
- 为 `annotation` 输出低置信 cluster、marker score margin、Unknown 比例。
- 为 `deg` 输出每簇支持 marker 数、显著 DEG 数、marker 热图可用性。
- 将证据卡写入 manifest 和 report。

验收：

- 每个关键模块的 summary 里有结构化 `review_evidence`。
- 前端可展示“通过 / 警告 / 需复核”状态。
- 不把警告当失败；只影响人工审批建议。

### 阶段 4：Agent 候选分支对比与采纳记录

目标：让目标驱动分析真正闭环。

任务：

- 扩展 `CandidateScore` 或新增候选摘要字段，保存多维评分。
- 在 `routes/branches.py` 增加候选对比 API。
- 在 `templates/project_detail.html` 的 Agent 面板中展示候选 diff。
- 采纳分支时写入 `accepted_branch_report.json` 和 Markdown 摘要。
- `continue_goal_agent` 基于用户反馈生成下一轮搜索约束。

验收：

- 同一目标下多个候选可以并排比较。
- 用户能看出推荐候选为什么优于其他候选。
- 采纳后 `/current-context` 明确指向已采纳分支。

### 阶段 5：PipelineRun 结果聚合与失败续跑

目标：完整流程不再只是后台串行任务，而是一个可复核的实验运行。

任务：

- 新增 PipelineRun 详情页。
- 将子任务 manifest 聚合为 pipeline manifest。
- 失败时记录失败模块、输入文件、参数和可续跑建议。
- 支持从失败模块或任意 completed 子任务继续跑后续模块。

验收：

- PipelineRun 页面能看到每一步状态和结果链接。
- 失败后可定位具体模块和错误类型。
- 续跑不会覆盖已完成任务的结果。

### 阶段 6：依赖与健康检查

目标：降低部署和运行失败成本。

任务：

- 增加依赖声明文件。
- 增加 `/api/system/dependencies`，返回核心依赖和可选依赖状态。
- 首页或系统状态面板展示模块可用性。
- worker 捕获 ImportError 时返回更友好的缺失依赖提示。

验收：

- 新环境能按依赖文件安装核心功能。
- 未安装可选依赖时，对应模块显示“不可用/需安装”，而不是运行到中途才失败。

## 4. 推荐实施顺序

优先顺序：

1. 阶段 1：manifest。它是后续 report、审批证据、Agent 解释的共同底座。
2. 阶段 2：HTML 图库和项目 report。马上改善结果可读性。
3. 阶段 3：单细胞审批证据卡。直接服务“结果完善”。
4. 阶段 4：Agent 候选对比。直接服务“功能增加”和目标驱动调参。
5. 阶段 5：PipelineRun 聚合。增强完整流程执行体验。
6. 阶段 6：依赖健康检查。降低部署和长期维护风险。

## 5. 当前风险与注意事项

- 当前工作区已有大量未提交修改，实施前需要避免覆盖既有变更。
- `docs/codebase-map.md` 中测试数量描述已经落后于当前仓库，应在后续文档更新阶段同步修正。
- `tests/test_semantic_full.py` 会检查模块静态语义，新 helper 或新模块需要遵守既有模块形态约定。
- 结果完善不应只增加前端展示；必须落到可审计的 JSON/Markdown/HTML 文件。
- Agent 的自动评分只能作为推荐依据，不能替代人工生物学判断。

## 6. 第一阶段建议修改清单

如果审批通过，第一阶段建议只做小而稳的增量：

- 新增 `modules/reporting/result_manifest.py`
- 修改 `worker.py`
- 修改 `routes/results.py` 或新增轻量 API
- 修改 `templates/analysis_result.html`
- 新增 `tests/test_result_manifest.py`

第一阶段不改各个分析模块内部算法，只统一收集它们已有的 `summary`、`result_files` 和 `output_adata`，因此风险最低，也最利于后续所有增强。

## 7. 完成记录

已完成阶段：

- 阶段 1：任务级 `analysis_manifest.json`，接入单任务和 PipelineRun 子任务。
- 阶段 2：项目级 `project_report.md` 和离线 Plotly HTML 图库。
- 阶段 3：单细胞 QC/clustering/annotation/DEG 审批证据卡。
- 阶段 4：Agent 候选分支对比 API、前端摘要、采纳记录 JSON/Markdown。
- 阶段 5：PipelineRun 详情页、pipeline manifest/report、失败续跑 API。
- 阶段 6：`requirements.txt`、`/api/system/dependencies` 依赖健康检查。

新增主要文件：

- `modules/reporting/result_manifest.py`
- `modules/reporting/project_report.py`
- `modules/reporting/review_evidence.py`
- `modules/reporting/branch_report.py`
- `modules/reporting/pipeline_report.py`
- `modules/platform/system_health.py`
- `templates/pipeline_run_detail.html`
- `requirements.txt`

验证命令：

```bash
python -m py_compile modules/reporting/result_manifest.py modules/reporting/project_report.py modules/reporting/review_evidence.py modules/reporting/branch_report.py modules/reporting/pipeline_report.py modules/platform/system_health.py worker.py routes/api.py routes/results.py routes/branches.py models.py
pytest tests/test_review_evidence.py tests/test_result_manifest.py tests/test_project_report.py tests/test_branch_routes.py tests/test_pipeline_report.py tests/test_system_health.py tests/test_pipeline.py tests/test_models.py -q
```
