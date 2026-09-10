# 单细胞分析审查问题修复对照说明

- 日期：2026-09-10
- 对照报告：`docs/single-cell-code-review-2026-09-10.md`
- 工作方式：只读核对 + 逐条修复 + 每批修复后运行相关 pytest

## 1. 本次修复的实际情况

进入修复时，工作区中**已经存在一批未提交的修复**（覆盖审查报告 P0 的 9 项与多条安全项，另有 9 个测试文件被同步补充）。本轮工作分三步：

1. **检验现状**：逐条核对已有修复是否真正解决了对应问题（读 diff、重跑测试），确认没有仅"看起来改了"的情况；
2. **补齐剩余 P0**：审查报告的 12 项 P0 中，还有 3 项未被处理（`sc_cluster` 评分、LLM 截断、零证据簇复核），本轮全部修复；
3. **追加高价值 P1**：优先修复会造成崩溃、内存翻倍、并发干扰、隐私外发与入口不一致的问题。

所有结论均以代码与实际测试为准；下方每一项都标注了对应测试。

## 2. P0 修复对照（12/12）

| 编号 | 问题 | 修复位置与方式 | 验证 |
| --- | --- | --- | --- |
| P0-1 | scVI 分支引用未定义变量 `layer`，100% 失败 | `modules/batch_correct.py`：引入 `training_layer = 'counts'`，返回值改用该变量 | `test_batch_correct_evaluation.py::test_scvi_records_the_actual_training_layer` |
| P0-2 | 上传 pickle 在继承全部环境变量的子进程中反序列化，口令/密钥可被读出 | `modules/virtual_ko.py`：子进程环境白名单（剔除 `*_KEY/*_TOKEN/*PASSWORD`）；`modules/celloracle_worker.py`：禁用 pickle/Oracle/Links 载入，base GRN 仅接受表格；`routes/upload.py`/`routes/analysis.py`/`schemas.py`/模板同步收敛扩展名 | `test_virtual_ko.py` |
| P0-3 | pairwise DEG 把比较子集写成标准输出，下游静默丢细胞 | `modules/deg.py`：`output_adata` 保留完整 scoped 细胞，summary 记录 `n_cells_in_output/used_for_contrasts/excluded` | `test_deg_comparison_scope.py` |
| P0-4 | manifest `status` 硬编码 `completed`，失败任务留"已完成"审计 | `modules/reporting/result_manifest.py`：新增 `_result_status()` 推导状态并写入 `error` | `test_result_manifest.py::test_manifest_records_module_declared_failure` |
| P0-5 | doublet 阈值混用三种不同定义的量，列名却是 `doublet_fraction` | `modules/qc_reassess.py`：拆出 `doublet_fraction` / `mean_doublet_score` / `doublet_fraction_source`，新增 `doublet_score_cutoff` 参数 | `test_qc_semantics.py::test_qc_reassess_does_not_treat_mean_score_as_doublet_fraction` |
| P0-6 | 逐细胞（非独立）Spearman 的 p 值写进 CSV，与 manifest 声明矛盾 | `modules/functional_state.py`：`cell_level_descriptive` 层级 `pvalue = NaN` | `test_functional_state.py::test_cell_level_correlation_omits_nonindependent_pvalue` |
| P0-7 | 组成比例截断 Top-N 后重新归一化分母 | `modules/sc_figure_diagnostics.py`：补 `Other` 段，分母保持全部细胞 | `test_sc_figure_diagnostics.py` |
| P0-8 | marker 评分不校验表达尺度，raw counts 下 `0.35×均值` 直接饱和 | `modules/evaluators/sc_cluster.py`：显式解析尺度（counts → normalize+log1p；含负值 → 拒绝并报错），`total_score` 只由有界、尺度无关的 marker 证据构成 | `test_signature_scoring.py`（新增 3 例） |
| P0-9 | 评分逐基因 `except Exception: pass`，失败被当成 0 分 | `modules/evaluators/sc_cluster.py`：改为一次性取 marker 表达矩阵，任一步失败返回 `error` | 同上 |
| P0-10 | `max_tokens=4096` 与 30 簇中文输出结构性失配，截断即整任务失败 | `modules/llm_annotation.py`：按批大小动态预算（2048–16000）、响应契约失败自动减半重试并记录 `batch_splits`、Anthropic 分支补 `temperature=0` | `test_llm_annotation.py`（新增 4 例） |
| P0-11 | 零证据簇仅靠逐细胞投票恢复标签，却不标 `needs_review` | `modules/annotation.py`：`unknown_recovered_by_cell_vote` 强制 `needs_review=True` | `test_annotation_universal.py::test_cluster_review_flags_zero_evidence_cell_vote_recovery` |
| P0-12 | `convert_10x` 固定输出名，二次导入覆盖历史任务产物 | `modules/convert_10x.py`：三处输出改用 `_available_path()`；批次名在清洗后校验唯一性 | `test_sc_matrix_import.py` |

## 3. 本轮追加修复的 P1

| 编号 | 问题 | 修复 | 验证 |
| --- | --- | --- | --- |
| #31 | worker 预检 `adata` 存活到 `run()` 结束，峰值内存约 2× | `worker.py`：验证后 `del adata` | 全量回归 |
| #32 | 重复基因名令 `qc.py` 在分析末尾抛 `InvalidIndexError` | `modules/base.py`：`load_adata` 统一保证 `var_names` 唯一；`modules/qc.py`：回填前加唯一性守卫并给出提示 | `test_base.py`、全量回归 |
| #33 | QC 散点图三列独立过滤后按下标配对，错位且可能 ValueError | `modules/qc.py`：改为一次性构表 `dropna` 后取列（counts/genes/MT% 与 novelty 两处） | `test_qc_semantics.py` |
| #34 | 全基因 PCA 把稠密 `layers['scaled']` 写进每个中间 h5ad | `modules/dimred.py`：PCA 完成后删除该层并在 summary 记录 | `test_integration_full.py`（断言扩展） |
| #40 | worker 线程内 `plt.close('all')` 会关掉并发任务的画布 | `modules/annotation.py`：改为 `plt.close(fig_dotplot)` | 全量回归（无新增告警） |
| #46 | `n_propagation` 无界、非 TF 扰动基因只在 GRN 推断后才失败 | `modules/schemas.py` 加 1–5 界；`virtual_ko.py` 子进程启动前校验；`celloracle_worker.py` 在 GRN 推断前完成 TF/基因名校验，并修正误导性的 counts 回退日志 | `test_virtual_ko.py`、全量回归 |
| #52 | 绝对服务器路径随 `result_json` 回流给外部 AI | `modules/ai_tools.py`：`execute_tool` 统一脱敏（项目内 → 项目相对路径，服务器根目录 → 仅文件名），`_validate_project_path` 支持相对路径以保持工具可用 | `test_agent_tools.py`（新增 2 例） |
| #61 | 网页表单与 AI 路径的硬预检范围不一致（人工提交反而更宽松） | `modules/design_preflight.py` 新增共享 `PREFLIGHT_MODULES`，`routes/analysis.py` 与 `ai_tools.py` 共用 | `test_design_preflight.py`、`test_input_selection.py` |
| #62 | `sc_cell_deg`/`sc_pseudobulk_deg`/`functional_state` 未声明 `INPUT_REQUIRES`，输入下拉过度宽松 | 三个模块补 `['leiden']` / `['celltype']` | `test_input_selection.py` |
| 附加 | `save_matplotlib_figure` 在静态格式被全部取消时提前返回却不关闭画布 | `modules/base.py`：该分支关闭画布 | `test_base.py` |
| 附加 | 失败/不可用任务仍把上游 h5ad 登记为输出（存储浪费 + 误读） | `worker.py`：`_result_error` 为真时不复制 `output_adata` | 全量回归 |
| 附加 | `manual_map_csv` 可读取任意服务器路径 | `modules/annotation.py`：限定项目目录内、拒绝符号链接、要求 `.csv` | `test_annotation_refined_v2.py` |
| 附加 | `dimred` 用未净化的 obs 列名拼文件名，可越出 `plots/` | `modules/dimred.py`：文件名 token 清洗 | `test_dimred_filename_security.py` |
| 附加 | virtual_ko base GRN 白名单覆盖整个 `DATA_DIR`，可跨项目引用他人 pickle | `modules/virtual_ko.py`：收敛为当前项目 uploads + 表格格式白名单 | `test_virtual_ko.py` |

## 4. 本轮未修复（需要单独决策或更大改动）

审查报告中的以下条目**本轮有意未改动**，原因分三类：属于行为/设计变更需产品确认、会牵动跨模块流程、或需要真实大数据验证。

### 4.1 需要领域/产品确认（行为会变化）
- **#35 陈旧校正表示未清理**：重跑 dimred 后 `X_pca_harmony` 等仍取自旧 PCA。清理正确但会让用户必须重跑 `batch_correct`（scVI 可能数十分钟），需要确认交互预期。
- **#36 用 PCA 肘部维度截断 scVI/SysVI 潜空间**：需要确认不同表示的维度选择策略。
- **#37 marker 预览固定用首个分辨率**：主分辨率与预览不一致，修复会改变现有图与 `uns['marker_selection']` 的簇编号。
- **#4.2 默认 `min_samples_per_condition=2`**（proportion / neighborhood_da 在 n=2 时永远无法显著）：提高默认值会改变既有结论门槛。
- **4.1 富集双实现的 FDR 族统一**（cell_level 走 gseapy、pseudobulk 走 `run_ora_full`）：统一会改变历史结果数值，建议先确认课题口径。

### 4.2 较大的工程改动
- **#42 `task_artifacts` 无保留策略**：实测单项目已 122 GB。需要设计配额/清理/去重策略（不是单文件改动）。
- **#41 annotation logreg 每簇全量 copy**、**#45 celloracle Δ 矩阵稠密化**：需要重写热点循环并用真实数据做内存基准。
- **#56 ZIP 解压无体积上限 / 全平台无 `MAX_CONTENT_LENGTH`**：需要确定部署级上传上限。

### 4.3 安全/隐私的剩余面
- **#54 LLM 组织背景的标识符过滤仍偏弱**（仅拦邮箱与少数路径前缀）：建议改为受控枚举或白名单校验。
- **#59 `sc_batch` 把服务器绝对路径写入 obs/summary**：AI 侧已在出口脱敏（#52），但源头仍建议改为相对路径。
- **#60 `enrichr` 外传缺少审计字段**：建议在 audit 中登记 `external_transfer{service,endpoint,n_genes}`。

### 4.4 观察项
- matplotlib 在整份 `test_annotation_universal.py` 同进程运行时会出现 "More than 20 figures have been opened" 告警（单测隔离运行不触发）。`plt.close('all')` 此前部分掩盖了该现象；现已移除全局关闭，建议后续排查 `sc_figure_diagnostics`/`native_figures` 辅助函数的画布生命周期（P2）。

## 5. 测试结果

| 测试批次 | 结果 |
| --- | --- |
| 第一批已修复项对应测试（10 个文件） | 107 passed |
| 单细胞回归集（23 个文件） | 310 passed |
| dimred 集成（含 #34 新断言） | 4 passed |
| AI/安全/上传/schema 回归 | 126 passed（脱敏误伤 API URL 后修复并复验 45 passed） |
| `test_signature_scoring.py`（含新增尺度无关断言） | 16 passed |
| `TestAllModulesInputRequires` + DE 契约 + 输入选择 | 15 passed |
| 全量 `pytest tests/`（第一次） | 1072 passed / 8 skipped / **1 failed**（`TestAllModulesInputRequires::test_requires_matches_validate`：`sc_cell_deg` 声明了 `INPUT_REQUIRES=['leiden']` 却未在 `validate_input`/`run` 中体现） |
| 全量 `pytest tests/`（修复后终态） | **1074 passed / 8 skipped / 0 failed**（22 分 17 秒） |

**最终结果**：修复 `sc_cell_deg.validate_input`（补上分组列与 counts 层的显式前置校验，使声明与实现一致）后，全量套件复跑全绿。

> 说明：第一次全量运行暴露的这 1 个失败本身就是本轮 #62 改动的回归防护，已按"声明必须与实现一致"的原则修正（补 `validate_input`，而不是放宽测试），因此 `INPUT_REQUIRES` 既过滤了输入下拉，也在 worker 预检阶段给出可操作错误。


## 6. 变更文件（本轮）

- `modules/evaluators/sc_cluster.py`（P0-8/9，重写评分内核）
- `modules/llm_annotation.py`（P0-10）
- `modules/annotation.py`（P0-11、#40）
- `modules/qc.py`（#32、#33）
- `modules/base.py`（#32、画布关闭）
- `modules/dimred.py`（#34）
- `modules/celloracle_worker.py`、`modules/virtual_ko.py`、`modules/schemas.py`（#46）
- `modules/ai_tools.py`（#52）
- `modules/design_preflight.py`、`routes/analysis.py`（#61）
- `modules/sc_cell_deg.py`、`modules/sc_batch_export.py`、`modules/functional_state.py`（#62）
- `worker.py`（#31、失败任务不登记输出）
- 测试：`tests/test_signature_scoring.py`、`tests/test_llm_annotation.py`、`tests/test_annotation_universal.py`、`tests/test_agent_tools.py`、`tests/test_integration_full.py`
