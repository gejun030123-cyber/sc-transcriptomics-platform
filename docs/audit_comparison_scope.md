# 比较组与样本范围（Scope）逻辑审计 — 簇 Marker 及同类分析

> 审计日期：2026-08-27（本次会话）
> 问题类别：**“选不了比较双方”** 与 **“设不了在哪个样品/条件下比较”**，
> 以及同类“指定了比较却在别的口径上仍输出全部组”的逻辑不一致。

## 1. 问题描述

用户报告：**簇 Marker 分析（deg 模块）**
- 无法选择“哪个细胞（簇）与哪些细胞（簇）对比”——此前只能
  “每个簇 vs 其余全部（rest）”或“每个簇 vs 单个参考簇”；
- 无法设置在**哪个样品**中比较——此前全量数据一起检验。

审计发现同类问题还存在于其他模块（见下表），本次一并修复。

## 2. 审计结果（谁涉及这类问题）

| 模块 | 展示名 | 之前能否指定比较对 | 之前能否限定样品/范围 | 问题判定 |
|---|---|---|---|---|
| deg | 簇 Marker（探索性） | 否（仅 rest/单参考） | 否 | **核心问题，本次修复** |
| clustering（marker 预览） | 聚类分析 | 否（cluster-vs-rest 预览） | 否 | 预览性质，保留；如需簇内比较请用 deg/subcluster |
| subcluster | 子簇精细分析 | 否（重聚类后 vs rest） | 否 | **同类问题，本次增加范围限制** |
| proportion | 细胞比例分析 | 细胞级有 compare_groups，但样本级检验忽略它 | 否 | **口径不一致 + 无范围，本次修复** |
| sc_pseudobulk_deg | 样本级 pseudobulk DEG | 有（all_pairwise/vs_reference/comparisons） | 否 | **无范围限制，本次增加** |
| cell_communication | 细胞通讯 | 不适用（全细胞群间） | 否 | **无范围限制，本次增加** |
| trajectory | 轨迹分析 | 不适用（全局轨迹） | 否 | **无范围限制，本次增加** |
| virtual_ko | 虚拟敲除 | 不适用 | 否 | **无范围限制，本次增加** |
| sc_timecourse | 单细胞时序动态 | 有比较/条件概念 | 否 | **无范围限制，本次增加** |
| sc_cell_deg | 细胞级探索性比较 | 有（comparisons/模式） | 有（selected_sample/selected_cluster） | 已是正确范式，作为对照 |
| bulk_deg | Bulk DEG | 有（group1/group2/comparisons） | 不适用（样本即行） | 正常 |
| bulk_heatmap | Bulk 热图 | 有（deg_comparison） | 有（sample_display_mode） | 正常 |
| annotation / qc_reassess / clustering 等结构性步骤 | — | 不适用 | 不适用 | **有意不开放范围**：这些步骤生成全量标签/分群，截断会导致下游注释不完整；样本特异分析请用 subcluster + 范围 |

## 3. 统一方案

### 3.1 统一的“分析范围”参数（新增到受影响模块）

每个受影响模块新增两个参数：

- `scope_key`：obs 列名（如 `sample_id`、`condition`、`timepoint`），留空 = 全部细胞；
- `scope_values`：勾选（动态多选）或逗号/分号/换行填写的取值列表。

实现：`BaseAnalysis.apply_scope(adata)`（modules/base.py）在模块 run() 的最早阶段
按该列取值子集化 AnnData，并校验：
- 列不存在 → 明确报错；
- 填了列但没填值 → 报错；
- 部分取值不存在 → 警告并忽略，全部不存在 → 报错。

### 3.2 deg（簇 Marker）新增比较方式

新增参数：
- `comparison_mode`：`reference`（默认，向后兼容，rest=其余全部）/ `pairwise`（全部两两）/ `custom`（自定义 A-vs-B）；
- `comparisons`：custom 模式的 `A-vs-B` 列表（每行或分号分隔）；
- `min_cells_per_group`：pairwise/custom 模式下每组最少细胞数，不足则跳过并记录。

行为：
- 旧参数（只传 `reference`，无 `comparison_mode`）行为完全不变；
- pairwise/custom 只对**参与比较的组**做检验、绘图、热图和导出（避免“选了比较，
  图里却出现全部组”）；
- 每对比较独立调用 `rank_genes_groups(groups=[A], reference=B)`，并在
  CSV 中新增 `comparison`、`control_group` 两列溯源。

### 3.3 proportion 的一致性修复

- `compare_groups`（A-vs-B）现在**同时**约束细胞级 pairwise 与样本级检验：
  样本级检验只运行指定对（Mann-Whitney），不再无视指定比较输出全部条件 Kruskal；
- 新增 `scope_key/scope_values` 范围限制。

## 4. 涉及文件

- modules/base.py（apply_scope 统一入口）
- modules/deg.py（比较方式 + 范围）
- modules/schemas.py（各模块参数定义）
- modules/proportion.py（范围 + 比较口径一致）
- modules/sc_batch_export.py（sc_pseudobulk_deg 范围）
- modules/cell_communication.py、modules/trajectory.py、modules/virtual_ko.py、
  modules/subcluster.py、modules/sc_timecourse.py（范围）
- modules/design_preflight.py（范围一致性检查 + deg 分组建议）
- modules/ai_tools.py（deg 比较方式建议）
- templates/analysis_select.html（scope_key/scope_values 动态 UI）
- tests/test_deg_comparison_scope.py（新增回归测试）

## 5. 验证

- 新增测试 `tests/test_deg_comparison_scope.py`：10 项通过
  （custom 只跑指定对、pairwise 全两两、旧 reference 兼容、单参考、非法对、
  scope 限制/报错、apply_scope 助手、proportion 样本级比较口径）；
- 既有回归：`tests/test_semantic.py -k deg`（9 通过）、
  `tests/test_integration_full.py -k 'deg_progress or sc_pipeline_qc_to_deg'`（2 通过）、
  `tests/test_design_preflight.py + tests/test_grouping_semantics.py`（17 通过）。
- 注：`tests/test_sc_de_contract.py::test_go_route_binds_the_selected_deg_task_instead_of_scanning_csv`
  在本工作区**改动前即已失败**（已用 git stash 隔离验证），与本次修改无关。
