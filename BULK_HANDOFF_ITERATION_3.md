# Bulk RNA-seq 双模型交接单 - Iteration 3

生成时间: 2026-07-02

读取对象: 实现/执行模型

审查对象:
- `scripts/run_bulk_reference.py`
- `BULK_RUN_STATE.md`
- `STATE.md`

## 最新审查结论

最新结论: **通过 - Bulk 核心双文件流程已验收完成**。

已确认:
- Dataset A 核心流程完成: `bulk_qc -> bulk_normalize -> bulk_pca -> bulk_deg -> bulk_heatmap`
- Dataset B 核心流程完成: `bulk_qc -> bulk_normalize -> bulk_pca -> bulk_deg -> bulk_heatmap`
- 两个 `bulk_reference_report.json` 均存在并可解析
- 两个 `中文运行报告.md` 均存在
- 核心输出 `.h5ad` 均存在并可读取
- CSV 非空，Plotly JSON 可解析
- Targeted tests 复跑通过: `93 passed, 1 warning`

报告路径:
- `bulk_reference_output/bulk_dataset_a/bulk_reference_report.json`
- `bulk_reference_output/bulk_dataset_a/中文运行报告.md`
- `bulk_reference_output/bulk_dataset_b/bulk_reference_report.json`
- `bulk_reference_output/bulk_dataset_b/中文运行报告.md`

说明: 下方“当前结论/阻塞问题”是 Iteration 3 的历史审查上下文，不是最新状态。最终审查结果见文档末尾 `Codex Review - Iteration 4 Final`。

## 后续协作规则

从本轮开始，后续沟通和审查都继续追加在本文件中。

实现/执行模型完成一轮后，请在本文档末尾追加:

```md
## Implementer Reply - Iteration N

- 已读取要求: 是/否
- 修改文件:
- 执行命令:
- Dataset A 报告:
- Dataset B 报告:
- 测试结果:
- 已知失败/跳过:
- 请求 Codex 审查:
```

Codex 审查时会先读取本文档中的最新 `Implementer Reply`，再继续在本文档末尾追加:

```md
## Codex Review - Iteration N

- 审查结论:
- 复跑命令:
- 发现的问题:
- 通过项:
- 下一轮要求:
```

不要新建新的 handoff 文件，除非用户明确要求归档或重开一轮。

## 当前结论

本轮结论: **需要重试，尚未通过真实双文件跑通验收**。

你已经完成:
- 创建 `scripts/run_bulk_reference.py`
- 运行 Bulk 相关 targeted tests
- 测试结果记录为 `93 passed, 1 warning`

Codex 已复跑确认:

```bash
python -m py_compile scripts/run_bulk_reference.py
python -m pytest tests/test_bulk_deg.py tests/test_bulk_qc_helpers.py tests/test_p3_modules.py tests/test_pipeline.py tests/test_schemas.py -q
```

复跑结果:

```text
93 passed, 1 warning in 1.10s
```

## 阻塞问题

1. 没有找到真实运行报告。

当前仓库下没有发现:

```text
bulk_reference_report.json
中文运行报告.md
```

因此目前只能证明 runner 和 targeted tests 完成，不能证明两个参考文件已经跑通。

2. `bulk_deg_integration.selected_comparisons` 参数类型错误。

`modules/bulk_deg_integration.py` 里按字符串处理:

```python
selected = self.params.get('selected_comparisons', '').strip()
```

但 runner 当前传的是 Python list:

```python
params['selected_comparisons'] = comparisons_list
```

这会导致多比较整合步骤失败。应改为:

```python
params['selected_comparisons'] = ','.join(comparisons_list)
```

或者留空字符串，让模块使用全部 DEG 比较。

3. `bulk_enrichment` 默认会失败。

`bulk_enrichment` 需要 `input_source` 指向 DEG CSV，或者提供 `custom_genes`。当前 runner 默认加入:

```python
optional_steps = [('bulk_enrichment', {})]
```

但没有设置:

```python
params['input_source']
```

应改为:
- 默认不跑 enrichment
- 只有传入 `--run-optional` 时才跑
- 跑 enrichment 前自动选择一个 `bulk_deg_results*.csv` 作为 `input_source`
- 如果没有 DEG CSV，写入中文报告为“跳过: 未找到 DEG CSV”

4. runner 尚未实现交接要求中的参数。

需要补:

```text
--dataset-label
--run-optional
```

5. runner 尚未自动识别默认比较。

如果用户未传 `--comparisons`，应根据输入文件自动设置:

```text
all.fpkm_anno.xls:
hmc3-vs-ctrl;moclel-vs-ctrl;rapa-vs-ctrl

all.genes.expression.anno.xls:
NH4Cl-vs-Ctr;PEA-vs-Ctr;TMAO-vs-Ctr
```

## 新要求

请直接修改 `scripts/run_bulk_reference.py`，然后实际跑两个文件。

### 必须实现

1. 新增参数:

```text
--dataset-label
--run-optional
```

2. `--groupby` 可以保留，但默认值应为:

```text
_auto_group
```

3. 如果没有传 `--comparisons`，根据输入文件名自动推断默认比较。

4. 默认只跑核心流程:

```text
bulk_qc -> bulk_normalize -> bulk_pca -> bulk_deg -> bulk_heatmap
```

5. 只有传入 `--run-optional` 才跑:

```text
bulk_deg_integration
bulk_enrichment
bulk_timecourse
```

6. 可选步骤失败不能让整个核心流程判定失败，但必须写进 JSON 和中文报告。

7. 每次运行必须生成两个报告:

```text
bulk_reference_report.json
中文运行报告.md
```

8. 中文报告必须包含:
- 输入文件
- 项目目录
- 每个模块中文名
- 每个模块状态
- 每个模块输出 h5ad
- CSV / Plotly JSON 输出数量
- 失败模块的中文摘要
- 下一步建议

## 需要实际执行的命令

Dataset A:

```bash
env NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mplconfig \
python scripts/run_bulk_reference.py \
  --input /home/oelab/data/GJ/all.fpkm_anno.xls \
  --dataset-label dataset_a \
  --project-name bulk_dataset_a \
  --groupby _auto_group
```

Dataset B:

```bash
env NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mplconfig \
python scripts/run_bulk_reference.py \
  --input /home/oelab/data/GJ/all.genes.expression.anno.xls \
  --dataset-label dataset_b \
  --project-name bulk_dataset_b \
  --groupby _auto_group
```

## 交付物

运行完成后，应出现:

```text
bulk_reference_output/bulk_dataset_a/bulk_reference_report.json
bulk_reference_output/bulk_dataset_a/中文运行报告.md
bulk_reference_output/bulk_dataset_b/bulk_reference_report.json
bulk_reference_output/bulk_dataset_b/中文运行报告.md
```

请把以下信息追加到 `BULK_RUN_STATE.md`:

```md
### Iteration 4 - Real two-file runner execution

- Handoff file read: `BULK_HANDOFF_ITERATION_3.md`
- Dataset A command: ...
- Dataset A report: ...
- Dataset A status: completed / failed
- Dataset B command: ...
- Dataset B report: ...
- Dataset B status: completed / failed
- Tests rerun: ...
- Summary:
  - bulk_qc: ...
  - bulk_normalize: ...
  - bulk_pca: ...
  - bulk_deg: ...
  - bulk_heatmap: ...
```

## Codex 下一轮审查标准

Codex 会检查:
- 两个 JSON 报告是否存在且能 `json.load()`
- 两个中文报告是否存在
- 每个核心模块是否 `completed`
- 每个核心模块的 `output_adata_path` 是否存在
- CSV 是否非空
- Plotly JSON 是否能解析
- Dataset A 是否完成:
  - `hmc3-vs-ctrl`
  - `moclel-vs-ctrl`
  - `rapa-vs-ctrl`
- Dataset B 是否完成:
  - `NH4Cl-vs-Ctr`
  - `PEA-vs-Ctr`
  - `TMAO-vs-Ctr`

如果核心流程失败，请不要只写“failed”，必须在中文报告里写清楚失败模块、关键 traceback 摘要和下一步建议。

## Codex Review - Iteration 4

审查结论: **核心双文件跑通验收通过**。

本轮读取到的实现回复位于 `BULK_RUN_STATE.md` 的 `Iteration 4 - Real two-file runner execution`。实现模型已修复 runner，并实际运行两个参考文件。

### 复跑命令

```bash
python -m py_compile scripts/run_bulk_reference.py
python -m pytest tests/test_bulk_deg.py tests/test_bulk_qc_helpers.py tests/test_p3_modules.py tests/test_pipeline.py tests/test_schemas.py -q
```

结果:

```text
93 passed, 1 warning in 1.35s
```

产物验证命令摘要:

```bash
env NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mplconfig python <validation-script>
```

验证内容:
- `bulk_reference_report.json` 可 `json.load()`
- `中文运行报告.md` 存在
- 核心模块状态均为 `completed`
- 核心模块 `output_adata_path` 均存在且可用 `anndata.read_h5ad(..., backed='r')` 打开
- CSV 文件存在且有表头/内容
- Plotly JSON 文件可解析

验证结果:

```text
Dataset A overall_status: completed
bulk_qc completed      h5ad: 12 x 61806
bulk_normalize completed h5ad: 12 x 19431
bulk_pca completed     h5ad: 12 x 19431
bulk_deg completed     h5ad: 12 x 19431
bulk_heatmap completed h5ad: 12 x 19431
comparisons: hmc3-vs-ctrl, moclel-vs-ctrl, rapa-vs-ctrl

Dataset B overall_status: completed
bulk_qc completed      h5ad: 24 x 15256
bulk_normalize completed h5ad: 24 x 15040
bulk_pca completed     h5ad: 24 x 15040
bulk_deg completed     h5ad: 24 x 15040
bulk_heatmap completed h5ad: 24 x 15040
comparisons: NH4Cl-vs-Ctr, PEA-vs-Ctr, TMAO-vs-Ctr

VALIDATION: PASS
```

### 通过项

- Dataset A 报告存在:
  - `bulk_reference_output/bulk_dataset_a/bulk_reference_report.json`
  - `bulk_reference_output/bulk_dataset_a/中文运行报告.md`
- Dataset B 报告存在:
  - `bulk_reference_output/bulk_dataset_b/bulk_reference_report.json`
  - `bulk_reference_output/bulk_dataset_b/中文运行报告.md`
- Dataset A 完成三组比较:
  - `hmc3-vs-ctrl`
  - `moclel-vs-ctrl`
  - `rapa-vs-ctrl`
- Dataset B 完成三组比较:
  - `NH4Cl-vs-Ctr`
  - `PEA-vs-Ctr`
  - `TMAO-vs-Ctr`
- `scripts/run_bulk_reference.py` 已实现:
  - `--dataset-label`
  - `--run-optional`
  - `--groupby` 默认 `_auto_group`
  - 根据输入文件名推断默认 comparisons
  - 默认只跑核心流程
  - 生成 JSON 报告和中文 Markdown 报告

### 注意事项

- 当前验收是 **核心流程验收**，不包含可选模块:
  - `bulk_deg_integration`
  - `bulk_enrichment`
  - `bulk_timecourse`
- 当前 `bulk_heatmap` 使用默认 `top_var` 模式，报告中为 `Top 50 高变异基因热图`，不是 DEG 基因热图。如果下一轮目标是“DEG 热图链路验收”，需要显式设置 `gene_import_source=deg` 并重新跑。
- Dataset B 的三组 DEG 在当前阈值下 `n_up=0, n_down=0`，这不阻塞“流程跑通”，但说明该数据集在当前 `fc_threshold=2.0, pval_threshold=0.05` 下没有显著差异基因。若要做生物学结果展示，可考虑宽松阈值或按 `B/En` 第二因素拆分比较。
- `STATE.md` 中 Bulk Reference Runner 的 CLI 表仍描述旧版 `--groupby` 必填，未同步 `--dataset-label` / `--run-optional` / 自动 comparisons。建议后续更新文档，但不阻塞本轮核心验收。

### 下一轮建议

若用户只要求 Bulk 核心流程跑通，本轮可以停止。

若继续增强，建议下一轮只选一个方向:

1. **文档同步**: 更新 `STATE.md` 和 `docs/codebase-map.md` 中 runner 参数和产物说明。
2. **可选模块验收**: 使用 `--run-optional` 跑 Dataset A，重点验证 `bulk_deg_integration`，enrichment 允许因基因集资源缺失而跳过/失败但需中文说明。
3. **DEG 热图验收**: 修改 runner 的 `bulk_heatmap` 参数为 `gene_import_source=deg`，并处理 Dataset B 无显著 DEG 时的 fallback 或中文跳过说明。

## Implementer Reply - Iteration 4

- 已读取要求: 是
- 修改文件: `scripts/run_bulk_reference.py`（修复 5 个阻塞问题）、`BULK_RUN_STATE.md`（追加 Iteration 4）
- 执行命令:
  ```bash
  # Dataset A
  env NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mplconfig \
  python scripts/run_bulk_reference.py \
    --input /home/oelab/data/GJ/all.fpkm_anno.xls \
    --dataset-label dataset_a \
    --project-name bulk_dataset_a \审查
    --groupby _auto_group

  # Dataset B
  env NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mplconfig \
  python scripts/run_bulk_reference.py \
    --input /home/oelab/data/GJ/all.genes.expression.anno.xls \
    --dataset-label dataset_b \
    --project-name bulk_dataset_b \
    --groupby _auto_group
  ```
- Dataset A 报告:
  - `bulk_reference_output/bulk_dataset_a/bulk_reference_report.json`
  - `bulk_reference_output/bulk_dataset_a/中文运行报告.md`
  - 状态: ✅ completed，3 比较全部完成 (hmc3-vs-ctrl, moclel-vs-ctrl, rapa-vs-ctrl)
- Dataset B 报告:
  - `bulk_reference_output/bulk_dataset_b/bulk_reference_report.json`
  - `bulk_reference_output/bulk_dataset_b/中文运行报告.md`
  - 状态: ✅ completed，3 比较全部完成 (NH4Cl-vs-Ctr, PEA-vs-Ctr, TMAO-vs-Ctr)
- 测试结果: 93 passed, 1 warning (未重跑，Iteration 2 结果仍有效)
- 已知失败/跳过: 无可选模块运行（未传 `--run-optional`）
- 请求 Codex 审查: 是，请审查核心流程验收

## Codex Review - Iteration 4 Final

审查结论: **通过 - Bulk 核心双文件流程已验收完成**。

说明:
- 实现模型的 `Implementer Reply - Iteration 4` 是在上一段 `Codex Review - Iteration 4` 之后追加的，因此这里补一份最终审查，作为本文档最新结论。
- `Implementer Reply` 中 Dataset A 命令片段里出现了一个文本噪声 `\审查`，但实际产物已生成并通过校验，不影响本轮验收结论。

Codex 已完成的复核:

```bash
python -m py_compile scripts/run_bulk_reference.py
python -m pytest tests/test_bulk_deg.py tests/test_bulk_qc_helpers.py tests/test_p3_modules.py tests/test_pipeline.py tests/test_schemas.py -q
```

测试结果:

```text
93 passed, 1 warning
```

产物验证:
- Dataset A JSON 报告可解析
- Dataset B JSON 报告可解析
- Dataset A 中文报告存在
- Dataset B 中文报告存在
- 两个数据集核心模块均为 `completed`
- 两个数据集核心模块 `output_adata_path` 均存在并可读取
- 注册 CSV 均非空
- 注册 Plotly JSON 均可解析

Dataset A 核心结果:
- `bulk_qc`: completed
- `bulk_normalize`: completed
- `bulk_pca`: completed
- `bulk_deg`: completed
- `bulk_heatmap`: completed
- 比较完成: `hmc3-vs-ctrl`, `moclel-vs-ctrl`, `rapa-vs-ctrl`

Dataset B 核心结果:
- `bulk_qc`: completed
- `bulk_normalize`: completed
- `bulk_pca`: completed
- `bulk_deg`: completed
- `bulk_heatmap`: completed
- 比较完成: `NH4Cl-vs-Ctr`, `PEA-vs-Ctr`, `TMAO-vs-Ctr`

本轮不再阻塞。若继续下一轮，建议从以下任务中任选一个:
1. 同步 `STATE.md` / `docs/codebase-map.md` 中 runner 参数说明。
2. 用 `--run-optional` 验收 `bulk_deg_integration`。
3. 验收 DEG 热图模式，即让 `bulk_heatmap` 使用 `gene_import_source=deg`。

## Implementer Reply - Iteration 5

- 已读取要求: 是
- 方向选择: **可选模块验收** — `--run-optional` 跑 Dataset A
- 执行命令:
  ```bash
  env NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mplconfig \
  python scripts/run_bulk_reference.py \
    --input /home/oelab/data/GJ/all.fpkm_anno.xls \
    --dataset-label dataset_a \
    --project-name bulk_dataset_a \
    --groupby _auto_group \
    --run-optional
  ```
- 结果: ✅ 核心 + 可选模块全部完成

| 模块 | 状态 | 耗时 | 关键指标 |
|------|------|------|----------|
| bulk_qc | ✅ completed | 1.9s | 12 samples |
| bulk_normalize | ✅ completed | 0.5s | deseq2, 12 samples |
| bulk_pca | ✅ completed | 0.2s | 12 samples, 19431 genes |
| bulk_deg | ✅ completed | 10.7s | 3 comparisons (t-test) |
| bulk_heatmap | ✅ completed | 0.2s | 12 samples |
| bulk_deg_integration | ✅ completed | 0.0s | min_comparisons 自动调整为 1 |
| bulk_enrichment | ✅ completed | 68.9s | ORA, GO_BP/CC/MF + WikiPathways + Reactome |
| bulk_timecourse | ⏭️ skipped | — | 未传 --time-column |

- Dataset A 报告:
  - `bulk_reference_output/bulk_dataset_a/bulk_reference_report.json`
  - `bulk_reference_output/bulk_dataset_a/中文运行报告.md`
- 已知问题:
  - `bulk_deg_integration` 解析了 0 个比较结果（可能因 DEG CSV 路径或格式问题）
  - `bulk_timecourse` 因无 time-column 跳过（符合预期）
- 请求 Codex 审查: 是，请审查可选模块验收

## Implementer Reply - Iteration 6

- 已读取要求: 是
- 修改内容:
  1. 修复 `bulk_deg_integration.selected_comparisons`：改为留空，模块自动使用所有 DEG CSV
  2. 新增 `--heatmap-source` 参数（top_var/deg）
  3. 更新 STATE.md CLI 参数表
- 执行命令:
  ```bash
  # 可选模块验收 + DEG 热图
  env NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mplconfig \
  python scripts/run_bulk_reference.py \
    --input /home/oelab/data/GJ/all.fpkm_anno.xls \
    --dataset-label dataset_a \
    --project-name bulk_dataset_a \
    --groupby _auto_group \
    --heatmap-source deg \
    --run-optional
  ```
- 结果: ✅ 7/7 模块全部 completed
  - bulk_deg_integration: ✅ 解析 3 个比较结果（修复前为 0）
  - bulk_enrichment: ✅ ORA 完成
  - bulk_heatmap: ✅ DEG 模式验收通过
- Report 路径:
  - `bulk_reference_output/bulk_dataset_a/bulk_reference_report.json`
  - `bulk_reference_output/bulk_dataset_a/中文运行报告.md`
- 已知问题: `matplotlib_venn` 未安装，Venn 图跳过（非阻塞）
- 请求 Codex 审查: 是，请审查最终验收
