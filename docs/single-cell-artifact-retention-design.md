# 任务产物（task_artifacts）清理策略设计

- 日期：2026-09-10
- 触发：审查发现单项目 `results/task_artifacts` 达 122 GB，无配额/清理/去重
- 目标：在不破坏"历史任务结果不可变、可追溯"这一既有保证的前提下，把增量与存量都压到可接受水平

## 1. 现状测量（只读实测）

| 项目 | task_artifacts 占用 |
| --- | --- |
| `c5bf72d7-220` | **122 GB** |
| `c75e7fc1-21a` | 40 GB |
| `ac674b57-61d` | 130 MB |
| 合计 | **约 162 GB** |

`c5bf72d7-220` 内部构成（`3380` 个文件，其中只有 `59` 个是 h5ad）：

| 模块 | 目录数 | h5ad 数 | h5ad 合计 | 说明 |
| --- | --- | --- | --- | --- |
| sc_cell_go | 21 | 21 | 44.3 GB | **直通**：模块不读 h5ad，却把上游文件整份复制 21 次 |
| functional_state | 16 | 16 | 33.6 GB | 每次重跑一份新产物 |
| annotation | 5 | 5 | 10.5 GB | 重跑历史 |
| sc_pseudobulk_deg | 4 | 4 | 8.4 GB | **直通**：返回 `output_adata=input_path` |
| qc | 4 | 4 | 6.1 GB | 重跑历史 |
| dimred / hvg | 2 / 2 | 4 | 8.3 GB | 重跑历史 |
| 其余（normalize/clustering/batch_correct/convert_10x） | 各 1 | 5 | 约 9.6 GB | 单次产物 |

关键事实：

1. **体积几乎全部来自 `output_adata.h5ad`**：59 个 h5ad ≈ 121 GB，其余 3321 个文件（CSV/图/manifest）合计不到 1 GB；
2. **全部 58 个产物目录都属于 `completed` 任务**，该项目没有分支（`analysis_branches` 为空），因此当前没有需要保护的引用；
3. **约 53 GB 是纯直通复制**：`sc_cell_go`(44.3) + `sc_pseudobulk_deg`(8.4) 的"输出"就是它们的输入；
4. 采样文件头指纹显示多个 2.255837128 GB 文件内容一致 → 同一上游被反复复制；
5. 全仓共有 **9 处** 模块把 `input_path` 直接作为 `output_adata` 返回：
   `sc_cell_go:1936`、`sc_pseudobulk_deg`/`sc_csv_export`（`sc_batch_export.py:630,1152`）、`sc_cell_deg:573`、`cell_communication.py:61,120,143`、`virtual_ko.py:381`、`convert_10x.py:91`、`batch_correct.py:979`、`bulk_deg_integration.py:79,133`；
6. 唯一写 `intermediate/*.h5ad` 的地方是 `base.py:508 save_output`，使用 `adata.write_h5ad()` **原地截断**（同一 inode）——这一点决定了硬链接方案必须先改成原子写（见 §4.2）。

## 2. 必须保持的不变量

1. **历史任务结果仍可读**：结果页的表格、图、manifest 必须继续存在（它们是小文件，不在清理范围内）；
2. **产物不可变**：重跑同一模块不得改写旧任务已登记的产物字节；
3. **分支保护**：`analysis_branches` 中 `accepted=1` 或未被软删除（`deleted=0`）的 `output_adata_path` 永不清理；
4. **运行安全**：任一任务/流水线运行期间不得删除该项目的任何产物；
5. **可审计**：每次删除都有记录（谁、何时、哪个文件、多大、为什么）；
6. **失败不影响分析**：清理逻辑异常绝不能导致任务失败。

## 3. 策略总览（三层，按收益/风险排序）

```
L1 源头减量   ← 收益最大、改动最小（约 20 行）
   1a 直通模块不复制（引用上游）
   1b 中间产物原子写 + h5ad 快照硬链接（0 字节）
   1c 单任务快照上限
L2 保留策略   ← 清理存量（GC 脚本 + 规则）
   每 (项目,模块) 保留最新 K 次成功产物；只删 output_adata.h5ad
L3 配额兜底   ← 防再次失控
   项目/磁盘水位阈值 + 系统状态页可见性
```

## 4. L1 源头减量

### 4.1 直通模块不复制（1a）

`worker._snapshot_task_artifacts` 已知 `result['output_adata']`，只需再传入 `input_path`：

```python
if os.path.realpath(output_path) == os.path.realpath(input_path):
    result['output_adata'] = output_path          # 记录引用，不复制
    result['artifact_mode'] = 'reference'          # 写入 manifest 便于审计
else:
    result['output_adata'] = copy_one(...)         # 原逻辑
```

- 依据：模块自己声明"我没有改动数据"，因此它的输出**就是**输入；再复制一份没有任何信息量。
- 代价：若上游是 `intermediate/x.h5ad`（可被后续重跑改写），旧任务登记的输出路径可能随之变化。补偿措施：manifest 记录 `artifact_mode='reference'` 与上游路径，且下游链式分析本来也读取同一文件。
- 收益（该项目）：`sc_cell_go` 44.3 GB + `sc_pseudobulk_deg` 8.4 GB → **约 53 GB 不再产生**。

### 4.2 硬链接快照 + 中间产物原子写（1b）

```python
# base.py save_output：原子写，重跑时换 inode，旧硬链接内容不被改写
tmp = f'{output_path}.tmp-{os.getpid()}'
adata.write_h5ad(tmp)
os.replace(tmp, output_path)
```

```python
# worker.copy_one：仅对平台生成的 intermediate 大文件用硬链接
def snapshot(source):
    if hardlink_enabled and source.startswith(intermediate_dir) and os.path.getsize(source) > 10 * 1024**2:
        try:
            os.link(source, target)   # 0 字节，且原子写保证不可变
            return target
        except OSError:
            pass                       # EXDEV/权限 → 回退复制
    shutil.copy2(source, target)
```

- 收益：每次"真实新产物"的快照额外占用从 2.2 GB 降到 0；删除 intermediate 时链接计数保护数据。
- 风险控制：跨文件系统自动回退；`uploads/`（上传原件，命名策略未逐一审计）与 `plots/`（`savefig` 原地覆盖）**仍走复制**；开关默认开，可用环境变量关闭。
- 附带收益：原子写本身修掉了"写入过程中崩溃留下半截 intermediate"的隐患。

### 4.3 单任务快照上限（1c）

快照前检查源文件大小：超过 `SC_ARTIFACT_MAX_TASK_GB`（默认 20 GB）时**跳过快照**，`output_adata_path` 直接记录中间产物路径，并在 manifest 写入 `artifact_snapshot_skipped: true` 与原因。宁可降低"不可变"强度，也不让单任务把磁盘打满。

## 5. L2 保留策略（GC）

### 5.1 索引文件（避免遍历）

每次快照时向 `<project>/results/task_artifacts/_index.jsonl` 追加一行：

```json
{"task_id":"...","module":"sc_cell_go","path":".../output_adata.h5ad","bytes":2255837128,
 "kind":"h5ad_snapshot","mode":"copy|hardlink|reference","created_at":"2026-09-10T18:16:27","status":"completed"}
```

GC 只读索引 + DB（`analysis_tasks`、`analysis_branches`），不做 `du`/`find` 全量遍历；索引缺失时才回退扫描。

### 5.2 规则

| 规则 | 内容 |
| --- | --- |
| **R0 范围** | **只删除文件名为 `output_adata.h5ad` 的文件**，且路径必须严格匹配 `<project>/results/task_artifacts/<module>/<task_id>/`（realpath + 段数校验）。CSV/图/manifest 一律不动，保证结果页完整 |
| **R1 保留最新** | 每个 `(project, module)` 保留最新 `K` 个成功任务的 `output_adata.h5ad`（`K = SC_ARTIFACT_KEEP_PER_MODULE`，默认 **1**，可设 2 以保留上一版对比） |
| **R2 失败/取消** | 对应任务 `status != 'completed'` 的产物**无条件**删除（UI 不提供它的链式输入） |
| **R3 分支保护** | 任何 `analysis_branches` 中 `deleted=0` 的 `output_adata_path` 及其所属目录永不删除；`accepted=1` 的记录额外在日志中标注 |
| **R4 运行保护** | 该项目存在运行中任务（`worker._active_futures` 或 `status='running'` 的 DB 记录）时，整个 GC 跳过 |
| **R5 引用保护** | 保留最新 K 个之外的候选，若其路径仍是某个**更新**任务 `result_json`/manifest 中记录的输入（`artifact_mode='reference'` 的引用源），则一并保留 |
| **R6 审计** | 每次实际删除追加 `_gc_log.jsonl`：时间、路径、bytes、规则、任务 ID、操作者（cli/auto）；同时 `logger.warning` 一行汇总 |

### 5.3 触发时机

- **手动（默认）**：`python scripts/sc_artifacts_gc.py --project <pid> [--all] [--keep N] [--apply]`，**不加 `--apply` 即 dry-run**，打印将删除的清单与可回收字节；
- **任务后自动（可选）**：`SC_ARTIFACTS_AUTO_GC=1` 时，任务完成回调里做一次轻量检查（项目产物 > `SC_ARTIFACT_MAX_PROJECT_GB`，默认 60 GB 才真正执行），失败只告警；
- **只报告不删除（默认）**：`SC_ARTIFACT_AUTO_GC` 默认 `0`，避免内部平台出现"结果悄悄消失"的意外。

### 5.4 删除后 UI 的行为（已核实）

`routes/analysis.build_input_options` 会 `os.path.isfile` 过滤，`_path_has_module_requirements` 对不存在的文件返回 False，因此被清理的旧 h5ad 只会从下拉框消失，旧任务的表格/图/manifest 仍正常展示，**不需要额外兼容层**。

## 6. L3 配额与可见性

- 快照前检查：项目产物合计 > `SC_ARTIFACT_MAX_PROJECT_GB`（默认 60）先触发一次 GC；仍超限则跳过快照并告警；
- 全局水位：`DATA_DIR` 所在分区可用空间 < `SC_ARTIFACT_DISK_FREE_PCT`（默认 15%）时，暂停所有新快照并写 `progress(-1, ...)` 告警；
- `/api/system/status` 增加只读字段：`artifacts_bytes`、`artifacts_reclaimable_bytes`、`artifacts_index_age`（供运维判断是否需要跑 GC）。

### 关于内容级去重（暂不实施）

内容寻址（`_blobs/<sha256>`）+ 引用计数能进一步合并"同一内容多模块"的副本，但：
- 本项目的重复**主要来自直通复制**，L1 已消除；
- sha256 一个 2.2 GB 文件需数秒 I/O，且 GC 需要维护引用计数，复杂度明显上升；
- 触发条件：当 L1+L2 落地后仍观察到"同一内容被 3 个以上模块各自快照"时再评估（符合 AGENTS.md"出现明确扩展触发条件才引入复杂架构"）。

## 7. 立即回收估算（`c5bf72d7-220`）

| 动作 | 释放 |
| --- | --- |
| R2：失败任务产物 | 0 GB（实测失败任务未留下目录） |
| R1：每模块保留最新 1 个 | 约 **99 GB**（121 GB → 约 22 GB 的"每模块最新一次"） |
| L1a：直通模块不再快照 | 长期增量再减约 **53 GB/轮** |
| L1b：硬链接 | 其余模块每次快照的 2.2 GB → **约 0** |

即：**只跑一次 GC 即可回收约 99 GB；L1 落地后同类项目的增量近乎为零**（仅保留每个模块最新一次的中间产物，约 20 GB/项目）。

## 8. 实施阶段

| 阶段 | 改动 | 风险 | 验证 |
| --- | --- | --- | --- |
| **P1 源头减量** | `worker._snapshot_task_artifacts` 加 `input_path` 参数与直通判断；`base.save_output` 原子写；`copy_one` 中间产物硬链接（带回退与开关） | 低：不改分析结果，只改文件如何落盘 | `test_pipeline.py` 增加"直通模块不产生副本"与"硬链接快照内容不变"断言；跑一次完整流水线对比 `du` |
| **P2 GC** | `scripts/sc_artifacts_gc.py` + `_index.jsonl` 写入 + 审计日志 | 中：涉及删除。默认 dry-run、只删 `output_adata.h5ad`、有 R0–R6 约束 | 先在 `c5bf72d7-220` dry-run 人工核对清单，再 `--apply`；核对该项目历史任务页表格/图仍可打开 |
| **P3 配额与可见性** | 快照前配额检查、水位暂停、`/api/system/status` 字段 | 低 | `test_system_health.py` 扩展；构造超限项目验证跳过与告警 |

建议先做 P1（约 20 行 + 2 个测试），把增量压到 0，再考虑 P2 清理存量——顺序反过来的话，清理完仍会继续长回来。

## 9. 配置项

| 环境变量 | 默认 | 含义 |
| --- | --- | --- |
| `SC_ARTIFACT_HARDLINK` | `1` | 中间产物快照使用硬链接（跨文件系统自动回退复制） |
| `SC_ARTIFACT_KEEP_PER_MODULE` | `1` | 每个 (项目,模块) 保留的最新成功产物个数 |
| `SC_ARTIFACTS_AUTO_GC` | `0` | `0`=只报告；`1`=超阈值自动清理 |
| `SC_ARTIFACT_MAX_TASK_GB` | `20` | 单任务快照大小上限，超过则跳过快照 |
| `SC_ARTIFACT_MAX_PROJECT_GB` | `60` | 项目产物软上限，触发 GC |
| `SC_ARTIFACT_DISK_FREE_PCT` | `15` | 分区可用空间低于该比例时暂停新快照 |

## 10. 未纳入本次设计

- `intermediate/`（该项目 23 GB）的清理：它是流水线工作区，删除会影响直通引用与 UI 的中间结果入口，建议 P1 落地后单独评估"只保留最近一次 `{module}_output.h5ad`"；
- 上传原件（`uploads/` 4.2 GB）与图表（`plots/` 265 MB）：体积小且是用户资产，不纳入自动清理。
