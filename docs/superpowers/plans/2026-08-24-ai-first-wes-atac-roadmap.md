# AI 主工作台、WES 与 ATAC 扩展项目规划

日期：2026-08-24  
状态：待确认  
建议周期：16–20 周（不含 scATAC，可并行配置 3–4 人）  
建议定位：从“模块选择型转录组网站”升级为“AI 驱动的多组学分析工作台”

使用场景约束：本平台首先服务实验室内部科研人员，由小型团队维护并运行少量并发任务。路线图中的企业级调度、完整事件溯源、多租户隔离和大规模集群能力均为按需升级项，不作为内部 MVP 的默认前置条件；实现复杂度遵循仓库根目录 `AGENTS.md`。

## 1. 结论与范围

本轮不重写现有单细胞和 Bulk RNA-seq 能力，而是在其上增加一个通用工作流内核，并把项目首页改造成以 AI 对话为主体的全屏工作台。

建议采用以下范围：

1. AI 对话从项目详情页右下角的浮层，升级为项目默认入口和主要操作界面。
2. 保留现有单细胞、Bulk RNA-seq 参数页面，作为“专家模式”和 AI 推荐参数的人工复核入口。
3. WES v1 同时支持：
   - germline 单样本/小队列 SNV、Indel；
   - 肿瘤-正常配对 somatic SNV、Indel；
   - tumor-only 作为受限模式，明确降低证据等级并要求 PoN/群体频率资源。
4. ATAC v1 先支持 bulk ATAC-seq：FASTQ 到 QC、比对、peak、consensus peak、差异可及性、motif 和报告。
5. scATAC/10x Multiome 放入 v2；它与 bulk ATAC 的数据模型、QC 和下游算法不同，不应在首版混为一个流程。
6. WES 首版定位为科研分析，不输出临床诊断或自动分级结论。

## 2. 当前基础与关键缺口

### 2.1 可以直接复用的能力

当前平台已经具备：

- Flask 项目、任务和结果路由；
- `AnalysisTask`、`ResultFile`、`PipelineRun`；
- 后台任务、进度、日志和失败状态；
- 模块注册表和参数 schema；
- AI 的 OpenAI-compatible/Anthropic-compatible 接入；
- 只读工具自动执行、写操作确认；
- Agent session、goal、step、branch 和候选评分；
- 任务 manifest、项目报告、结果图库和科研图导出；
- 项目目录边界、路径校验和符号链接拒绝。

这些能力应继续使用，不另建一套相互独立的 WES/ATAC 网站。

### 2.2 必须先消除的四个约束

| 现状 | 影响 | 改造方向 |
| --- | --- | --- |
| AI 对话是 `project_detail.html` 中 460×680 px 的浮层 | 无法承载样本设计、执行计划、实时日志和多结果证据 | 新建全屏 AI workspace；原浮层只保留为快捷入口 |
| 聊天历史是进程内 LRU，服务重启后丢失 | 无法形成可审计的长期科研会话 | 将 session、message、tool proposal/execution 持久化到数据库 |
| 流程类型硬编码为 `sc` / `bulk` | WES、ATAC 无法合法注册和校验 | 建立 assay/workflow registry，支持 `scrna`、`bulk_rna`、`wes`、`bulk_atac`、后续 `scatac` |
| Pipeline 每一步必须返回 `output_adata` | BAM、CRAM、VCF、BED、bigWig 无法串联 | 改为通用 `primary_output + artifacts + metrics + provenance` 结果契约 |

另一个运行层缺口是：当前 `ThreadPoolExecutor` 适合 Python 分析任务，但不适合直接承担多样本 WES/ATAC 重计算。后台线程应只负责启动、监控和收集外部工作流，重计算交给 Nextflow 的 local/Slurm 执行器与容器。

## 3. 产品目标与非目标

### 3.1 产品目标

- 用户进入项目后首先看到 AI 工作台，可以用自然语言完成“上传识别 → 设计检查 → 参数建议 → 确认 → 执行 → 结果解释 → 追问”。
- 每个 AI 结论都能回到具体任务、参数、文件、QC 指标和软件版本。
- AI 不直接拼 shell 命令；它只能调用后端注册、校验和授权过的工具。
- 单细胞、Bulk RNA、WES、Bulk ATAC 共用项目、资产、任务、报告和审计模型。
- 支持 OpenAI-compatible、Anthropic-compatible 以及提供 OpenAI-compatible API 的本地 vLLM/Ollama 类服务。
- 大文件和人类个体变异默认不发给外部模型；模型只接收结构化摘要。
- 流程可取消、失败可诊断、可从成功步骤断点续跑。

### 3.2 首版非目标

- 临床诊断、ACMG 自动定级、用药决策或医院 LIS 对接；
- 任意自然语言生成并执行 shell/Python/R 代码；
- 云厂商全托管部署；
- 大队列 germline 联合分型服务；
- scATAC、10x Multiome 和 RNA–ATAC 联合建模；
- WGS、甲基化、蛋白组学。

## 4. AI 主工作台设计

### 4.1 页面信息架构

新增默认路由：`/projects/<pid>/workspace`。

桌面端布局：

```text
┌──────────────────────────────────────────────────────────────────────────────┐
│ 项目 / assay / 参考版本 / 当前模型 / 运行环境 / 安全状态                    │
├───────────────┬─────────────────────────────────────┬────────────────────────┤
│ 会话与数据资产 │ AI 对话与执行时间线（主体，约 60%）  │ 项目上下文/任务/结果     │
│ 220–260 px     │                                     │ 320–400 px，可折叠       │
│               │ 目标、建议、审批、进度、证据卡       │ 样本表、QC、文件、日志   │
├───────────────┴─────────────────────────────────────┴────────────────────────┤
│ 多行输入框 / 附件 / 快捷动作 / 停止生成 / 发送                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

中心区必须是视觉和交互主体。右侧详情用于查看当前消息引用的任务或文件，不能抢占对话区。小屏时左右栏折叠为抽屉。

### 4.2 对话中的结构化卡片

不能只显示纯文本气泡。至少提供以下卡片：

1. **数据识别卡**：样本数、assay、文件配对、参考基因组、缺失字段。
2. **实验设计检查卡**：比较组、重复、肿瘤-正常配对、capture BED、阻断项和警告。
3. **分析计划卡**：阶段、方法、关键参数、预计资源、预计输出。
4. **执行审批卡**：将运行的 workflow、版本、输入、参数差异、外部数据策略和确认按钮。
5. **实时任务卡**：当前阶段、进度、已用时间、资源、日志摘要、取消和重试。
6. **结果证据卡**：核心 QC、图表缩略图、结论、限制、引用文件。
7. **候选比较卡**：沿用现有 branch 评分和采纳机制。
8. **下载交付卡**：报告、manifest、参数、软件版本和结果包。

### 4.3 会话能力

- 一个项目可有多个命名会话，例如“WES 肿瘤配对分析”“ATAC 差异可及性复核”。
- 会话可切换模型，但每条消息记录 provider、model、配置版本和时间。
- 支持服务端流式输出，前端采用 SSE；断线重连不能重复提交工具。
- 支持停止模型生成，但停止生成不等同于取消已批准的工作流。
- 支持引用项目文件、任务、结果和图，但默认不把原始大文件传给模型。
- 会话摘要用于控制上下文长度；不保存或展示模型隐式推理链。

### 4.4 “任意大模型”的边界

新增模型配置 `ModelProfile`：

| 字段 | 作用 |
| --- | --- |
| `provider_type` | `openai_compatible` / `anthropic_messages` |
| `base_url`、`model` | API 地址与模型名；密钥仍走服务端 secret |
| `supports_tools` | 是否支持原生 tool call |
| `supports_streaming` | 是否支持流式返回 |
| `supports_json_schema` | 是否能可靠输出结构化计划 |
| `context_window` | 上下文预算 |
| `data_policy` | `external_summary_only` / `local_sensitive_ok` |
| `enabled_tools` | 此 profile 可使用的工具白名单 |

不宣称所有模型都能完成同样的 Agent 工作。没有原生 tool call 的模型只能做解释和建议，不能执行分析；界面要明确显示能力徽标。

## 5. 通用数据与工作流内核

### 5.1 数据资产模型

新增 `data_assets`：

```text
id, project_id, sample_id, assay_type, role,
artifact_kind, format, file_path, size_bytes, checksum,
reference_build, metadata_json, parent_asset_ids,
created_by_task_id, created_at
```

核心 `artifact_kind`：

```text
fastq, bam, cram, bai, crai, vcf, gvcf, tbi,
bed, narrow_peak, broad_peak, bigwig,
count_matrix, h5ad, h5mu, table, figure, report, manifest
```

新增 `sample_manifests`，保存 assay、schema 版本、结构化内容、原始 CSV 路径、校验结果和版本。运行任务只能引用一个已校验的 manifest 版本。

### 5.2 通用结果契约

用以下契约替代 Pipeline 对 `output_adata` 的强依赖，同时保留兼容字段：

```json
{
  "status": "completed",
  "primary_output": {
    "path": "...",
    "kind": "vcf",
    "format": "vcf.gz"
  },
  "artifacts": [],
  "summary": {},
  "metrics": {},
  "provenance": {
    "workflow": "nf-core/sarek",
    "workflow_version": "pinned-release",
    "container_digests": [],
    "reference_asset_ids": []
  }
}
```

兼容策略：现有模块返回 `output_adata` 时，适配器自动转换为 `primary_output.kind=h5ad`；旧任务和现有测试不需要一次性迁移。

### 5.3 工作流注册表

新增声明式 `WORKFLOW_REGISTRY`，每个工作流包含：

- assay 类型；
- 输入 manifest schema；
- 可接受的输入资产类型；
- 参数 schema 与条件字段；
- preflight 校验器；
- executor 类型；
- 工作流版本；
- 结果收集器；
- AI 可见的能力、前置条件与风险文案。

现有 `MODULE_REGISTRY` 通过 `PythonModuleWorkflowAdapter` 接入，不删除。

### 5.4 执行器

定义统一接口：

```text
prepare(run) -> LaunchSpec
launch(spec) -> external_run_id
poll(external_run_id) -> status/progress
cancel(external_run_id)
resume(run_id)
collect(run_id) -> AnalysisResult
```

首版实现：

- `PythonModuleExecutor`：复用现有单细胞/Bulk RNA 任务；
- `NextflowExecutor`：运行 WES/Bulk ATAC；
- profile：`local + Docker/Apptainer`，预留 `Slurm + Apptainer`。

安全要求：

- 后端使用参数数组启动固定命令，禁止 `shell=True`；
- pipeline 名、版本、参数和挂载目录由白名单生成；
- 每次运行有独立 run/work/results 目录；
- 记录 PID/Nextflow run name，支持进程组取消；
- 使用幂等键避免用户重复点击造成双任务；
- Nextflow 版本和 nf-core release 必须固定，不能在生产运行时追踪 `dev`；
- `-resume` 只复用同一工作流版本、输入 checksum 和参数签名匹配的缓存。

Nextflow 原生提供中间结果跟踪、断点续跑以及 local/Slurm/云等执行层抽象，适合作为重计算后端，而不是把多进程生信命令塞入 Flask 线程。

### 5.5 参考资源中心

新增 `reference_assets`：

```text
id, species, assembly, asset_type, version, file_path,
checksum, source, metadata_json, status, created_at
```

需要管理：

- FASTA、FAI、sequence dictionary、BWA/BWA-MEM2 index；
- known sites/dbSNP、群体频率资源、PoN；
- WES capture target BED 和 padding BED；
- VEP cache 与对应 FASTA；
- GTF/TSS、chromosome sizes、blacklist；
- ATAC genome index。

所有运行前检查 assembly 一致性。GRCh37/GRCh38、chr/no-chr 命名和 capture BED 不一致必须阻断，不能只显示警告。

### 5.6 目录规划

```text
data/projects/<project_id>/
  uploads/                  # 小文件上传与现有兼容入口
  manifests/                # 版本化样本表
  workflow_runs/<run_id>/
    launch/                 # 固定后的配置与 samplesheet
    logs/                   # stdout/stderr/trace/timeline
    work/                   # Nextflow work，可按保留策略清理
    results/                # 原始工作流输出
  artifacts/                # 平台登记后的交付工件
  intermediate/             # 现有 h5ad 兼容
  plots/
  results/
  branches/
```

大 FASTQ 不建议通过浏览器复制上传。参考现有 `SC_BATCH_SOURCE_ROOTS`，新增管理员配置的 `GENOMICS_SOURCE_ROOTS`；只登记经过 allowlist 和 canonical path 校验的只读文件，所有输出仍写入项目目录。

## 6. WES v1 规划

### 6.1 输入模式

支持：

- paired-end FASTQ；
- 已排序并带 index 的 BAM/CRAM；
- 已压缩并带 index 的 VCF，从 annotation 阶段进入。

manifest 最少字段：

```text
patient_id, sample_id, role, sex, lane,
fastq_1, fastq_2, bam, bai, cram, crai,
reference_build, capture_bed_id
```

`role` 为 `normal`、`tumor` 或 `germline`。同一 somatic patient 必须能确定唯一匹配关系；样本 ID、read group、BAM header 不一致时阻断。

### 6.2 预检查

- FASTQ R1/R2、lane 和样本配对完整；
- 文件可读、checksum 可计算、压缩格式可验证；
- BAM/CRAM 排序、index、header 和参考版本正确；
- tumor-normal 配对完整，tumor-only 明确确认；
- capture BED、FASTA、known sites、VEP cache 版本匹配；
- 样本重名、角色冲突、混合参考版本阻断；
- 运行空间和预计 CPU/内存/临时存储满足要求。

### 6.3 主流程

推荐以固定版本 `nf-core/sarek` 为执行后端。Sarek 面向 WES/WGS/targeted sequencing，支持 germline、single-tumor、tumor-normal 和从 FASTQ/BAM/CRAM/VCF 不同阶段进入。

| 平台阶段 | 主要内容 | 默认实现 | 核心输出 |
| --- | --- | --- | --- |
| `wes_ingest` | manifest 与参考校验 | 平台 preflight | 固化 samplesheet、检查报告 |
| `wes_preprocess` | FASTQ QC、比对、排序、重复标记、BQSR、覆盖统计 | Sarek/GATK 工具链 | CRAM/BAM、index、MultiQC、coverage |
| `wes_germline_call` | germline SNV/Indel | HaplotypeCaller GVCF 工作流 | gVCF、filtered VCF |
| `wes_somatic_call` | tumor-normal/tumor-only SNV/Indel | Mutect2 + FilterMutectCalls | filtered somatic VCF、stats |
| `wes_annotate` | consequence annotation | VEP 本地 cache/offline | annotated VCF、可查询表 |
| `wes_report` | 样本、覆盖、变异和限制摘要 | 平台 report adapter | HTML/Markdown/JSON/CSV |

Germline 与 somatic 不共用 caller：HaplotypeCaller 用于 germline SNP/Indel，Mutect2 用于 somatic SNV/Indel。VEP 默认本地 cache/offline，避免把个体变异坐标发送到公共数据库。

### 6.4 结果界面

WES 结果页至少包含：

- 样本身份与配对状态；
- reads、比对率、重复率、污染估计、平均 target coverage、≥20×/≥30× 覆盖比例、on-target rate；
- germline/somatic PASS 数量和 Ti/Tv 等集合级指标；
- 可筛选变异表：gene、consequence、HGVS、depth、VAF、FILTER、population AF；
- 每个变异回链到原始 VCF 记录和 annotation 版本；
- 下载 CRAM/BAM、VCF、index、QC、manifest、完整 provenance；
- “科研用途，不构成临床诊断”提示。

### 6.5 WES v1.1（不阻塞首版）

- CNV：CNVkit/同类方法；
- MSI：需专门工具和适用输入检查；
- TMB：必须绑定 capture territory、PASS 过滤规则和版本，不给出脱离 panel/外显子区域定义的数值；
- cohort joint genotyping；
- 自定义 annotation 插件和已审核知识库；
- IGV.js 局部证据查看，但不把 BAM 发送给 LLM。

## 7. Bulk ATAC-seq v1 规划

### 7.1 输入和实验设计

manifest 最少字段：

```text
sample_id, condition, replicate, batch,
fastq_1, fastq_2, control_sample_id,
reference_build
```

差异可及性要求每组至少 2 个独立生物学重复，建议 3 个；不满足时允许做 QC、peak 和描述性比较，但必须阻断正式推断或标为探索性结果。技术 lane 不能被误当作生物学重复。

### 7.2 主流程

推荐固定版本 `nf-core/atacseq` 负责原始处理和 QC，再由平台封装差异可及性与结果解释。

| 平台阶段 | 主要内容 | 核心输出 |
| --- | --- | --- |
| `atac_ingest` | manifest、重复、参考、blacklist 检查 | 固化 samplesheet、设计检查 |
| `atac_qc_align` | raw/trim QC、比对、过滤、去重复、线粒体读段统计 | BAM、MultiQC、insert size、复杂度 |
| `atac_signal_peaks` | Tn5 shift、signal track、MACS2 peak | bigWig、narrowPeak/broadPeak、summits |
| `atac_consensus` | replicate/condition consensus peak 与计数 | consensus BED、peak count matrix |
| `atac_differential` | 设计矩阵、contrast、差异可及性 | 完整结果表、火山图、MA、热图、PCA |
| `atac_motif` | peak annotation、motif enrichment | motif 表、富集图、peak-to-gene 表 |
| `atac_report` | QC gate、统计限制和交付包 | HTML/Markdown/JSON/CSV |

### 7.3 QC 与结果卡

核心指标：

- raw/trim read quality；
- mapping rate、duplicate rate、mitochondrial fraction；
- insert-size/nucleosome periodicity；
- TSS enrichment；
- FRiP；
- peak 数、blacklist fraction、library complexity；
- replicate correlation、PCA 和可选 IDR/一致性证据。

QC gate 按参考注释和组织类型配置，不在代码中写一个对所有项目都适用的单阈值。界面可引用 ENCODE 的参考等级，但允许用户在项目里记录例外理由。

### 7.4 ATAC v2：scATAC / Multiome

单独立项，预计再增加 4–6 周：

- 10x fragments/peak matrix 导入；
- barcode 级 TSS、fragment、nucleosome、doublet QC；
- LSI/UMAP/聚类；
- peak calling、gene activity、motif deviation；
- RNA–ATAC label transfer 和联合可视化；
- h5ad/h5mu 或 Arrow/RDS 资产策略。

进入 v2 前先做 SnapATAC2、ArchR、Signac 的技术选型 ADR，不能仅因现有平台是 Python 就忽略成熟 R 生态，也不能让不同对象格式在同一任务链中隐式转换。

## 8. AI 工具扩展

### 8.1 新增只读工具

- `list_assay_capabilities`
- `inspect_project_assets`
- `validate_sample_manifest`
- `inspect_reference_assets`
- `inspect_experimental_design`
- `get_workflow_run_status`
- `inspect_wes_qc`
- `inspect_atac_qc`
- `query_variant_summary`
- `query_peak_summary`
- `compare_workflow_runs`

### 8.2 需要确认的写工具

- `create_manifest_from_assets`
- `run_workflow`
- `resume_workflow`
- `cancel_workflow`
- `run_downstream_analysis`
- `accept_workflow_baseline`

审批卡必须显示：

- workflow 与固定版本；
- 样本数和输入角色；
- 参考资源版本；
- 关键参数与默认值差异；
- 资源和磁盘估算；
- 是否使用外部模型/服务；
- 预计产生的主要工件；
- 阻断项、警告和可撤销性。

取消运行属于有副作用操作，也需要二次确认；查看状态、QC 和结果摘要仍自动执行。

### 8.3 模型上下文策略

模型可见：

- manifest 摘要；
- 结构化 QC 和统计结果；
- 工件元数据、文件名和任务 lineage；
- 用户选择的小规模表格切片；
- 报告摘要和图表描述。

外部模型默认不可见：

- FASTQ/BAM/CRAM 原始内容；
- 全量 VCF 和个体变异坐标；
- 直接身份信息；
- 服务器绝对路径和 secret；
- 未经脱敏的临床字段。

## 9. 数据库与 API 改造

### 9.1 新表

- `chat_sessions`
- `chat_messages`
- `tool_proposals`
- `tool_executions`
- `data_assets`
- `sample_manifests`
- `reference_assets`
- `workflow_runs`
- `workflow_artifacts`

现有 `analysis_tasks` 和 `result_files` 保留。`workflow_runs` 可通过 `task_id`/`pipeline_run_id` 连接旧任务，避免一次迁移全部历史数据。

### 9.2 关键 API

```text
GET/POST /api/projects/<pid>/chat/sessions
GET      /api/chat/sessions/<sid>/messages
POST     /api/chat/sessions/<sid>/messages
GET      /api/chat/sessions/<sid>/events

GET/POST /api/projects/<pid>/assets
GET/POST /api/projects/<pid>/manifests
POST     /api/projects/<pid>/manifests/<mid>/validate

GET/POST /api/projects/<pid>/workflow-runs
GET      /api/workflow-runs/<rid>
GET      /api/workflow-runs/<rid>/events
POST     /api/workflow-runs/<rid>/cancel
POST     /api/workflow-runs/<rid>/resume
GET      /api/workflow-runs/<rid>/artifacts
```

所有 mutation API 使用项目权限、CSRF/Token、幂等键和审计日志。当前平台若继续处理人类 WES 数据，正式部署前必须增加用户认证、项目级权限和访问审计，不能只依赖可选的 AI API Token。

## 10. 代码落点建议

```text
routes/
  workspace.py
  workflows.py
  assets.py

templates/
  workspace.html

static/
  css/workspace.css
  js/workspace/
    chat.js
    sessions.js
    tool_cards.js
    run_timeline.js
    artifact_panel.js

modules/
  llm/
    providers/base.py
    providers/openai_compatible.py
    providers/anthropic_messages.py
    context_builder.py
  workflows/
    contracts.py
    registry.py
    service.py
    executors/python_module.py
    executors/nextflow.py
    collectors/common.py
    wes/manifest.py
    wes/preflight.py
    wes/sarek.py
    wes/collector.py
    atac/manifest.py
    atac/preflight.py
    atac/nfcore_atacseq.py
    atac/collector.py
  assets/
    registry.py
    references.py
    checksums.py
```

应把 `project_detail.html` 中大段内联聊天 CSS/JS 抽出，不继续在同一个模板内叠加 WES/ATAC 功能。

## 11. 分阶段实施计划

### 阶段 0：范围冻结与基线（第 1 周）

交付：

- 确认 WES 三种模式、ATAC v1=bulk；
- 确认运行环境：单机 Docker/Apptainer 或 Slurm；
- 确认 GRCh38 首选及需要兼容的参考；
- 建立最小公开测试数据和现有回归基线；
- 完成通用结果契约、资产模型、敏感数据策略 ADR。

验收：范围、接口、数据字典和科学验证数据均有版本化文档。

### 阶段 1：AI 全屏工作台（第 2–4 周）

交付：

- workspace 路由和三栏响应式界面；
- 持久化会话、消息、tool proposal 和 execution；
- SSE 流式消息、停止生成和重连；
- 模型 profile 与能力徽标；
- 现有 AI 工具、确认卡和 Agent branch 面板迁入工作台；
- 原项目详情页的浮层改为“打开 AI 工作台”。

验收：服务重启后会话不丢；OpenAI-compatible 与 Anthropic-compatible 各通过一套契约测试；重复点击不会重复执行工具。

### 阶段 2：通用工作流与资产内核（第 4–6 周）

交付：

- data asset、manifest、reference asset、workflow run；
- 通用 AnalysisResult；
- Python module compatibility adapter；
- Nextflow executor POC、日志、取消和 resume；
- workflow/artifact 页面和 API；
- 系统健康检查增加 Nextflow、容器、参考资源、磁盘。

验收：同一 Pipeline 可串联非 h5ad 工件；运行中重启 Web 服务后仍能恢复监控；失败任务可以从缓存续跑。

### 阶段 3：WES MVP（第 7–11 周）

交付：

- WES manifest、上传/服务器资产登记、preflight；
- 固定版本 Sarek profile；
- germline、tumor-normal、tumor-only；
- VEP offline annotation；
- QC/variant 结果页、AI 只读解释工具、交付报告；
- 公开 truth/test 数据验证。

验收：三种输入入口至少各有一条 E2E；样本角色互换、reference/capture 不一致能被阻断；所有输出可回溯到 workflow/container/reference 版本。

### 阶段 4：Bulk ATAC MVP（第 12–15 周）

交付：

- ATAC manifest 和实验设计 preflight；
- 固定版本 nf-core/atacseq；
- QC、alignment、signal、peak、consensus peak；
- 差异可及性、motif/peak annotation；
- AI QC 解释和结果交付报告。

验收：有重复的公开数据完成 E2E；无生物学重复时统计推断被阻断或明确降级；bigWig、peak、count matrix、差异表和 QC 均登记为工件。

### 阶段 5：AI 闭环与跨模块体验（第 16–17 周）

交付：

- AI 能从用户目标选择 assay workflow；
- 自动读取 manifest/QC/任务状态，生成下一步建议；
- 失败原因解释、resume 提案、运行对比；
- 结果证据卡和文件引用；
- 项目级多组学结果索引，但不做未经定义的跨组学统计整合。

验收：预定义的 20 条自然语言场景均产生正确工具、参数和确认边界。

### 阶段 6：加固与试点（第 18–20 周）

交付：

- 用户认证、项目权限、审计、CSRF、下载授权；
- 磁盘配额、清理策略、任务超时、资源限制；
- 科学验证、性能压测、故障演练；
- 管理员部署手册、用户手册、示例项目；
- 1–2 个真实试点项目的反馈修订。

验收：安全、科学、E2E、恢复性验收全部通过后才标记生产可用。

## 12. 人员与工作量

建议最低配置：

| 角色 | 投入 | 主要职责 |
| --- | --- | --- |
| 后端/工作流工程师 | 1 人全程 | 通用内核、Nextflow、任务状态、资产/API |
| 前端/全栈工程师 | 1 人，前 8 周为主 | AI workspace、流式事件、卡片、结果浏览 |
| 生信工程师 | 1 人全程 | WES/ATAC profile、参考资源、科学验证、报告 |
| QA/DevOps | 0.5 人 | 容器、Slurm/单机、E2E、恢复、安全和发布 |
| 科研负责人 | 每周评审 | 分析范围、阈值、结果解释和试点验收 |

若只有 1 名全栈生信开发者，建议按 24–32 周估算，并严格保持 WES → Bulk ATAC → scATAC 的顺序。

## 13. 测试与验收体系

### 13.1 软件测试

- schema/manifest/preflight 单元测试；
- provider 与 tool call 契约测试；
- workflow collector 的 fixture 测试；
- API 权限、幂等、路径和参数注入测试；
- SSE 断线重连与重复事件测试；
- executor 的启动、取消、失败、resume 和服务重启恢复测试；
- 现有 scRNA/Bulk RNA 全量回归。

### 13.2 科学验证

WES：

- nf-core/Sarek 官方 test profile 做安装冒烟；
- germline 使用 NIST Genome in a Bottle 的高置信区 truth set，在固定 capture 区域用标准比较工具建立版本基线；
- somatic 使用公开配对/合成 spike-in 小数据验证 caller、filter、VAF 和样本配对；
- 每次工作流升级都与已批准版本比较，不用任意的单一 precision/recall 数值作为所有实验的硬门槛。

ATAC：

- nf-core/atacseq 官方 test 数据做安装冒烟；
- 使用有生物学重复的公开数据验证 TSS、FRiP、nucleosome pattern、replicate correlation、peak 和差异方向；
- QC gate 记录参考注释版本和阈值来源。

### 13.3 产品验收指标

- 90% 的标准流程可以从 AI 工作台完成，不需要切到命令行；
- 100% 写操作有明确确认和审计记录；
- 100% 结果文件有 task/workflow/reference lineage；
- 服务重启不丢聊天、运行状态和审批记录；
- 重复提交不会产生两个等价运行；
- 失败任务能显示阶段、错误摘要、日志位置和可行恢复动作；
- 外部模型调用日志证明默认未传输 FASTQ/BAM/CRAM/全量 VCF。

## 14. 主要风险与控制

| 风险 | 控制措施 |
| --- | --- |
| WES 范围过大 | v1 只做 SNV/Indel + annotation；CNV/MSI/TMB 分到 v1.1 |
| 把 bulk ATAC 与 scATAC 混做 | v1 明确 bulk；scATAC 单独数据模型和 ADR |
| 当前 h5ad 契约阻塞新工件 | 先交付通用 AnalysisResult/asset 层 |
| Flask worker 被重任务拖垮 | Web 线程只编排；Nextflow + container/HPC 执行 |
| 参考资源错配导致静默错误 | reference asset ID、checksum、assembly 和 contig preflight 硬阻断 |
| LLM 幻觉参数或解释 | schema 白名单、确定性 preflight、工具确认、证据卡和引用 |
| 人类基因组隐私 | 本地/离线处理、最小化模型上下文、RBAC、审计和下载授权 |
| pipeline/tool 自动更新破坏复现 | 固定 release、container digest、reference/annotation 版本 |
| 磁盘快速增长 | 配额、run/work/results 分层、可审计清理和归档策略 |

## 15. 开工前需要确认的四个产品决策

以下是本规划采用的默认选择；如无反对即可据此进入实施：

1. **ATAC 首版范围**：bulk ATAC-seq；scATAC/Multiome 为 v2。
2. **WES 首版范围**：germline + tumor-normal + 受限 tumor-only；CNV/MSI/TMB 为 v1.1。
3. **运行后端**：Nextflow；开发环境 local + Docker/Apptainer，生产预留 Slurm + Apptainer。
4. **参考策略**：优先 GRCh38，允许登记 GRCh37，但单次 workflow 绝不混用。

## 16. 第一批可执行任务（确认后前两周）

1. 新建三份 ADR：通用工件契约、Nextflow 执行、敏感基因组数据与 LLM 边界。
2. 为现有 `output_adata` 增加兼容适配测试，先不改现有模块内部实现。
3. 建立 `chat_sessions/chat_messages/tool_proposals/tool_executions` migration 和模型类。
4. 新建 `/projects/<pid>/workspace`，迁移现有聊天、确认卡和 Agent 面板。
5. 将内联 CSS/JS 拆分到 `static/css` 和 `static/js/workspace`。
6. 建立 provider contract 和两类现有 API adapter 的回归测试。
7. 建立 `data_assets/sample_manifests/reference_assets/workflow_runs` 数据字典与 migration 草案。
8. 用 Nextflow 官方 hello/test 工作流完成启动、日志、取消、resume POC。
9. 下载/固定 nf-core/sarek 与 nf-core/atacseq 的目标 release，记录 container/reference 需求。
10. 选择小型公开 WES、ATAC fixtures，建立离线 CI 与完整 E2E 两级测试策略。

第一批 Definition of Done：全屏 AI 工作台能完整复用现有工具；聊天持久化；通用结果契约测试通过；Nextflow POC 可启动、查看、取消并 resume；未开始实现未经确认的 WES/ATAC 业务参数。

## 17. 技术依据

- [nf-core/sarek](https://nf-co.re/sarek/latest/docs/usage) 面向 WGS/WES/targeted sequencing 的 germline/somatic 预处理、变异检测和注释，并支持从 FASTQ、BAM/CRAM、VCF 等阶段进入。
- [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360036365812-HaplotypeCaller) 用于 germline SNP/Indel，并可通过 GVCF 工作流进入联合分型。
- [GATK somatic short variant workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels) 和 [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/9570422171291-Mutect2) 支持 tumor-normal 与 tumor-only somatic SNV/Indel。
- [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html) 支持本地 cache/offline 注释，适合减少敏感变异坐标的外发。
- [nf-core/atacseq](https://nf-co.re/atacseq/2.1.2/docs/output/) 输出 alignment、MultiQC、peak、FRiP、consensus peak 和下游结果所需工件。
- [ENCODE ATAC-seq 数据标准](https://www.encodeproject.org/atac-seq/) 提供 FRiP、TSS enrichment 和 nucleosome pattern 等 QC 参考。
- [Nextflow](https://nextflow.io/) 提供持续 checkpoint、resume 和多执行平台抽象。
- [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) 提供人类小变异 benchmark 参考数据。
