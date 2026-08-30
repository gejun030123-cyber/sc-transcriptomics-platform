# WES（全外显子组测序）平台补充方案

**状态：** 评审修订稿 v0.2  
**日期：** 2026-08-24  
**定位：** 在现有单细胞 / Bulk RNA-seq 平台上补充科研级 WES 能力  
**使用场景：** 实验室内部科研人员、小型维护团队、少量并发任务  
**适用范围：** 人类短读长、双端测序的 germline SNV/InDel 与肿瘤 somatic SNV/InDel  
**上位设计：** `docs/superpowers/plans/2026-08-24-ai-first-wes-atac-roadmap.md`

**阶段命名约定：** 后续开发、验收和汇报只使用 P0–P3。旧稿中的 M0–M7
为早期拆分方式，已由第 13 节的四阶段计划取代，避免把 “M3” 与当前 P2/P3
混用。

---

## 0. 评审结论与本版决策

原 v0.1 的科学方向可保留，但不能直接按“给现有 `BaseAnalysis` 增加 10 个模块”的方式实施。WES 具有多样本配对、多输入工件、长任务、并行分支、参考资源强约束和人类遗传数据敏感等特点，必须先补充通用资产与外部工作流执行能力。

本版作出以下调整：

| 事项 | v0.2 决策 |
| --- | --- |
| 平台架构 | WES 接入 `WORKFLOW_REGISTRY`，不直接塞入现有线性 `PIPELINE_ORDER` |
| 实施复杂度 | 按实验室内部使用规模采用轻量增量实现；企业级调度、完整事件溯源和多租户能力按实际触发条件延期 |
| 执行后端 | 固定版本 nf-core/sarek + Nextflow；Web worker 只负责启动、监控、取消、恢复和收集 |
| v1 核心范围 | FASTQ/BAM/CRAM/VCF 分阶段入口；germline 与 tumor-normal/tumor-only SNV/InDel；离线注释；QC；报告；完整 provenance |
| 默认 germline caller | GATK HaplotypeCaller gVCF 工作流；DeepVariant 仅作为独立备选，不做并集/交集 |
| 默认 somatic caller | Mutect2 + 官方过滤链；不默认组合 Strelka2/VarScan2 |
| 首期移出 | CNV、SV、MSI、TMB、突变特征、ACMG 自动分级、大队列 joint genotyping |
| tumor-only | 支持受限模式，必须明确确认并降低证据等级；matched normal 为默认推荐 |
| 注释器 | VEP offline/cache 为默认；COSMIC 等许可资源按部署环境可选 |
| 参考基因组 | 新项目仅启用一个经过验证的 GRCh38 bundle；GRCh37/hg19 是独立 bundle 和独立验证，不是运行时随意切换 |
| 安全 | 真实人类 WES 试点前必须完成认证、项目权限、下载审计和数据保留策略 |

---

## 1. 背景、目标与边界

### 1.1 平台价值

现有平台主要分析表达层面的生物学状态。WES 可补充：

- germline SNV/InDel 与家系遗传模式候选；
- 肿瘤 somatic SNV/InDel、VAF 与候选驱动基因；
- 与现有单细胞 / Bulk RNA-seq 项目共享样本、任务、结果、报告和审计框架；
- 为后续表达—基因型联合解释保留标准化资产接口。

WES 不能被描述为对所有基因组变异的完整检测。短读长外显子捕获对非编码区、重复区、低覆盖外显子、复杂结构变异和多数非外显子断点存在天然盲区。

### 1.2 v1 产品目标

1. 建立可复现、可取消、可恢复、可审计的 WES 外部工作流执行能力；
2. 支持 FASTQ、已处理 BAM/CRAM、已调用 VCF 三类合法入口；
3. 支持 germline 单样本/小家系和 somatic tumor-normal/tumor-only；
4. 交付 filtered VCF、离线注释表、QC、日志、参数、参考资源和软件版本；
5. 每个结论能回链至样本 manifest、原始 VCF 记录、过滤状态和注释版本；
6. 对公开 truth set 完成分层 precision/recall 验证。

### 1.3 v1 非目标

- 临床诊断、用药决策或医院 LIS 对接；
- 自动输出平台认可的 ACMG/AMP 五级结论；
- 大队列人群分析、队列级关联分析或长期 joint-genotyping 服务；
- CNV、SV、MSI、TMB、突变特征的生产级交付；
- WGS、RNA variant calling、RNA/WES 自动联合建模；
- 将 FASTQ、BAM、CRAM 或全量 VCF 发送给外部大模型。

---

## 2. 设计原则

1. **样本设计先于运行：** 没有通过校验的版本化 manifest，不能启动工作流。
2. **资产而非路径：** 任务引用 `asset_id` 和 checksum；绝对路径只是服务端实现细节。
3. **入口分层：** FASTQ、BAM/CRAM、VCF 从各自合法阶段进入，不重复处理也不跳过必要校验。
4. **默认单一可信 caller：** 多 caller 只有在预定义合并规则并通过独立验证后才能成为正式模式。
5. **参考资源不可混用：** FASTA、索引、known sites、群体资源、capture BED 和 VEP cache 必须属于同一 assembly bundle。
6. **计算与 Web 解耦：** Web 线程不直接承担多小时重计算。
7. **科研结果有边界：** 自动注释、证据聚合和致病性分类是三个不同层次，不混称“诊断”。
8. **失败必须显式：** 必需工件缺失、登记失败、checksum 不一致或结果收集不完整均使任务失败。
9. **敏感数据最小暴露：** AI 默认只读取结构化摘要和经授权的小规模表格切片。
10. **内部使用优先：** 在满足科学正确性、数据安全和任务可恢复基本要求的前提下，优先复用现有 Flask/SQLite 架构和简单可维护实现，不为假设中的企业级规模提前建设复杂内核。

---

## 3. 平台前置改造

### 3.1 当前约束

当前平台存在以下 WES 阻断项：

- `analysis_type` 仅支持 `sc` / `bulk`；
- Pipeline 以单个 `input_path` 串联，并要求每一步返回 `output_adata`；
- `PIPELINE_DEPS` 不能表达 OR 依赖、并行分支和多工件汇合；
- `ThreadPoolExecutor` 的运行状态仅在进程内，Web 重启后不能可靠恢复监控；
- `ResultFile` 仅表达展示文件，不能表达 BAM+BAI、VCF+TBI 等成组工件及 lineage；
- 大文件上传、用户认证、项目级授权和下载审计尚不足以承载真实人类 WES。

### 3.2 必须新增的数据模型

#### `data_assets`

```text
id, project_id, sample_id, assay_type, role,
artifact_kind, format, file_path, size_bytes, checksum,
reference_build, metadata_json, parent_asset_ids,
created_by_run_id, created_at
```

首期 `artifact_kind`：

```text
fastq, bam, bai, cram, crai, vcf, gvcf, tbi,
bed, interval_list, table, figure, report, manifest, log
```

#### `sample_manifests`

保存 schema 版本、原始 CSV/JSON、结构化内容、校验报告、checksum、状态和版本。已经被任务使用的 manifest 不允许原地覆盖，只能创建新版本。

#### `reference_assets`

```text
id, species, assembly, bundle_version, asset_type,
file_path, checksum, source, license_note, metadata_json, status
```

#### `workflow_runs` 与 `workflow_artifacts`

记录 workflow release、executor、external run ID、输入资产、参数签名、状态、资源、时间、日志、工件和 provenance。现有 `AnalysisTask`/`ResultFile` 保留为兼容视图。

### 3.3 通用结果契约

```json
{
  "status": "completed",
  "primary_output": {
    "asset_id": "asset_xxx",
    "kind": "vcf",
    "format": "vcf.gz"
  },
  "artifacts": [],
  "metrics": {},
  "summary": {},
  "warnings": [],
  "provenance": {
    "workflow": "nf-core/sarek",
    "workflow_version": "3.10.0",
    "nextflow_version": "pinned-in-deployment-lock",
    "container_digests": [],
    "reference_bundle_id": "",
    "input_manifest_id": "",
    "parameter_signature": ""
  }
}
```

现有模块通过兼容适配器将 `output_adata` 转为 `primary_output.kind=h5ad`，不要求一次性迁移旧任务。

---

## 4. 工作流与执行架构

### 4.1 工作流注册表

新增 `WORKFLOW_REGISTRY`。每个工作流声明：

- assay 与支持的分析模式；
- manifest schema 和允许的输入工件组合；
- preflight 校验器；
- 参数 schema 与不可修改的安全参数；
- workflow release、profile 和容器 digest；
- executor、结果收集器和必需工件清单；
- AI 可见能力、风险说明和审批要求。

建议注册三个用户级工作流，而不是十个可任意拼接的小模块：

| workflow key | 用途 | 合法入口 |
| --- | --- | --- |
| `wes_germline` | 单样本/小家系 germline SNV/InDel | FASTQ、BAM/CRAM |
| `wes_somatic` | tumor-normal 或受限 tumor-only SNV/InDel | FASTQ、BAM/CRAM |
| `wes_annotate_only` | 已有 VCF 的标准化与离线注释 | VCF.gz + TBI |

网页仍可把内部阶段显示为 ingest、preprocess、call、filter、annotate、report，但用户不能拼出科学上非法的顺序。

### 4.2 执行器

```text
prepare(run) -> LaunchSpec
launch(spec) -> external_run_id
poll(external_run_id) -> status/progress/resources
cancel(external_run_id)
resume(run_id)
collect(run_id) -> AnalysisResult
```

首期采用：

- Nextflow + nf-core/sarek；
- local + Docker/Apptainer profile，预留 Slurm + Apptainer；
- workflow 固定为经验证 release，禁止使用 `dev` 或浮动 `latest`；
- 容器按 digest 固定；
- `-resume` 仅在 workflow版本、输入 checksum、参考 bundle 和参数签名完全一致时使用；
- 每次运行使用独立的 run/work/results/logs 目录；
- 记录 PID、进程组或 Nextflow run name，支持真正取消，而不是只修改数据库状态。

当前内部版已落地上述接口的轻量实现：`prepare` 会生成可审阅的 Sarek
samplesheet 和 launch bundle；`launch`、`GET run`、`cancel`、`resume` 和日志 tail
均由 Flask API 提供。`WES_EXECUTOR_ENABLED` 默认关闭，管理员确认 Nextflow、固定
profile、容器和 reference bundle 已部署后才开启；关闭时仍可完整审阅 bundle，但不会
启动外部进程。运行状态只保留当前实验室规模所需的 `prepared/running/completed/
failed/cancel_requested/cancelled/interrupted` 等状态，不提前建设调度中心或运行事件表。
Sarek 的 `step`、`tools` 和 `genome` 由 workflow key、入口类型与管理员配置
固定生成：germline 使用 HaplotypeCaller，somatic 使用 Mutect2，VCF 入口使用
`annotate + vep`；FASTQ 从 `mapping` 开始，BAM/CRAM 从 `variant_calling`
开始。网页参数不能覆盖这些 caller 选择。Nextflow 26 下布尔参数通过
服务端生成的 `parameters.json` 以 `-params-file` 传入，启动时再校验文件
与固定参数一致。germline/somatic calling 不得只登记 capture BED ID，还必须
提供真实、存在且非符号链接的 `capture_bed_path`。对 tumor-only
manifest，preflight 还必须同时收到非空 `tumor_only_reason` 和布尔值
`tumor_only_confirmed=true`；字符串 `"true"` 不视为有效确认。
默认 profile 是 Sarek 的 `docker` profile（仍由 Nextflow local executor 执行），
可由管理员切换到已验证的 Apptainer 或机构 profile。

Web worker 仅启动和监控外部运行。服务重启后根据 `external_run_id` 重新关联，不能把“存在 checkpoint 文件”视为已经实现恢复。

### 4.3 建议目录

```text
data/projects/{pid}/
├── manifests/
├── workflow_runs/{run_id}/
│   ├── launch/             # 固化 samplesheet、参数、profile、checksums
│   ├── logs/               # stdout/stderr/trace/timeline/report
│   ├── work/               # Nextflow work；按保留策略清理
│   └── results/            # 工作流原始输出
├── artifacts/              # 平台登记后的交付工件
├── results/                # 兼容现有结果页
└── plots/
```

---

### 4.4 Python 包与外部工具边界

WES 不在 Python 中重写变异检测器。平台 Python 层只负责输入契约、轻量预检查、任务状态、资产登记和结果索引；大文件计算交给经过验证的外部工具。

| 能力 | 平台侧参考 | 责任边界 |
| --- | --- | --- |
| BAM/CRAM/FASTA header 与索引检查 | `pysam`（可选依赖） | 读取 header、read group、参考序列和索引状态；不替代 `samtools` |
| VCF/BCF 记录读取 | `pysam.VariantFile`；`cyvcf2` 可作为高吞吐可选项 | 查询 FILTER、DP、AD、VAF 等；不替代 `bcftools` 标准化 |
| FASTQ 完整性 | Python `gzip`、`seqkit`/`fastp` | Python 只做快速抽样检查，完整 QC 由外部工具完成 |
| 变异检测与比对 | GATK、BWA-MEM2、samtools、Picard | 通过 Nextflow/container 执行，不进入 Flask 请求线程 |
| 注释与汇总 | VEP offline/cache、MultiQC | 记录 cache、软件版本和 checksum |

首期 `pysam` 只增强内容级 preflight；未安装时仍可执行 manifest 和路径检查，并明确返回能力缺失警告。`cyvcf2` 不作为首期强制依赖，避免重复引入 VCF 解析栈。

---

## 5. 输入、manifest 与预检查

### 5.1 支持的输入入口

| 输入 | 必需伴随文件 | 起始阶段 | 关键要求 |
| --- | --- | --- | --- |
| paired-end FASTQ.gz | R1+R2；可多 lane | preprocessing | gzip 可验证，lane/read group 完整 |
| BAM | BAI | variant calling | coordinate sorted、header/read group、参考版本可确认 |
| CRAM | CRAI + 对应参考 | variant calling | 可解码、索引和参考 MD5 一致 |
| VCF.gz | TBI | annotation | bgzip、排序、索引、assembly 可确认 |

每个样本只能选择一种入口类型。BAM/CRAM 入口不默认重新比对；如果预处理状态不满足要求，preflight 阻断并提示从 FASTQ 重跑或导入合格的 recalibrated alignment。

### 5.2 manifest 最小字段

```text
patient_id, sample_id, role, matched_normal_id, sex,
library_id, lane, input_type,
fastq_1, fastq_2, bam, bai, cram, crai, vcf, tbi,
reference_bundle_id, capture_bed_id,
capture_kit, capture_kit_version, capture_lot,
dna_source, preservation, umi, expected_depth,
tumor_purity, notes
```

约束：

- `role` 仅允许 `germline`、`normal`、`tumor`；
- tumor-normal 通过 `patient_id + matched_normal_id` 明确配对，不能依赖文件顺序；
- 一个 normal 可按项目规则关联多个同患者 tumor，但必须显式登记；
- lane 是技术分片，不得被识别为生物学样本；
- `sex` 使用平台定义枚举，同时保留 unknown/other，不根据文件名猜测；
- tumor-only 必须记录缺少 matched normal 的原因和用户确认；
- 家系另附版本化 PED，样本 ID 必须与 manifest、VCF/BAM header 一致。

### 5.3 preflight 硬阻断

- 文件不存在、不可读、checksum 失败或索引缺失；
- FASTQ R1/R2、lane 或 read group 不完整；
- 一个样本同时声明多种入口；
- 样本 ID 重复、角色冲突或 tumor-normal 配对不明确；
- BAM/CRAM 未 coordinate sort、header sample 与 manifest 不一致；
- VCF 未 bgzip/排序/index 或无法确认 assembly；
- FASTA、known sites、capture BED、PoN、gnomAD、VEP cache assembly 不一致；
- `chr` 命名、sequence dictionary 或 contig 集不一致；
- capture BED 越界、为空或与 bundle 不兼容；
- 可用磁盘、临时空间或内存低于运行前估算；
- 真实人类数据项目未启用要求的权限与审计策略。

### 5.4 preflight 警告

- tumor-only、无 PoN、未知肿瘤纯度；
- FFPE、高龄样本或未知 DNA 质量；
- capture kit/lot 未登记；
- 实际深度低于项目预期；
- normal 来源可能存在肿瘤浸润或克隆性造血风险；
- 外部 BAM 的预处理工具和参数无法完整追溯。

警告必须进入审批卡、结果页和最终报告，不能只写日志。

---

## 6. 参考资源管理

### 6.1 GRCh38 bundle

首个生产 bundle 必须选择一个确定的参考序列集合，例如经验证的 GRCh38/hs38DH 方案；不能使用“GRCh38 primary + decoy 或 hs38DH”这样的运行时二选一描述。

bundle 至少包含：

- reference FASTA、FAI、sequence dictionary、BWA/BWA-MEM2 index；
- BQSR known sites；
- Mutect2 AF-only germline resource；
- 经验证的 PoN（若有）；
- VEP offline cache 与对应 FASTA；
- ClinVar、gnomAD 等注释资源；
- capture target BED、分析 padding interval、callable territory 定义；
- 每个文件的来源、版本、发布日期、checksum、许可说明和构建命令。

GRCh37/hg19 如果后续支持，必须作为独立 bundle 完成全套 benchmark，不允许在同一项目中混用或通过 liftover 后假定等价。

### 6.2 capture BED

原始厂商 BED、用于 calling 的 padded interval 和用于 QC/报告的 target interval 分开登记。禁止把 padded interval 直接当作 TMB 分母或覆盖率分母。

### 6.3 Panel of Normals

- PoN 必须与参考、capture kit、建库类型和预处理链兼容；
- 首期允许导入已验证 PoN；
- 平台自建以不少于约 40 个兼容 normal 为建设目标，并通过留出正常样本评估伪影抑制效果；
- 不把“gnomAD + 严格过滤”描述为 PoN 的等价替代；
- 每个 PoN 记录来源样本集合、排除规则、workflow版本和 checksum。

---

## 7. 科学分析流程

### 7.1 公共预处理

FASTQ 入口：

```text
raw FastQC
→ 可选 fastp trimming（由接头/文库/UMI信息决定）
→ BWA-MEM2 alignment + 完整 read groups
→ coordinate sort
→ duplicate marking
→ BQSR
→ alignment/coverage/fingerprint QC
→ CRAM/BAM + index
```

原则：

- trimming 不是无条件步骤；默认保留原始 QC，只有检测到接头或文库规范要求时才启用；
- UMI 文库必须进入专门的 UMI-aware 分支，不能先按普通重复序列处理；
- BQSR known sites、reference 和 interval 必须来自同一 bundle；
- 中间清洗 FASTQ 不作为长期交付物，默认在运行验证完成后按保留策略清理；
- CRAM 可作为长期 alignment 交付格式，但必须保留 reference bundle 关联。

#### 7.1.1 胚系与肿瘤分支边界

平台采用“公共预处理、独立 calling contract”，不复制两套 FASTQ 处理代码，也不把两种分析压缩成同一套 caller 和过滤参数。

| 环节 | Germline | Somatic |
| --- | --- | --- |
| 样本设计 | 单个 germline 样本或小家系 | tumor-normal 为默认；tumor-only 为受限模式 |
| 默认 caller | HaplotypeCaller GVCF | Mutect2 |
| 分型模型 | 二倍体 genotype 与联合分型 | 肿瘤纯度、亚克隆和正常细胞混入下的体细胞概率 |
| 专用资源 | germline known sites、群体频率和家系信息 | matched normal、AF-only germline resource、PoN、F1R2 和 contamination evidence |
| 核心结果 | gVCF、genotyped/filtered VCF、GQ/DP/遗传模式 | unfiltered/filtered VCF、Mutect2 stats、tumor/normal AD/DP/VAF、污染和方向性模型 |
| 主要风险 | 覆盖不足、基因型错误、家系关系错误 | 胚系残留、低纯度、FFPE/方向性伪影、样本配对错误 |

平台约束：

- workflow 启动前必须明确选择 `wes_germline` 或 `wes_somatic`，不得在运行中隐式切换 caller；
- 一个 matched normal 可以单独进入 germline run，同时作为同患者 somatic run 的对照，但两次运行具有独立 run ID、过滤状态、结果表和报告；
- 公共 alignment 资产只有在 reference、capture、read group 和预处理版本一致且通过 preflight 时才能复用；
- germline VCF 和 somatic VCF 分别登记 artifact kind 和证据字段，不生成一张混合“所有变异表”；
- `GQ/DP/遗传模式` 与 `tumor-normal AD/VAF/污染证据` 分属不同解释语义，前端筛选和 AI 摘要不得混用；
- somatic 不以固定 VAF 阈值替代 Mutect2 过滤；tumor-only 必须保留胚系残留风险提示。

### 7.2 Germline SNV/InDel

默认流程：

```text
HaplotypeCaller -ERC GVCF（逐样本）
→ GenomicsDBImport（小家系/小样本组）
→ GenotypeGVCFs
→ 经过验证的 variant filtering
→ normalization
→ VEP offline annotation
```

规则：

- 单样本仍保留 gVCF，便于未来在兼容 bundle 下重新联合分型；
- 家系 joint genotyping 只在同参考、同 capture 和兼容预处理样本间执行；
- 小样本不强行套用依赖大队列训练的 VQSR；在 M0/M2 通过 truth set 固化适合该工作流的 caller-specific 或 hard-filter 策略；
- `FilterVariantTranches` 只有在明确生成相应数值评分并验证 tranche 阈值时才能使用；
- DeepVariant 是独立备选 workflow。不得因为存在 GPU 就与 HaplotypeCaller 自动取并集或交集；
- 家系模块可输出 Mendelian error、候选遗传模式和共分离证据，但不自动宣布致病性。

### 7.3 Somatic SNV/InDel

默认流程：

```text
Mutect2（tumor-normal 或 tumor-only）
  + AF-only germline resource
  + compatible PoN（推荐）
  + F1R2 evidence
→ GetPileupSummaries / CalculateContamination
→ LearnReadOrientationModel
→ FilterMutectCalls
→ normalization
→ VEP offline annotation
```

规则：

- matched normal 为默认模式；
- tumor-only 必须二次确认，并在变异表、摘要和报告显著标记“胚系残留和技术伪影风险升高”；
- 无 PoN 的 tumor-only 不作为默认自动运行路径，只能由专家模式显式批准；
- 不设置跨项目通用的 `VAF ≥ 0.05` calling 阈值。caller 使用经过验证的统计过滤；VAF阈值只能作为查看/导出筛选，并保留原始 PASS 状态；
- VerifyBamID2/Somalier 类工具用于样本身份和群体污染辅助检查；somatic filtering 仍使用 Mutect2 工作流的 contamination 估计；
- FFPE和氧化损伤风险必须进入 orientation/context artifact 评估；
- Strelka2 可作为独立验证性 workflow，只有完成预定义合并算法和 truth-set 比较后才考虑 ensemble；
- VarScan2 不进入 v1 默认路径。

### 7.4 标准化与注释

默认 VEP 使用本地 cache/offline 模式，不把个体变异坐标发送到公共服务。输出至少包括：

- normalized chromosome/position/ref/alt；
- gene、transcript、consequence、HGVS c./p.；
- caller FILTER、DP、AD、VAF/GQ 等原始证据；
- gnomAD population AF；
- ClinVar significance、review status 和版本；
- somatic 可选公开驱动基因标签；
- 注释 cache、plugin、FASTA 和命令版本。

COSMIC、OncoKB 等资源按许可和部署场景配置，不在默认镜像中假定可自由再分发。

### 7.5 ACMG 边界

v1 不输出“ACMG 分级表”。允许输出：

- ClinVar 原始分类及 review status；
- 自动可计算的候选证据字段；
- 家系共分离/Mendelian error；
- 供人工复核的 evidence worksheet。

任何 P/LP/VUS/LB/B 结论必须注明规则版本、疾病/基因上下文、证据来源、人工复核人和日期；此能力另行立项。

### 7.6 v1.1 候选能力

以下功能分别建立适用条件、truth set、验收标准后再开放：

| 能力 | 前置条件 |
| --- | --- |
| 肿瘤 CNV | capture-compatible normals；纯度/倍性模型；CNV truth/正交验证 |
| Germline gCNV | 同 capture/建库 cohort，建议至少 30 个兼容样本；独立 CNV benchmark |
| WES SV | 明确标记探索性；限定可检测事件；人工证据回看；不得宣称完整SV检测 |
| MSI | 工具和样本模式适用性检查；每个 capture kit 的位点数和阈值验证 |
| TMB | callable territory、覆盖门槛、PASS规则、胚系过滤、版本和校准固定 |
| 突变特征 | 足够突变数/队列规模；优先已知 signature assignment；输出稳定性指标 |
| ACMG辅助 | phenotype、inheritance、ClinGen规则版本和人工审核工作流 |

---

## 8. QC 与科学门槛

### 8.1 硬阻断与软告警

硬阻断用于输入完整性、身份、参考一致性和资源充分性；测序质量指标通常采用项目/文库特异告警，不使用一个阈值覆盖所有 capture kit、FFPE状态和肿瘤纯度。

### 8.2 核心 QC 指标

| 类别 | 指标 | v1 处理方式 |
| --- | --- | --- |
| 原始 reads | reads数、Q20/Q30、GC、接头、read length | 报告并与项目预期比较 |
| 比对 | mapping rate、proper pairs、insert size | 明显偏离参考/同批样本时告警 |
| 捕获 | on-target rate、mean target depth、≥20×/≥30×覆盖率 | 按 capture kit 和用途配置目标 |
| 文库 | duplicate rate、complexity、coverage uniformity | 不设统一 `≤20%` 硬门槛；结合深度/FFPE解释 |
| 身份 | fingerprint concordance、tumor-normal一致性、sex concordance | 不一致硬阻断 |
| 污染 | germline FREEMIX；somatic contamination estimate | 按模式分别解释；高风险阻断或降级 |
| Germline calls | PASS SNV/InDel数、Ti/Tv、Het/Hom、dbSNP overlap | WES Ti/Tv通常约3.0–3.3，但按target/padding/人群解释 |
| Somatic calls | PASS数、VAF/DP分布、filter原因、orientation/context artifacts | 不以单一VAF阈值替代caller过滤 |
| 家系 | Mendelian errors、亲缘关系、缺失率 | 关系不符硬阻断家系解释 |

肿瘤 `100×`、normal `60×` 可作为常见项目规划起点，但不是跨项目验收定律。最终门槛由样本类型、肿瘤纯度、目标 LoD 和 capture 设计决定，并写入项目审批卡。

### 8.3 报告中的 QC 状态

每个样本输出 `PASS / WARN / FAIL / NOT_APPLICABLE`，同时给出：

- 指标值、参考范围及范围来源；
- 是否影响 calling 或仅影响解释；
- 用户是否覆盖警告、覆盖理由和操作者；
- 受影响的下游结果。

---

## 9. 结果与平台展示

### 9.1 必需交付物

| 类别 | 工件 |
| --- | --- |
| 设计 | 固化 manifest、preflight 报告、审批记录 |
| QC | MultiQC、结构化 metrics JSON/CSV、样本身份/污染摘要 |
| Alignment | CRAM/BAM + index（按保留策略） |
| Germline | gVCF + index、filtered VCF + index、统计表 |
| Somatic | unfiltered/filtered VCF + index、Mutect2 stats、contamination/orientation模型 |
| 注释 | annotated VCF + index、可查询 Parquet/TSV、字段字典 |
| 报告 | HTML/Markdown、参数、软件、参考、限制、manifest |
| 运行 | Nextflow trace/timeline/report、stdout/stderr、container/reference provenance |

`FASTQ/BAM/VCF` 不应仅通过扩展 `VALID_RESULT_FILE_TYPES` 当作普通展示文件。它们应登记为带类型、索引、checksum、父子关系、敏感级别和保留策略的 `DataAsset/WorkflowArtifact`；`ResultFile` 只保存兼容展示入口。

### 9.2 WES 结果页

至少包含：

- 样本、角色、配对、入口类型和 reference/capture 状态；
- QC gate 和受影响结果；
- 可筛选变异表：gene、consequence、HGVS、DP、AD、VAF/GQ、FILTER、population AF；
- 每条变异回链到原始 VCF 行、caller和注释版本；
- 下载工件、index、manifest 和 provenance；
- tumor-only、低覆盖、无PoN等限制的显著提示；
- 固定文案：“科研用途，不构成临床诊断或治疗建议”。

默认前端只加载分页/过滤后的表格，不把全量 VCF 读入浏览器。

### 9.3 AI 工具

只读工具：

- `inspect_wes_manifest`
- `inspect_wes_preflight`
- `get_workflow_run_status`
- `inspect_wes_qc`
- `query_variant_summary`
- `get_wes_artifact_provenance`

写工具：

- `create_wes_manifest`
- `run_wes_workflow`
- `resume_wes_workflow`
- `cancel_wes_workflow`

所有写工具要求确认。审批卡显示 workflow版本、输入样本和角色、参考bundle、capture BED、关键参数、资源估算、主要工件、数据外发策略、阻断项和警告。

外部模型默认不可见 FASTQ/BAM/CRAM、全量 VCF、绝对路径和直接身份信息；只接收结构化摘要。变异坐标只有在本地模型策略或用户明确授权的小范围查询中可见。

---

## 10. 验证与验收

### 10.1 验证数据

| 路径 | 数据 | 用途 |
| --- | --- | --- |
| Germline | GIAB HG001/HG002 等具有 high-confidence truth 的外显子数据 | SNV/Indel precision/recall、入口和参考验证 |
| Somatic | HCC1395 tumor-normal truth set | somatic SNV/Indel、VAF分层和过滤验证 |
| 工程测试 | 小型公开/合成区间数据 | 快速 E2E、失败路径、升级回归 |

合成 1 Mb 数据只能用于工程冒烟测试，不能替代科学性能验证。

### 10.2 比较口径

- 使用规范化、haplotype-aware 工具（如 hap.py）比较；
- 评估区域为 truth high-confidence BED、capture BED 和平台 callable BED 的交集；
- 只对定义清楚的 PASS/filtered callset 评分；
- SNV 与 Indel 分开；somatic 再按 VAF、深度、纯度和样本模式分层；
- 报告 TP、FP、FN、precision、recall、F1 和不可评估区域；
- tumor-normal 与 tumor-only 独立验证，不能共享一个结论。

### 10.3 初始科学验收下限

以下为 v1 工程上线下限，M0 可在不降低科学要求的前提下根据 truth/capture 组合形成 ADR：

| 路径/层级 | Precision | Recall | 备注 |
| --- | ---: | ---: | --- |
| Germline SNV | ≥ 0.99 | ≥ 0.97 | high-confidence/capture/callable 交集 |
| Germline Indel | ≥ 0.98 | ≥ 0.90 | 单独报告长度与困难区域 |
| Somatic SNV，VAF ≥ 0.10 | ≥ 0.90 | ≥ 0.90 | HCC1395 tumor-normal 可评估区域 |
| Somatic Indel，VAF ≥ 0.10 | ≥ 0.80 | ≥ 0.80 | 低计数时同时报告置信区间 |

低 VAF 层单独报告性能曲线，不用一个阈值掩盖 LoD。若数据集、capture 或 truth 区域不足以支持上述评价，状态为“未验证”，不能以“流程跑通”替代。

### 10.4 工程验收矩阵

必须覆盖：

1. FASTQ、BAM/CRAM、VCF 三种入口；
2. germline 单样本、trio；
3. tumor-normal、tumor-only、无PoN tumor-only；
4. 多lane合并、重复sample ID、错误配对；
5. 缺index、错误排序、header/manifest冲突；
6. GRCh38/GRCh37混用、`chr`命名冲突、错误capture BED；
7. Web重启后恢复监控、失败后resume、取消进程组；
8. 重复点击幂等、防止重复工作流；
9. 必需工件缺失、登记失败和checksum改变必须使任务失败；
10. 权限、跨项目路径、下载授权和审计测试。

### 10.5 回归策略

- 固定 workflow release、容器 digest、reference bundle 和参数；
- 比较 normalized VCF 内容、分层指标和必需工件，不比较时间戳或压缩文件原始字节；
- 相同版本重复运行的 normalized callset 应一致；
- 升级后 precision/recall 下降超过预设容差或出现未解释变异差异时阻断发布；
- 所有允许差异进入版本化变更记录。

---

## 11. 安全、隐私与合规

真实人类 WES 数据试点前必须完成：

- 用户认证、项目级 RBAC 和最小权限；
- 上传、查看、查询、下载、删除和AI调用的审计日志；
- CSRF/API token、严格 CORS 和幂等键；
- 传输加密、受控存储、备份策略和密钥管理；
- 数据保留期限、用户可见清理策略和可验证删除；
- 参考目录只读、容器最小挂载、禁止任意 shell 和任意 host path；
- 敏感工件下载授权，避免仅凭可猜测 URL；
- 外部 LLM 默认 `summary_only`，禁止原始序列和全量变异外发；
- 报告脱敏和直接身份信息最小化。

如果这些能力未完成，只允许使用公开测试数据，不允许真实样本试运行。

---

## 12. 资源、性能与保留策略

### 12.1 运行估算

资源根据 FASTQ 总碱基数、样本数、入口阶段和 capture 区域动态估算。典型单样本 WES 可用作容量规划参考，但不能作为 SLA：

| 项目 | 初始规划值 |
| --- | --- |
| 单样本 preprocessing + calling | 8–16 vCPU、24–48 GB RAM |
| 典型运行时间 | 约 4–12 小时，取决于数据量、存储和scatter配置 |
| 原始 FASTQ | 常见约 15–30 GB/样本 |
| CRAM/BAM与中间文件 | 依深度和保留策略变化 |
| 临时/工作空间 | 运行前预留预计最终工件的 3–5 倍，并设置硬配额 |
| 参考资源 | 由实际 bundle 清单计算，不使用笼统固定数值 |

### 12.2 调度和配额

- WES 与 Python 模块使用独立并发池；
- Nextflow进程声明 CPU/内存/时限并由 executor 限制；
- 单机默认并发运行数从1开始，经压测后调整；
- 磁盘不足在启动前阻断，运行中持续监测；
- 支持取消、超时、失败诊断和安全resume。

### 12.3 保留策略

默认长期保留：manifest、报告、filtered VCF/gVCF及index、注释表、QC、provenance、必要alignment。默认可清理：trimmed FASTQ、scatter临时文件、Nextflow work缓存和未选择的中间VCF。

清理只能针对明确 `run_id` 和已登记工件执行；先生成清单，保留审计记录，不使用项目根目录递归删除。

---

## 13. P0–P3 实施阶段

| 阶段 | 当前范围 | 交付物 | 退出条件 |
| --- | --- | --- | --- |
| P0 合同与安全底座 | germline/somatic 分离合同、manifest、路径边界、配对与 tumor-only 门禁 | workflow registry、preflight、版本化 manifest、安全测试 | 错误角色、配对、路径和 capture 输入可硬阻断 |
| P1 执行器与工程冒烟 | 固定 Sarek/Nextflow、typed params、poll/cancel/resume/log、FASTQ/BAM/CRAM/VCF 入口 | LaunchSpec、不可变 launch bundle、三条官方微型数据冒烟记录 | 固定版本执行链跑通且退出码、日志、参数可追溯 |
| P2 平台补充与结果体验 | 本地 reference/capture catalog、每试剂盒 profile、启动资源门禁、工件收集、项目结果页 | reference 管理 CLI/API、结果 API/UI、checksum/provenance | 真实 bundle 与试剂盒资源 validated；必需 VCF/index/MultiQC 自动归档 |
| P3 科学验证与受控试点 | GIAB/HCC1395、评估 BED 交集、hap.py/som.py、故障/磁盘/安全演练、少量实验室试点 | 版本化 validation JSON/报告、回归基线、试点复盘 | 达到 10.3 阈值且无未解释阻断项；真实数据安全门槛通过 |

状态定义必须保守：代码和测试框架完成不等于科学验证完成；官方 tiny 数据跑通
只能记为工程冒烟。没有实际 capture kit、GIAB/HCC1395 输入或预定义评估区域时，
对应 P2/P3 项标记为 `blocked_by_data`，不能写成“已通过”。

CNV/SV/MSI/TMB/signature/ACMG辅助作为 v1.1 独立里程碑，不阻塞 SNV/InDel v1，也不能未经验证提前出现在生产报告中。

---

## 14. 风险与对策

| 风险 | 影响 | 对策 |
| --- | --- | --- |
| 复用现有线性Pipeline | 多输入/分支无法正确表达 | 采用WORKFLOW_REGISTRY和通用工件契约 |
| Web进程承担重计算 | 重启丢监控、资源失控 | Nextflow executor；持久化external run ID |
| reference/capture混用 | 系统性假阳性/假阴性 | ReferenceAsset、checksum、preflight硬阻断 |
| tumor-only胚系残留 | 假阳性和错误解释 | 二次确认、PoN、群体资源、降低证据等级、独立验证 |
| caller ensemble无规则 | 结果不可解释 | v1单一默认caller；备选caller独立benchmark |
| FFPE/低纯度 | 伪影和低VAF漏检 | 记录前分析元数据、artifact model、VAF分层性能 |
| 大文件占满磁盘 | 任务失败、平台不可用 | 启动前估算、硬配额、CRAM、run级清理 |
| 注释库更新漂移 | 同位点解释变化 | 版本化cache、checksum、报告中保留历史版本 |
| 自动ACMG被误用 | 科研结果被当作诊断 | v1禁用自动分级，仅输出证据字段和明确免责声明 |
| 人类遗传数据泄露 | 严重隐私风险 | 安全门槛前置；默认不向外部LLM发送变异数据 |

---

## 15. Definition of Done

- [ ] 上位路线图与本文档不存在执行器、数据模型和v1范围冲突；
- [ ] `wes_germline`、`wes_somatic`、`wes_annotate_only` 注册到 `WORKFLOW_REGISTRY`；
- [ ] FASTQ、BAM/CRAM、VCF入口均有正向和负向E2E测试；
- [ ] manifest、reference、capture、样本角色和配对preflight可硬阻断错误输入；
- [ ] 外部任务支持poll、cancel、resume和Web重启恢复；
- [ ] 所有必需工件有checksum、index、lineage、provenance和保留策略；
- [ ] Germline与somatic达到第10.3节分层验收下限；
- [ ] tumor-only有独立验证和明显限制说明；
- [ ] 结果页和报告不出现未经验证的CNV/SV/MSI/TMB/signature/ACMG结论；
- [ ] AI写操作需要确认，外部模型默认只接收结构化摘要；
- [ ] 真实数据试点前认证、RBAC、下载审计和数据删除策略通过安全测试；
- [ ] 用户手册、管理员部署手册、reference bundle文档和科学验证报告齐全。

---

## 16. 参考资料

- nf-core/sarek 3.10.0 usage：<https://nf-co.re/sarek/3.10.0/docs/usage>
- nf-core/sarek 3.10.0 parameters：<https://nf-co.re/sarek/3.10.0/parameters>
- GATK Mutect2：<https://gatk.broadinstitute.org/hc/en-us/articles/360036713111-Mutect2>
- GATK FilterMutectCalls：<https://gatk.broadinstitute.org/hc/en-us/articles/360041849811-FilterMutectCalls>
- GATK GenomicsDBImport：<https://gatk.broadinstitute.org/hc/en-us/articles/360047216891-GenomicsDBImport>
- GATK GermlineCNVCaller：<https://gatk.broadinstitute.org/hc/en-us/articles/4414586797467-GermlineCNVCaller>
- GATK germline callset QC：<https://gatk.broadinstitute.org/hc/en-us/articles/360035531572-Evaluating-the-quality-of-a-germline-short-variant-callset>
- GA4GH variant benchmarking：<https://www.ga4gh.org/product/variant-benchmarking-tools/>
- ACMG/AMP sequence variant interpretation：<https://pubmed.ncbi.nlm.nih.gov/25741868>
- ClinGen variant classification guidance：<https://www.clinicalgenome.org/tools/clingen-variant-classification-guidance/>
