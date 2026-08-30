# WES P2/P3 执行状态（2026-08-25）

状态只按 P0–P3 汇报。`engineering_complete` 表示代码与小型工程测试完成，
不等于真实 WES 科学验证通过。

## P2：平台补充与结果体验

| 项目 | 状态 | 证据/剩余项 |
| --- | --- | --- |
| 不同 capture kit 独立 profile | engineering_complete | 每个 profile 含独立 ID、版本、assembly、BED 角色、checksum 与状态 |
| 官方测试 BED | test_only | `nfcore-sarek-tiny-test` 已登记；35,000 bp，`production_allowed=false` |
| 真实试剂盒 | test_only | 已下载并登记 3 个官方 GRCh38 BED；仍需与实验室实际厂商/产品/版本核对后才能 validated |
| 本地 GRCh38 bundle 门禁 | validated | `gatk-grch38-lab-v1` 的 germline、somatic、annotate-only readiness 均为 `True []` |
| 本地资源实体 | validated | FASTA/FAI/DICT、Broad dbSNP 138/index、gnomAD/index、PoN/index、VEP 116 manifest 已落盘并完成 checksum 校验 |
| 结果工件收集 | engineering_complete | 主 VCF、TBI、MultiQC、pipeline provenance；缺失必需工件会使 run 失败 |
| WES 结果页 | engineering_complete | `/projects/<pid>/wes`；项目内下载复核 checksum |
| 用户端 sample manifest 向导 | engineering_complete | WES 页面按流程、capture kit、样本角色和输入路径自动生成 manifest；预检后才可准备 LaunchSpec |
| 用户端 SRA 转 FASTQ | engineering_complete | `.sra` 上传、磁盘门禁、异步 `fasterq-dump`、R1/R2 配对检查和 FASTQ data asset 登记 |

## P3：科学验证与受控试点

| 项目 | 状态 | 证据/剩余项 |
| --- | --- | --- |
| 评估区域生成 | engineering_complete | truth confident ∩ capture ∩ callable，空交集和 contig 冲突阻断 |
| benchmark contract | engineering_complete | germline 固定 hap.py；somatic 固定 som.py；参数数组、`shell=False` |
| 指标与门槛 | engineering_complete | SNP/Indel precision、recall、F1；按方案 10.3 判断 |
| Somatic VAF 分层 | engineering_complete | 只有绑定 truth/query checksum 的分层证据才允许判定 validated |
| 故障演练 | engineering_complete | Web 重启标 interrupted、`-resume`；磁盘门禁；进程组取消；缺失工件失败 |
| GIAB HG002 科学结果 | pending_web_validation | 平台不代跑；由实验室人员从网页提交公开数据后复核 |
| HCC1395 科学结果 | pending_web_validation | 平台不代跑；由实验室人员从网页提交 tumor-normal 数据后复核 |
| 真实样本试点 | pending_web_validation | 执行器继续关闭；网页试运行通过并完成人工复核后再启用 |

当前结论：P2/P3 的平台代码、资源 catalog 和可重复验证框架已经完成；资源 readiness
已通过。P3 的科学性能验收留给实验室人员从网页执行，因此 `.env` 继续保持
`WES_EXECUTOR_ENABLED=false`。

## 已补充的公开 GRCh38 capture BED（2026-08-25）

以下文件已写入实验室 WES references，并通过平台登记为 `test_only`。它们只用于
开发、公开数据试跑和候选 profile 对照，不能替代实验室实际使用的 kit manifest。

| profile | 文件 | SHA-256 | 检查结果 |
| --- | --- | --- | --- |
| `illumina-exome-panel-v1.2` | `capture/official_grch38/illumina-exome-panel-v1.2/Illumina_Exome_TargetedRegions_v1.2.hg38.bed` | `830532d6f8962cd594ad748887b68f8e357ec565abe3c659bfedaf655d0c39b2` | 213,726 intervals，排序通过 |
| `twist-human-core-exome` | `capture/official_grch38/twist-human-core-exome/Twist_Exome_Core_Covered_Targets_hg38.bed` | `d7bafeb53f8130b10724359425e900b91ed4d83ae7488c7a7404d96726c86483` | 192,262 intervals，排序通过 |
| `idt-xgen-exome-hyb-v2` | `capture/official_grch38/idt-xgen-exome-hyb-v2/xGen_Exome_Hyb_Panel_v2_targets_hg38.bed` | `9b18f157033c49380e146ab370258976aa0eaf2a48e4f466c05a5e9f4e41df3a` | 197,769 intervals，自然染色体排序通过 |

来源分别为 Illumina、Twist Bioscience 和 IDT 官方产品页面；平台按自然染色体顺序
检查（chr9 在 chr10 之前），原始文件不被改写。经过实验室版本确认后再提升状态。
