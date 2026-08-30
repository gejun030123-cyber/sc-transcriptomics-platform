# WES 执行环境冒烟测试记录（2026-08-24）

## 结论

Nextflow + Docker + nf-core/sarek 的本地执行链已在当前主机跑通
germline BAM/HaplotypeCaller、tumor-normal CRAM/Mutect2 和 VCF/VEP annotate-only
三个冒烟测试，生成 VCF、index、QC 统计和 MultiQC 报告。

本结论只证明软件链与容器执行环境可用，不代表真实 WES 生产方案已通过
科学验证。当前 `WES_EXECUTOR_ENABLED` 继续保持 `false`。

## 固定版本和入口

| 项目 | 值 |
| --- | --- |
| Nextflow | 26.04.6 |
| Pipeline | nf-core/sarek 3.10.0 |
| Sarek revision | `8ccac7ad37b05dd792447763bf9671b719824587` |
| Profile | `docker` |
| Genome key | `GATK.GRCh38` |
| Workflow | `wes_germline` / HaplotypeCaller |
| Input | Sarek 3.10.0 官方 mapped single BAM test samplesheet |
| Start step | `variant_calling` |

Nextflow 26 会将命令行 pipeline 参数解析为字符串，所以平台使用
`parameters.json` 通过 `-params-file` 传入布尔值 `"wes": true`，不再传入
裸 `--wes`。FASTQ 从 `mapping` 开始，BAM/CRAM 从 `variant_calling`
开始，避免对已比对数据重复执行 mapping。

## 结果证据

### Germline HaplotypeCaller

- Pipeline 状态：`Pipeline completed successfully`。
- 任务统计：16 succeeded，0 failed，峰值 9 CPU / 34 GB 内存。
- Filtered VCF：5 条变异记录，其中 2 条 `PASS`。
- VCF SHA-256：`74cd4abc54c7ef8d0f7e0b35a5f479c2986eee2f4119ce036a4084e606f059dc`。
- TBI SHA-256：`5d3ef11f38112405f92c660e1d13110ff6dfd48660971ca0f16ed9441466b8f5`。
- MultiQC SHA-256：`a7fad439965d1f8510c7718f0c979b73e575180a47c56412820f877a745c4c11`。

本机结果目录：

```text
/home/oelab/data/wes_smoke/germline_haplotypecaller_bam/results
```

本机执行日志：

```text
/home/oelab/data/wes_smoke/germline_haplotypecaller/logs/nextflow-bam.log
```

### Somatic tumor-normal Mutect2

- Pipeline 状态：`Pipeline completed successfully`。
- 任务统计：15 succeeded，0 failed，峰值 12 CPU / 32 GB 内存。
- 输入：公开 chr21 normal/tumor recalibrated CRAM test pair。
- Mutect2 VCF：40 条测试变异记录。该官方微型场景没有 PoN 和 germline
  resource，所有记录的 FILTER 为 `.`，不得用于科学性能判断。
- VCF SHA-256：`306af2bda5d3e9a4a8f632d14627124748fc833dee01da95d149be1e46767efd`。
- TBI SHA-256：`649887dc1aea9681571b69d26b24db9ca0edf0e994991db1e37beecda9ca63d7`。
- MultiQC SHA-256：`173d1f617512059ecfb6a18b17b4dec73859ad28402f7de7064adabc9c1a73c5`。

本机结果与日志：

```text
/home/oelab/data/wes_smoke/somatic_mutect2/results
/home/oelab/data/wes_smoke/somatic_mutect2/logs/nextflow.log
```

### VCF annotate-only / VEP

- Pipeline 状态：`Pipeline completed successfully`。
- 任务统计：2 succeeded，0 failed，峰值 4 CPU / 15 GB 内存。
- 输入：Sarek 官方空记录微型 VCF，因此产物为 0 变异记录。
- 输出头已包含 Ensembl VEP `CSQ` 及 LoFTEE 字段，证明 VEP 执行和
  VCF/TBI 产物链可用。
- VCF SHA-256：`5bb108fe8bcef4fe19514235275324957fc3e96e172c46ca59500a28f7aff4a0`。
- TBI SHA-256：`5c9e74e847b0ca3ca72ab2bcc803ac43efaff23fa701271af7d9208df054c08e`。
- MultiQC SHA-256：`17714c2cb165e695bbc43fa8774212cc67bfc4e86f42ff66e422b9d68047e40f`。
- 官方 CI 使用的 VEP 测试 cache 为微型非人类配置，不代表 GRCh38
  人类注释资源已验证。

本机结果与日志：

```text
/home/oelab/data/wes_smoke/annotate_vep/results
/home/oelab/data/wes_smoke/annotate_vep/logs/nextflow-resume.log
```

## 已发现并修复的平台问题

1. Nextflow 26 下裸 `--wes` 被 nf-schema 判定为字符串，现改为
   服务端固定的 JSON 类型参数。
2. BAM/CRAM 之前错误从 `mapping` 开始，现改为 `variant_calling`。
3. 之前只要有 `capture_bed_id` 就能准备 calling run；现在 germline/
   somatic calling 必须提供真实、存在且非符号链接的 `capture_bed_path`。
4. LaunchSpec 重新启动时会校验 `parameters.json` 与服务端固定参数
   完全一致，避免审阅后被替换。
5. P2 已增加按试剂盒独立登记的 capture profile、reference bundle
   readiness、运行前 checksum/assembly/path 联合门禁；不同用户可在 manifest
   中选择不同 `capture_bed_id`，不会修改全局 BED。
6. P2 已增加必需 VCF、TBI、MultiQC 与 pipeline provenance 的受限收集、
   checksum 登记、项目级结果页及下载时 checksum 复核。
7. 已登记 nf-core 官方 tiny BED profile `nfcore-sarek-tiny-test`，状态为
   `test_only`、`production_allowed=false`；它只用于工程验证，不能替代真实
   capture kit。
8. P3 已实现 truth confident BED ∩ capture BED ∩ callable BED、固定
   hap.py/som.py 命令、summary 解析、验收阈值和 validation JSON；重启、磁盘
   不足、进程组取消、缺失工件故障路径已有自动测试。

## 不得开启真实样本执行的当前阻断项

1. 尚未登记并验证实验室实际 capture kit 对应的 BED；当前只有
   `test_only` 官方微型 profile。
2. 尚未建立本地 GRCh38 reference bundle manifest，也未对 FASTA、known
   sites、capture BED 的 assembly/contig/checksum 做联合验证。
3. 尚未执行 GIAB germline 准确性验证和 HCC1395 tumor-normal 验证。
4. 带经验证 PoN/germline resource 的 somatic 过滤、GRCh38 人类 VEP
   cache 注释和 FASTQ 完整正向 E2E 仍需分别验证。
5. 当前主机尚未安装/固定 hap.py 与 som.py 可执行环境，也没有下载完整
   GIAB HG002 和 HCC1395 WES benchmark 输入；P3 科学指标仍为
   `not_validated`。

满足以上阻断项并记录 checksum 后，才可将 `WES_EXECUTOR_ENABLED`
改为 `true`。
