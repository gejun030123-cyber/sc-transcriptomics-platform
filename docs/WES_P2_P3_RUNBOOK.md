# WES P2/P3 实验室运行手册

本文只面向实验室内部科研部署。命令通过平台管理员在服务器执行；普通用户只在
项目 WES 页面选择已经登记的资源，不直接填写任意 shell 命令。

## 1. 多试剂盒处理规则

- 每个“厂商 + 产品 + 版本”使用独立、稳定的 `capture_kit_id`。
- calling BED 必填；vendor 原始 BED 和 QC BED 可选。三者不互相覆盖。
- 同一 run 只选择一个 profile；tumor/normal 应使用同一试剂盒。不同试剂盒
  不自动求交集，首期拆成独立 run。
- 官方 tiny BED 只能登记为 `test_only`。真实 BED 检查来源、assembly、contig、
  checksum 和区间后才可设为 `validated`。

登记官方测试或真实试剂盒：

```bash
python scripts/wes_admin.py register-capture \
  --id twist-human-core-exome-v2 \
  --name "Twist Human Core Exome" \
  --version v2 \
  --assembly GRCh38 \
  --calling-bed /允许的数据目录/references/twist_v2.calling.bed \
  --vendor-bed /允许的数据目录/references/twist_v2.vendor.bed \
  --qc-bed /允许的数据目录/references/twist_v2.qc.bed \
  --status test_only
```

完成公开 benchmark 和人工复核后，用返回的 asset ID 提升状态：

```bash
python scripts/wes_admin.py set-status ref_xxxxxxxxxxxxxxxx validated
```

### 用户上传自己的 BED

项目内用户可以通过 WES 页面上传自己的 `calling BED`，或调用同一个 JSON 接口：

```bash
curl -X POST \
  -F capture_kit_id=lab-exome-v1 \
  -F name='实验室实际外显子试剂盒' \
  -F version=v1 \
  -F assembly=GRCh38 \
  -F calling_bed=@/本地路径/lab-exome-v1.bed \
  http://平台地址/api/projects/<project_id>/wes/capture-kits/upload
```

上传文件只允许 `.bed`/`.bed.gz`，平台会检查文件大小、BED 区间、contig 命名、
自然染色体排序并计算 checksum。成功后状态固定为 `test_only`，不会自动启动流程；
管理员核对试剂盒来源、版本、assembly、实验室实际捕获方案和试运行结果后，才可用
返回的 `asset_id` 执行 `set-status ... validated`。不同版本必须使用不同的稳定
`capture_kit_id`，避免覆盖既有 profile。

### 用户端创建 sample manifest

项目 WES 页面中的“创建 WES 样本运行”向导会把表单自动转换为 manifest，并调用
`/api/projects/<project_id>/wes/preflight` 进行预检。用户不需要手写 JSON。测序公司
提供的 SampleSheet/样本信息表可作为填写依据。向导支持上传 FASTQ R1/R2、BAM+BAI
或 CRAM+CRAI，并将文件登记到项目目录；也可使用已放入项目目录或管理员配置的
`WES_SOURCE_ROOTS` 的文件路径。SRA 需先转换为 FASTQ，不能直接上传。

页面也提供“上传 SRA 并转换为 FASTQ”：管理员安装 SRA Toolkit 后，平台先检查磁盘
空间，再异步执行 `fasterq-dump --split-files`，检测 paired-end R1/R2、压缩为
`fastq.gz`，并将两个 FASTQ 登记为带 checksum 的 `data_assets`。转换任务可通过
`/api/projects/<project_id>/wes/sra/jobs/<job_id>` 查询；失败时保留 job 日志和原始
SRA，不能绕过配对检查直接进入 manifest。默认磁盘门禁至少保留 50 GB，并按 SRA
大小的 3 倍估算转换空间，可用 `.env` 中的 `WES_SRA_MIN_FREE_GB`、
`WES_SRA_DISK_FACTOR`、`WES_SRA_THREADS` 和 `WES_SRA_FASTERQ_BIN` 调整。

向导需要用户确认四件事：分析流程（胚系或肿瘤-正常）、已审核的 capture kit、
reference bundle，以及每个样本的文件路径和角色。肿瘤-正常 run 必须将 tumor 的
`matched_normal_id` 指向同一 `patient_id` 下的 normal。预检返回 `manifest_id` 后，
用户可以先生成可审阅的 LaunchSpec；只有管理员打开执行器且用户再次确认，才会提交
Nextflow。

## 2. 本地 GRCh38 bundle

Sarek 3.10.0 的 `GATK.GRCh38` 目录应放在
`WES_NEXTFLOW_IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/`。平台至少登记并校验
FASTA、FAI、DICT、dbSNP、gnomAD germline resource；somatic 还需要 PoN。
每个文件单独登记到同一个 `bundle_version`：

```bash
python scripts/wes_admin.py register-reference \
  --assembly GRCh38 --bundle gatk-grch38-lab-v1 \
  --type fasta --path /允许的数据目录/references/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  --status registered
```

重复命令登记 `fai`、`dict`、`dbsnp`、`germline_resource`、`pon`；完成联合验证后
逐项提升为 `validated`。VEP cache 另准备 `homo_sapiens/116_GRCh38`，并用一个
文本 manifest（包含 cache 版本、归档来源和归档 SHA-256）登记为
`vep_cache_manifest`。

`.env` 最终填写：

```dotenv
WES_NEXTFLOW_IGENOMES_BASE=/允许的数据目录/references/igenomes
WES_NEXTFLOW_VEP_CACHE=/允许的数据目录/references/vep_cache
WES_NEXTFLOW_PON=/允许的数据目录/references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz
WES_NEXTFLOW_GERMLINE_RESOURCE=/允许的数据目录/references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz
WES_REQUIRE_VALIDATED_REFERENCES=true
WES_MIN_FREE_GB=50
```

这些路径为空时允许准备和审阅 run，但禁止真正启动。不要在目录尚未完整时填入
占位路径。

本机内部部署的实际根目录为
`/home/oelab/AnaData/GJ/02_platform_runtime/wes`；当前 `gatk-grch38-lab-v1`
中的 FASTA/FAI/DICT、dbSNP 138 及其索引、gnomAD、PoN 及其索引、VEP 116
manifest 已完成文件完整性检查并登记为 `validated`，三种 WES workflow 的
reference readiness 均为通过。官方 Illumina、Twist、IDT capture profile 仅作 `test_only` 示例；
用户自己的 BED 通过上面的上传接口进入同样的审核流程。由于真实样本试运行由实验室
人员在网页完成，管理员配置阶段保持 `WES_EXECUTOR_ENABLED=false`，不要把资源
登记误认为科学结果已验证。

VEP archive 解压后生成 provenance manifest：

```bash
python scripts/wes_admin.py vep-manifest \
  --cache-root /home/oelab/AnaData/GJ/02_platform_runtime/wes/references/vep_cache \
  --output /home/oelab/AnaData/GJ/02_platform_runtime/wes/references/vep_cache/vep_116_GRCh38.manifest.json \
  --release 116 \
  --archive /home/oelab/AnaData/GJ/02_platform_runtime/wes/references/vep_cache_archive/homo_sapiens_vep_116_GRCh38.tar.gz
```

再将该 manifest 以 `vep_cache_manifest` 登记到同一个 bundle，并把 manifest
metadata 中的 `cache_root` 设为 `.env` 的 `WES_NEXTFLOW_VEP_CACHE`。

## 3. P3 benchmark

评估区域必须是 truth confident、当前 capture 和平台 callable 三者交集：

```bash
python scripts/wes_admin.py evaluation-bed \
  --bed /refs/GIAB/HG002_GRCh38_confident.bed.gz \
  --bed /refs/capture/twist_v2.calling.bed \
  --bed /results/HG002.callable.bed \
  --output /results/validation/HG002.evaluation.bed
```

Germline 使用 hap.py，somatic 使用 som.py；两者都应固定版本和容器 digest。

```bash
python scripts/wes_admin.py benchmark \
  --mode germline \
  --truth-vcf /refs/GIAB/HG002_GRCh38_truth.vcf.gz \
  --query-vcf /results/HG002.filtered.vcf.gz \
  --evaluation-bed /results/validation/HG002.evaluation.bed \
  --reference-fasta /refs/GRCh38/Homo_sapiens_assembly38.fasta \
  --output-prefix /results/validation/HG002.sarek-3.10.0
```

输出 `*.validation.json` 保存输入 checksum、固定参数、stderr/stdout 尾部、
SNV/Indel 指标和门槛结论。Somatic 的 VAF≥0.10 结果还要传入一个 JSON 证据
文件：至少包含 `stratum`、`truth_vcf_checksum`、`query_vcf_checksum` 和
`filter_expression`，且 checksum 必须对应本次输入；命令同时使用
`--stratum vaf_ge_0.10 --stratum-evidence evidence.json`。缺少这份证据时，
即使总指标较高，状态仍保持 `not_validated`。

## 4. 开启执行器前检查

1. 实际 capture profile 与 bundle 所有必需资产均为 `validated`。
2. `.env` 四个本地资源路径均存在，且都位于 `WES_SOURCE_ROOTS`。
3. GIAB 与 HCC1395 达到方案 10.3 的 SNV/Indel 分层阈值。
4. 登录或受控内网、项目权限、下载审计、备份和保留策略已落实。
5. 以上完成后才把 `WES_EXECUTOR_ENABLED=true`，并先执行一个公开数据 run。
