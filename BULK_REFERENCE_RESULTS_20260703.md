# Bulk RNA-seq 两个参考文件运行结果

生成时间: 2026-07-03

输出目录: `bulk_reference_output_codex_20260703/`

本次运行核心流程:

```text
bulk_qc -> bulk_normalize -> bulk_pca -> bulk_deg -> bulk_heatmap
```

运行参数:
- DEG 方法: `t-test`
- 标准化方法: `deseq2`
- 分组列: `_auto_group`
- 热图基因来源: `deg`
- 可选模块: 本次未运行 `bulk_enrichment`, `bulk_timecourse`, `bulk_deg_integration`

## 总体结论

| 数据集 | 输入文件 | 状态 | 样本数 | QC 后基因数 | Normalize 后基因数 | PCA PC1/PC2 | DEG 比较数 | DEG 热图基因数 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Dataset A | `/home/oelab/data/GJ/all.fpkm_anno.xls` | completed | 12 | 61,806 | 19,431 | 25.91% / 18.32% | 3 | 50 |
| Dataset B | `/home/oelab/data/GJ/all.genes.expression.anno.xls` | completed | 24 | 15,256 | 15,040 | 73.32% / 3.22% | 3 | 0 |

说明:
- Dataset A 在默认阈值下检出多组显著差异基因。
- Dataset B 在默认阈值 `padj < 0.05` 且 `|log2FC| >= 1` 下 3 个比较均未检出显著上下调基因，因此 DEG 热图基因数为 0。这是统计结果，不是运行失败。
- Dataset B 的 QC 检测到一个离群样本提示: `Ctr_B_2`，但没有自动移除样本。

## Dataset A 结果

项目目录:

```text
bulk_reference_output_codex_20260703/dataset_a_core/
```

报告:
- JSON: `bulk_reference_output_codex_20260703/dataset_a_core/bulk_reference_report.json`
- 中文报告: `bulk_reference_output_codex_20260703/dataset_a_core/中文运行报告.md`

模块状态:

| 模块 | 状态 |
|---|---|
| bulk_qc | completed |
| bulk_normalize | completed |
| bulk_pca | completed |
| bulk_deg | completed |
| bulk_heatmap | completed |

差异比较:

| 比较 | 上调基因数 | 下调基因数 | Top 5 基因 |
|---|---:|---:|---|
| hmc3-vs-ctrl | 59 | 100 | CCN4, PARM1, COL1A2, AFP, FDXR |
| moclel-vs-ctrl | 37 | 73 | COL1A2, DCN, FN1, COL3A1, SPP1 |
| rapa-vs-ctrl | 73 | 73 | GC, SULF2, SLC31A1, HMGN2, CETN2 |

主要输出:
- DEG 总表: `bulk_reference_output_codex_20260703/dataset_a_core/results/bulk_deg_all_comparisons.csv`
- PCA 图: `bulk_reference_output_codex_20260703/dataset_a_core/plots/bulk_pca.json`
- 多比较火山图: `bulk_reference_output_codex_20260703/dataset_a_core/plots/bulk_deg_volcano_multi.json`
- DEG 热图: `bulk_reference_output_codex_20260703/dataset_a_core/plots/bulk_heatmap.json`
- 最终 h5ad: `bulk_reference_output_codex_20260703/dataset_a_core/intermediate/bulk_heatmap_output.h5ad`

## Dataset B 结果

项目目录:

```text
bulk_reference_output_codex_20260703/dataset_b_core/
```

报告:
- JSON: `bulk_reference_output_codex_20260703/dataset_b_core/bulk_reference_report.json`
- 中文报告: `bulk_reference_output_codex_20260703/dataset_b_core/中文运行报告.md`

模块状态:

| 模块 | 状态 |
|---|---|
| bulk_qc | completed |
| bulk_normalize | completed |
| bulk_pca | completed |
| bulk_deg | completed |
| bulk_heatmap | completed |

差异比较:

| 比较 | 上调基因数 | 下调基因数 | Top 5 基因 |
|---|---:|---:|---|
| NH4Cl-vs-Ctr | 0 | 0 | CTPS1, EIF3C, DDIT4, GARS1, DNAJA1 |
| PEA-vs-Ctr | 0 | 0 | DDIT4, ENSG00000287856, HK2, EGLN1, CLASP1 |
| TMAO-vs-Ctr | 0 | 0 | UPF1, TULP4, FSD2, SPRY2, DDIT4 |

主要输出:
- DEG 总表: `bulk_reference_output_codex_20260703/dataset_b_core/results/bulk_deg_all_comparisons.csv`
- PCA 图: `bulk_reference_output_codex_20260703/dataset_b_core/plots/bulk_pca.json`
- 多比较火山图: `bulk_reference_output_codex_20260703/dataset_b_core/plots/bulk_deg_volcano_multi.json`
- DEG 热图: `bulk_reference_output_codex_20260703/dataset_b_core/plots/bulk_heatmap.json`
- 最终 h5ad: `bulk_reference_output_codex_20260703/dataset_b_core/intermediate/bulk_heatmap_output.h5ad`

## 产物校验

- Plotly JSON: 59 个，全部可解析
- h5ad: 10 个，全部可读取
- 两个中文运行报告均已生成
- 两个 `bulk_deg_all_comparisons.csv` 均已生成
