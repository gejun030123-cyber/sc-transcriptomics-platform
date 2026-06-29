# 集成测试与语义断言测试覆盖补充设计

## 概述

增量补充集成测试和语义断言测试的覆盖范围，不改动已有通过的测试，在现有 4 个测试文件中追加新测试类。

## 约束

- **环境**：仅使用已安装的库（scanpy 1.11, omicverse 2.1, plotly 6.7, scipy, sklearn, inmoose）
- **omicverse 绕过**：QC 用 batch_key='' + 较大合成数据；Normalize 用原始 count 数据
- **不可测模块**（缺少依赖）：annotation, trajectory, cell_communication, bulk_enrichment, bulk_timecourse, bulk_deg_integration, convert_10x — 仅做静态语义扫描
- **语言**：测试注释用中文，变量/函数名用英文（遵循 CLAUDE.md）

## 覆盖矩阵

### 可测模块（集成 + 语义）

| 模块 | 深度集成断言 | 运行时语义 |
|------|-------------|-----------|
| qc | cells_before/after 一致性、cell_cycle_available、输出 QC 列 | summary JSON 可序列化 |
| normalize | n_cells/n_genes 与输出一致、名称保持 | summary JSON 可序列化 |
| hvg | n_hvgs == highly_variable.sum()、force_include 生效 | |
| dimred | pca_variance_ratio_top5、UMAP 2D | plotly JSON 可解析 |
| clustering | 多分辨率 n_clusters_X.X、louvain | JSON 可序列化 |
| qc_reassess | n_clusters/n_low_quality、CSV 列结构 | |
| deg | CSV 列(gene/logfc/pval/padj/cluster)、完整结果 | groups 一致性、CSV 行数 |
| proportion | chi2/p_value/n_batches/n_groups | CSV 行数 |
| batch_correct | summary 字段 | |
| bulk_qc | summary 字段 | |
| bulk_normalize | summary 字段 | |
| bulk_pca | n_components/pc1_variance_pct、X_pca + variance_ratio | |
| bulk_deg | comparison/per_comparison | n_comparisons 一致性 |
| bulk_heatmap | 需 DEG 前置 | |

### 仅静态扫描的模块

annotation, trajectory, cell_communication, bulk_enrichment, bulk_timecourse, bulk_deg_integration, convert_10x

## 各文件新增内容

### test_integration_full.py（新增 ~500 行）

#### 深度模块断言（8 个 SC 类 + 2 个 Bulk 类）

**TestSCModuleQCDeepAssertions**（3 个方法）：
- `test_qc_filter_reduces_cells` — cells_before == cells_after + cells_removed
- `test_qc_summary_has_cell_cycle_info` — cell_cycle_available/s_genes_found/g2m_genes_found 存在
- `test_qc_output_adata_has_qc_columns` — 输出含 total_counts/n_genes_by_counts/pct_counts_mt/novelty_score

**TestSCModuleNormalizeDeepAssertions**（2 个方法）：
- `test_normalize_summary_matches_output` — summary n_cells/n_genes 与 h5ad 一致
- `test_normalize_preserves_cell_names` — obs_names/var_names 不变

**TestSCModuleHVGDeepAssertions**（2 个方法）：
- `test_hvg_summary_n_hvgs_matches_actual` — n_hvgs == highly_variable.sum()
- `test_hvg_force_include_genes` — 强制基因标记为 highly_variable

**TestSCModuleDimredDeepAssertions**（2 个方法）：
- `test_dimred_output_has_pca_variance` — n_pcs/pca_variance_ratio_top5/embedding_method
- `test_dimred_umap_coords_2d` — X_umap.shape[1]==2, X_pca.shape[1]==10

**TestSCModuleClusteringDeepAssertions**（2 个方法）：
- `test_clustering_multi_resolution_summary` — n_clusters_0.3/n_clusters_0.8、resolutions 列表
- `test_clustering_louvain_method` — louvain 方法可正常工作

**TestSCModuleDEGDeepAssertions**（4 个方法）：
- `test_deg_csv_columns_present` — CSV 含 gene/logfc/pval/pval_adj/cluster
- `test_deg_summary_groups_match_output` — n_groups == len(groups)
- `test_deg_with_ttest_method` — t-test 方法正常
- `test_deg_full_results_csv` — 完整 DEG 结果 CSV 非空

**TestSCModuleProportionDeepAssertions**（2 个方法）：
- `test_proportion_summary_has_chi2_and_pval` — chi2/p_value/n_batches/n_groups
- `test_proportion_csv_files_present` — cell_counts 和 cell_proportions CSV

**TestSCModuleQCReassessDeepAssertions**（2 个方法）：
- `test_qc_reassess_summary_has_cluster_stats` — n_clusters/n_low_quality/low_quality_clusters
- `test_qc_reassess_csv_has_cluster_columns` — CSV 含 cluster/n_cells/low_quality/low_reasons

**TestBulkModulePCADeepAssertions**（2 个方法）：
- `test_bulk_pca_summary_fields` — n_components/pc1_variance_pct/pc2_variance_pct/dimred_method
- `test_bulk_pca_output_has_pca_obsm` — X_pca.shape[1]==n_comps, uns['pca']['variance_ratio']

**TestBulkModuleDEGDeepAssertions**（2 个方法）：
- `test_bulk_deg_summary_has_comparison_info` — method/comparison 信息
- `test_bulk_deg_output_has_deg_columns` — CSV 含 pval/padj/logfc 列

#### 端到端全流程（1 个类 2 个方法）

**TestEndToEndFullPipeline**：
- `test_sc_pipeline_qc_to_deg` — QC→Normalize→HVG→Dimred→Clustering→DEG，验证最终输出
- `test_bulk_pipeline_qc_to_deg` — Bulk QC→Normalize→PCA→DEG

#### 跨模块数据完整性（1 个类 3 个方法）

**TestCrossModuleDataIntegrity**：
- `test_intermediate_files_are_valid_h5ad` — 每步 output 可读
- `test_cell_names_preserved_through_pipeline` — QC→Norm→HVG 细胞名一致
- `test_result_files_are_accessible` — 所有 result_files 文件存在非空

#### Progress 回调（1 个类 3 个方法）

**TestProgressCallbackTracking**：
- `test_qc_progress_starts_near_zero_ends_at_100`
- `test_normalize_progress_starts_near_zero_ends_at_100`
- `test_deg_progress_starts_near_zero_ends_at_100`

#### Worker 集成（1 个类 2 个方法）

**TestWorkerIntegration**：
- `test_module_registry_has_all_expected_modules` — 21 个模块全注册
- `test_all_modules_can_be_instantiated` — 所有模块可实例化

### test_semantic.py（新增 ~400 行）

**TestResultFilesConsistency**（2 个方法）：
- `test_clustering_result_files_valid_structure` — 每个 result_file 有完整结构
- `test_normalize_result_files_valid_structure` — 每个 result_file 有完整结构

**TestSummaryJSONSerializable**（4 个方法）：
- `test_qc_summary_json_serializable` — 递归检查无 NaN/Inf/numpy 类型
- `test_normalize_summary_json_serializable` — 同上
- `test_clustering_summary_json_serializable` — 同上
- `test_deg_summary_json_serializable` — 同上

**TestSCDEGSummaryConsistency**（2 个方法）：
- `test_deg_groups_count_matches_list` — n_groups == len(groups)
- `test_deg_csv_rows_consistent_with_summary` — CSV 行数 == total_deg_genes

**TestBulkDEGSummaryConsistency**（2 个方法）：
- `test_bulk_deg_n_comparisons_consistent` — n_comparisons == len(comparisons)
- `test_bulk_deg_per_comparison_has_up_down` — per_comparison 有 n_up/n_down

**TestOutputFileContent**（2 个方法）：
- `test_plotly_json_files_are_parseable` — Plotly JSON 含 data/layout
- `test_csv_files_have_header_and_rows` — CSV 有表头和数据行

**TestInputValidation**（2 个方法）：
- `test_clustering_validate_input_rejects_without_umap` — 无 X_umap 时返回错误
- `test_clustering_validate_input_passes_with_umap` — 有 X_umap 时返回 None

### test_semantic_full.py（新增 ~300 行）

**TestModuleErrorHandling**（2 个方法）：
- `test_no_bare_except_pass` — except:pass 检测（import 保护除外），warning 级别
- `test_error_returns_have_error_key` — 错误返回有 'error' 键，warning 级别

**TestModuleProgressCalls**（3 个方法）：
- `test_all_modules_call_progress` — 所有模块调用 self.progress()
- `test_progress_starts_early_ends_at_100` — 首次 ≤10% 结束 100%
- `test_progress_monotonically_increasing` — 单调递增（-1 除外），warning 级别

**TestModuleParamsAccess**（2 个方法）：
- `test_no_bare_params_subscript` — 无裸 self.params[]，warning 级别
- `test_params_get_has_default_for_required_types` — int/float(params.get()) 有默认值

**TestModuleImports**（1 个方法）：
- `test_heavy_deps_imported_inside_run` — scanpy/omicverse/plotly 在 run() 内导入，warning 级别

**TestResultFilePathSafety**（2 个方法）：
- `test_save_plotly_json_uses_safe_filenames` — 文件名仅含安全字符
- `test_csv_filenames_are_unique_per_module` — 同模块内 CSV 文件名唯一，warning 级别

**TestModuleDisplayInfo**（3 个方法）：
- `test_all_modules_have_display_name` — DISPLAY_NAME 存在且非空
- `test_all_modules_have_description` — DESCRIPTION 存在且非空
- `test_module_names_unique` — MODULE_NAME 唯一

**TestSaveOutputConsistency**（1 个方法）：
- `test_save_output_module_name_matches_declaration` — save_output() 参数与 MODULE_NAME 一致

### test_integration.py（新增 ~200 行）

**TestWorkerModuleInterface**（2 个方法）：
- `test_worker_can_parse_module_result` — Clustering 结果可 JSON 序列化
- `test_worker_result_files_file_paths_exist` — Dimred result_files 文件存在非空

**TestFilterPipeline**（3 个方法）：
- `test_apply_filters_removes_cells` — score >= 0.5 正确过滤
- `test_apply_filters_with_missing_column` — 缺失列不崩溃
- `test_apply_filters_between_operator` — between [5,15] 正确过滤

**TestSchemasIntegration**（2 个方法）：
- `test_all_registered_modules_have_schemas` — 所有注册模块有 PARAM_SCHEMAS（convert_10x 除外）
- `test_schema_metadata_matches_registry` — MODULE_DISPLAY_MAP 与 DISPLAY_NAME 一致或相关

## omicverse 绕过策略

| 模块 | 问题 | 绕过方式 |
|------|------|---------|
| qc | scrublet 对小数据 batch 处理 IndexError | 用 batch_key='' 禁用批次处理，或用无 batch 列的数据 |
| normalize | ov.pp.preprocess mode='shiftlog' IndexError | 用原始 count 数据（未经 normalize 的 AnnData） |

## 测试数据构造

所有测试复用现有的 `_make_sc_anndata()` 和 `_make_bulk_tsv()` helper，不新增 conftest.py。

## 预期结果

- 新增约 1400 行测试代码
- 新增约 50 个测试方法
- 集成测试覆盖：12/21 模块有运行时深度断言（原 8 个仅验证结构）
- 语义断言覆盖：14/21 模块有运行时语义验证 + 21/21 模块有静态扫描
- 所有新增测试应通过（warning 级别的仅记录不 fail）
