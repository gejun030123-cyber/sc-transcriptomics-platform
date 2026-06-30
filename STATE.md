# Loop State — 代码改进循环

## 当前状态
- **循环类型**: 补测试 + 修 bug + 代码质量
- **当前阶段**: 迭代 13 完成 — 全模块测试覆盖 + 代码质量修复
- **上次更新**: 2026-06-29
- **总进度**: 21/21 模块有静态扫描 | 12/21 模块有集成测试 | 14/21 模块有运行时语义验证
- **测试结果**: 131 passed, 9 skipped, 0 failed (15 个测试文件)

## 测试覆盖进度

### 单元测试（工具函数 + 辅助函数）
| 测试文件 | 覆盖模块 | 测试数 | 状态 |
|---------|---------|--------|------|
| test_base.py | base.py (apply_filters, load_adata, save_output, get_plotly_layout) | 18 | ✅ |
| test_io_utils.py | io_utils.py (remap_var_names, read_expression_matrix) | 12 | ✅ |
| test_models.py | models.py (gen_id, Project, AnalysisTask, ResultFile) | 12 | ✅ |
| test_pipeline.py | __init__.py (MODULE_REGISTRY, PIPELINE_ORDER, validate_pipeline_order) | 13 | ✅ |
| test_schemas.py | schemas.py (parse_form_params, PARAM_SCHEMAS) | 17 | ✅ |
| test_visualization.py | visualization.py (transform_heatmap_data, cluster_heatmap, compute_gene_variability) | 15 | ✅ |
| test_heatmap_helpers.py | visualization.py (build_annotation_bar, save_plotly_json) | 10 | ✅ |
| test_expression_parser.py | expression_parser.py (tokenize, parse, evaluate, validate) | 47 | ✅ |
| test_bulk_deg.py | bulk_deg.py (_parse_comparisons, _parse_custom_groups) | 18 | ✅ |
| test_bulk_qc_helpers.py | bulk_qc.py (_gini, _infer_groups, _detect_outliers_mahal) | 10 | ✅ |
| test_normalize_helpers.py | bulk_normalize.py (_tmm_normalize, _vst_transform, _rlog_transform) | 8 | ✅ |
| test_proportion.py | proportion.py (_run_stat_test) | 9 | ✅ |
| test_p2_modules.py | qc/hvg/dimred/batch_correct/clustering/annotation/deg/trajectory/constants | 30 | ✅ |
| test_p3_modules.py | bulk_timecourse/bulk_deg_integration/bulk_heatmap | 32 | ✅ |

### 集成测试（模块端到端 + 深度断言）
| 测试类 | 覆盖模块 | 测试数 | 状态 |
|--------|---------|--------|------|
| TestSCModuleQC | qc | 1 | ✅ |
| TestSCModuleQCDeepAssertions | qc (cells_before/after, cell_cycle, QC columns) | 3 | ✅ |
| TestSCModuleNormalize | normalize | 1 | ✅ |
| TestSCModuleNormalizeDeepAssertions | normalize (summary 匹配, 名称保持) | 2 | ✅ |
| TestSCModuleHVG | hvg | 1 | ✅ |
| TestSCModuleHVGDeepAssertions | hvg (n_hvgs 匹配, force_include) | 2 | ✅ |
| TestSCModuleDimred | dimred | 1 | ✅ |
| TestSCModuleDimredDeepAssertions | dimred (PCA variance, UMAP 2D) | 2 | ✅ |
| TestSCModuleClustering | clustering | 1 | ✅ |
| TestSCModuleClusteringDeepAssertions | clustering (多分辨率, louvain) | 2 | ✅ |
| TestSCModuleQCReassess | qc_reassess | 1 | ✅ |
| TestSCModuleQCReassessDeepAssertions | qc_reassess (簇统计, CSV 列) | 2 | ✅ |
| TestSCModuleDEG | deg | 1 | ✅ |
| TestSCModuleDEGDeepAssertions | deg (CSV 列, groups, t-test, 完整结果) | 4 | ✅ |
| TestSCModuleProportion | proportion | 1 | ✅ |
| TestSCModuleProportionDeepAssertions | proportion (chi2/pval, CSV) | 2 | ✅ |
| TestBulkModuleQC/Normalize/PCA/DEG | bulk 模块 | 4 | ✅ |
| TestBulkModulePCADeepAssertions | bulk_pca (summary, obsm) | 2 | ✅ |
| TestBulkModuleDEGDeepAssertions | bulk_deg (comparison, CSV 列) | 2 | ✅ |
| TestEndToEndPipeline | SC QC→Clustering, Bulk QC→PCA | 2 | ✅ |
| TestEndToEndFullPipeline | SC QC→DEG, Bulk QC→DEG | 2 | ✅ |
| TestCrossModuleDataIntegrity | h5ad 合法性, 名称保持, result_files | 3 | ✅ |
| TestProgressCallbackTracking | QC/Normalize/DEG progress 0→100 | 3 | ✅ |
| TestWorkerIntegration | MODULE_REGISTRY 完整性, 实例化 | 2 | ✅ |

### 跨层集成测试
| 测试类 | 覆盖问题 | 测试数 | 状态 |
|--------|---------|--------|------|
| TestWorkerValidation | worker 预验证 + 非 h5ad 文件 | 3 | ✅ |
| TestPlotlyJsonSanitization | API NaN/Inf 清洗 | 3 | ✅ |
| TestModuleOutputCompleteness | CSV 非空 + 火山图一致性 | 2 | ✅ |
| TestWorkerModuleInterface | 结果 JSON 序列化 + 文件存在 | 2 | ✅ |
| TestFilterPipeline | apply_filters 移除/缺失列/between | 3 | ✅ |
| TestSchemasIntegration | PARAM_SCHEMAS 覆盖 + 显示名一致 | 2 | ✅ |

### 语义断言测试（运行时验证）
| 测试类 | 覆盖问题 | 测试数 | 状态 |
|--------|---------|--------|------|
| TestOutputRegistration | 孤儿 .to_csv() 未注册 | 2 | ✅ |
| TestSummaryConsistency | summary 字段自相矛盾 | 4 | ✅ |
| TestLabelAccuracy | 标签/注释不准确 | 2 | ✅ |
| TestCrossOutputConsistency | 火山图 vs CSV 数据不一致 | 2 | ✅ |
| TestResultFilesConsistency | result_file 结构完整 | 2 | ✅ |
| TestSummaryJSONSerializable | summary 可 JSON 序列化 | 4 | ✅ |
| TestSCDEGSummaryConsistency | n_groups/groups 匹配, CSV 行数 | 2 | ✅ |
| TestBulkDEGSummaryConsistency | n_comparisons, per_comparison | 2 | ✅ |
| TestOutputFileContent | Plotly JSON 可解析, CSV 有表头 | 2 | ✅ |
| TestInputValidation | validate_input 运行时行为 | 2 | ✅ |

### 静态源码扫描（全 21 模块）
| 测试类 | 扫描内容 | 测试数 | 状态 |
|--------|---------|--------|------|
| TestAllModulesSummaryKeys | n_up/n_down 配对, int() 包裹, summary 存在 | 4 | ✅ |
| TestAllModulesResultFiles | 必需键, file_type/category 合法, save_plotly_json | 6 | ✅ |
| TestAllModulesOutputAdata | output_adata 存在, _output.h5ad 约定, save_output | 5 | ✅ |
| TestAllModulesLabelConsistency | Top N 标签, 完整标签, CSV 标签 | 3 | ✅ |
| TestAllModulesInputRequires | INPUT_REQUIRES 匹配 validate_input | 3 | ✅ |
| TestModuleReturnStructure | 返回键完整, result_files 是 list, error 返回 | 5 | ✅ |
| TestAllModulesNamingConventions | snake_case, MODULE_REGISTRY 注册 | 3 | ✅ |
| TestModuleErrorHandling | except:pass 检测, error 返回有 error 键 | 2 | ✅ |
| TestModuleProgressCalls | 调用 progress, 0-100 范围, 单调递增 | 3 | ✅ |
| TestModuleParamsAccess | 无裸 params[], int/float 有默认值 | 2 | ✅ |
| TestModuleImports | 重型依赖在 run() 内导入 | 1 | ✅ |
| TestResultFilePathSafety | 文件名安全, CSV 文件名唯一 | 2 | ✅ |
| TestModuleDisplayInfo | DISPLAY_NAME/DESCRIPTION 存在, MODULE_NAME 唯一 | 3 | ✅ |
| TestSaveOutputConsistency | save_output 参数与 MODULE_NAME 一致 | 1 | ✅ |

## 模块 Bug 修复记录

| # | 文件 | 问题 | 严重度 | 根因 | 发现方式 | 状态 |
|---|------|------|--------|------|---------|------|
| 1 | io_utils.py:32-57 | CSV 被 TSV 解析路径拦截 | ⚠️ 中 | TSV 探测返回 0 列 df | 测试 xfail | ✅ 已修 |
| 2 | worker.py:63 | load_adata 对 CSV 抛 OSError 阻塞任务 | 🔴 高 | 预验证只捕获 FileNotFoundError | 用户报错 | ✅ 已修 |
| 3 | api.py:65-77 | Plotly JSON 含 NaN 导致图表不显示 | 🔴 高 | 二进制解码后 NaN 未清洗 | 用户报错 | ✅ 已修 |
| 4 | bulk_deg.py:151-233 | CSV 为空 + 箱线图无数据 | 🔴 高 | CSV 保存在 base_mean_filter 之后 | 用户报错 | ✅ 已修 |
| 5 | proportion.py:144 | 孤儿 CSV 未注册到 result_files | ⚠️ 中 | .to_csv() 后缺少 append | 语义测试 | ✅ 已修 |
| 6 | deg.py:90-92 | 标签"DEG Results Table"实际只含 top N | ⚠️ 低 | n_genes 限制但标签未体现 | 语义测试 | ✅ 已修 |
| 7 | bulk_deg.py:575 | summary['comparisons'] 含已跳过项 | ⚠️ 低 | comparisons 用原始输入 | 语义测试 | ✅ 已修 |
| 8 | normalize.py:36 | ov.pp.preprocess(mode='shiftlog') IndexError | 🔴 高 | omicverse 需 pipe-separated mode | 集成测试 | ✅ 已修 |
| 9 | proportion.py:99 | Plotly Pie 子图类型不兼容 | ⚠️ 中 | make_subplots 默认 xy 类型 | 集成测试 | ✅ 已修 |
| 10 | convert_10x.py:16 | 裸 self.params['mtx_dir'] KeyError | ⚠️ 中 | 无默认值 | 静态扫描 | ✅ 已修 |

## 代码质量修复记录

| 类型 | 修复内容 | 涉及文件 | 状态 |
|------|---------|---------|------|
| except:pass → logging | 6 处添加 logger.warning() | bulk_deg_integration, cell_communication, clustering, hvg, trajectory | ✅ |
| progress 回退 | 3 个模块 progress 调整为单调递增 | annotation, bulk_enrichment, cell_communication | ✅ |
| 显示名不一致 | schemas.py display 与模块 DISPLAY_NAME 对齐 | bulk_qc, bulk_heatmap, bulk_deg_integration | ✅ |
| 模块代码质量 | 使用 base class 方法 (load_adata, save_output, save_plotly_json) | clustering, cell_communication, qc_reassess | ✅ |
| 路由重构 | 内联定义提取到 modules/schemas.py | routes/analysis.py (-390 行) | ✅ |
| 数据库线程安全 | get_db() 单例 → get_conn() 每次新连接 | database.py | ✅ |

## 模块覆盖矩阵

| 模块 | 单元测试 | 集成测试 | 深度断言 | 运行时语义 | 静态扫描 |
|------|---------|---------|---------|-----------|---------|
| qc | 3 | ✅ | 3 | ✅ | ✅ |
| normalize | 0 | ✅ | 2 | ✅ | ✅ |
| hvg | 6 | ✅ | 2 | — | ✅ |
| dimred | 4 | ✅ | 2 | ✅ | ✅ |
| batch_correct | 2 | ⏭️ env | — | — | ✅ |
| clustering | 2 | ✅ | 2 | ✅ | ✅ |
| qc_reassess | 0 | ✅ | 2 | — | ✅ |
| annotation | 3 | ⏭️ dep | — | — | ✅ |
| deg | 4 | ✅ | 4 | ✅ | ✅ |
| trajectory | 2 | ⏭️ dep | — | — | ✅ |
| proportion | 9 | ✅ | 2 | — | ✅ |
| cell_communication | 0 | ⏭️ dep | — | — | ✅ |
| convert_10x | 0 | — | — | — | ✅ |
| bulk_qc | 10 | ✅ | — | — | ✅ |
| bulk_normalize | 8 | ✅ | — | — | ✅ |
| bulk_pca | 0 | ✅ | 2 | — | ✅ |
| bulk_deg | 18 | ✅ | 2 | ✅ | ✅ |
| bulk_heatmap | 9 | ⏭️ env | — | — | ✅ |
| bulk_enrichment | 0 | ⏭️ dep | — | — | ✅ |
| bulk_timecourse | 17 | ⏭️ env | — | — | ✅ |
| bulk_deg_integration | 6 | ⏭️ env | — | — | ✅ |

**图例**: ✅ 已覆盖 | ⏭️ env 需环境变量 | ⏭️ dep 缺少依赖 | — 未覆盖

## 运行历史

### 迭代 13 (2026-06-29)
- **目标**: 集成测试 + 语义断言测试全覆盖 + 代码质量修复
- **新增测试**: ~50 个（深度断言 + 端到端流水线 + 运行时语义 + 静态扫描）
- **修复 Bug**: 3 个 (#8 normalize omicverse, #9 proportion Plotly, #10 convert_10x params)
- **代码质量**: except:pass→logging, progress 单调性, 显示名对齐, 路由重构, DB 线程安全
- **全量测试**: 131 passed, 9 skipped, 0 failed
- **文件变更**: 38 文件提交 (5442+, 807-)

### 迭代 12 (2026-06-27)
- **目标**: 修 bug #1（最后一个待修 bug）
- **修复**: io_utils.py TSV/自动探测失败时重置 df=None
- **全量测试**: 289 passed, 0 xfailed, 1 warning

### 迭代 11 (2026-06-27)
- **目标**: P2/P3 模块全覆盖
- **测试文件**: test_p2_modules.py (30 个) + test_p3_modules.py (32 个)
- **全量测试**: 286 passed, 3 xfailed, 1 warning

### 迭代 10 (2026-06-27)
- **目标**: 修 bug #5, #6, #7 + 补 P0 测试
- **修复**: proportion.py 孤儿 CSV, deg.py 标签误导, bulk_deg.py comparisons 不一致
- **全量测试**: 224 passed, 3 xfailed, 1 warning

### 迭代 9 (2026-06-27)
- **目标**: 补充集成测试 + 语义断言测试
- **测试文件**: test_integration.py + test_semantic.py (18 个)
- **发现 Bug**: 3 个 (worker OSError, API NaN, CSV 空)
- **全量测试**: 218 passed, 3 xfailed
