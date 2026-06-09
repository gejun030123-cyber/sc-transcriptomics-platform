import os
import json
from flask import Blueprint, render_template, request, redirect, url_for, flash
from models import Project, AnalysisTask
from config import Config

analysis_bp = Blueprint('analysis', __name__)

SC_MODULE_LIST = [
    {'name': 'qc', 'display': '质控', 'desc': 'MT/ribo/hb 过滤 + 双细胞检测 + 细胞周期评分 + 复杂度过滤'},
    {'name': 'preprocess', 'display': '预处理', 'desc': '标准化，选择高变异基因'},
    {'name': 'dimred', 'display': '降维分析', 'desc': 'PCA, UMAP'},
    {'name': 'batch_correct', 'display': '批次校正', 'desc': 'Harmony, ComBat, SysVI'},
    {'name': 'clustering', 'display': '聚类分析', 'desc': 'Leiden 聚类'},
    {'name': 'qc_reassess', 'display': 'QC 重新评估', 'desc': '聚类后检查 doublet 和 QC 指标，标记低质量簇'},
    {'name': 'annotation', 'display': '细胞注释', 'desc': '基于 Marker 的细胞类型注释'},
    {'name': 'deg', 'display': '差异表达', 'desc': '差异表达基因分析'},
    {'name': 'trajectory', 'display': '轨迹分析', 'desc': '拟时序分析'},
    {'name': 'proportion', 'display': '比例分析', 'desc': '细胞比例分析'},
]

BULK_MODULE_LIST = [
    {'name': 'bulk_qc', 'display': '数据质控', 'desc': '文库大小、基因检测、离群值过滤'},
    {'name': 'bulk_normalize', 'display': '数据标准化', 'desc': 'DESeq2 / CPM / 分位数标准化'},
    {'name': 'bulk_deg', 'display': '差异表达分析', 'desc': '组间差异基因检测（火山图、MA图）'},
    {'name': 'bulk_pca', 'display': 'PCA / UMAP', 'desc': 'PCA 和 UMAP 降维可视化'},
    {'name': 'bulk_heatmap', 'display': '热图分析', 'desc': 'Top 差异基因热图、样本相关性热图'},
    {'name': 'bulk_enrichment', 'display': '通路富集', 'desc': 'GO/KEGG/WikiPathways 通路富集分析（ORA / GSEA）'},
    {'name': 'bulk_timecourse', 'display': '时序分析', 'desc': '多时间点差异基因检测 + 轨迹聚类'},
    {'name': 'bulk_deg_integration', 'display': '多组差异整合', 'desc': '多组比较结果整合：Upset 图、一致性评分、logFC 矩阵分析'},
]

MODULE_LIST = SC_MODULE_LIST + BULK_MODULE_LIST

MODULE_DISPLAY_MAP = {m['name']: m['display'] for m in MODULE_LIST}

SC_MODULE_NAMES = {m['name'] for m in SC_MODULE_LIST}
BULK_MODULE_NAMES = {m['name'] for m in BULK_MODULE_LIST}

STATUS_MAP = {
    'pending': '等待中',
    'running': '运行中',
    'completed': '已完成',
    'failed': '失败',
    'processing': '处理中',
}

PARAM_SCHEMAS = {
    'qc': [
        {'key': 'mito_perc', 'label': '最大线粒体比例', 'type': 'number', 'default': 0.2, 'step': 0.01, 'help': '过滤线粒体基因比例高于此阈值的细胞。人类样本建议 0.1-0.2，小鼠可放宽至 0.25。过高保留低质量细胞，过低丢失应激细胞。'},
        {'key': 'nUMIs', 'label': '最小 UMI 数', 'type': 'number', 'default': 500, 'help': '每个细胞的最小 UMI 总数。低于此值的细胞被视为碎片或死细胞。常用范围 500-1000。'},
        {'key': 'detected_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 250, 'help': '每个细胞检测到的最小基因数。低于此值的细胞可能为低质量或空液滴。常用范围 200-500。'},
        {'key': 'max_detected_genes', 'label': '最大检测基因数（0 = 不限制）', 'type': 'number', 'default': 0, 'help': '每个细胞检测到的最大基因数。高于此值的细胞可能是双细胞或聚合物。设为 0 表示不限制。建议范围 5000-8000，根据数据分布调整。'},
        {'key': 'ribo_perc', 'label': '最大核糖体比例 %（0 = 不过滤）', 'type': 'number', 'default': 0, 'step': 1.0, 'help': '过滤核糖体蛋白基因比例高于此值的细胞。核糖体比例过高可能反映细胞应激或人为扩增。设为 0 表示不过滤。建议范围 30-50%。'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch', 'help': 'adata.obs 中标识实验批次的列名。用于分批次运行 Scrublet 双细胞检测。'},
    ],
    'preprocess': [
        {'key': 'n_top_genes', 'label': '高变异基因数量', 'type': 'number', 'default': 2000, 'help': '选择的高变异基因数量。2000 为标准值，适合大多数分析。基因数过少会丢失生物学信号，过多会引入噪声。'},
        {'key': 'target_sum', 'label': '标准化目标总数', 'type': 'number', 'default': 10000, 'help': '每个细胞标准化后的总计数目标。10000 为 scanpy 默认值，设为 None 则中位数标准化。'},
        {'key': 'batch_key', 'label': '批次列名（可选）', 'type': 'text', 'default': '', 'help': '用于 HVG 选择的批次校正。留空则不进行批次校正。设置后会在选择高变异基因时考虑批次效应。'},
    ],
    'dimred': [
        {'key': 'n_comps', 'label': 'PCA 主成分数量', 'type': 'number', 'default': 50, 'help': 'PCA 主成分数量。通常 30-50 即可捕获大部分方差。过多会引入噪声维度。可通过肘部图选择。'},
        {'key': 'use_mde', 'label': '使用 MDE（加速 UMAP）', 'type': 'checkbox', 'default': False, 'help': '使用 Minimum Dystortion Embedding 替代标准 UMAP，速度更快但结果略有差异。'},
    ],
    'batch_correct': [
        {'key': 'method', 'label': '校正方法', 'type': 'select', 'options': ['harmony', 'combat', 'sysvi'], 'default': 'harmony', 'help': '批次校正方法。Harmony：速度快，推荐首选。ComBat：适用于已知批次的情况。SysVI：深度学习方法，适合复杂批次效应，需要 GPU。'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch', 'help': 'adata.obs 中标识批次的列名。必须存在且包含至少 2 个不同的批次值。'},
        {'key': 'n_pcs', 'label': '主成分数量', 'type': 'number', 'default': 50, 'help': '用于批次校正的主成分数量。通常与降维步骤使用的 n_comps 一致。'},
        {'key': 'max_epochs', 'label': '最大迭代轮数（仅 SysVI）', 'type': 'number', 'default': 200, 'help': 'SysVI 方法的最大训练轮数。仅对 SysVI 方法有效。200 为默认值，复杂数据可能需要更多。'},
    ],
    'clustering': [
        {'key': 'resolutions', 'label': 'Leiden 分辨率（逗号分隔）', 'type': 'text', 'default': '0.6,0.8,1.0', 'help': 'Leiden 聚类分辨率，多个值用逗号分隔。值越大聚类越细（簇越多）。0.4-0.6 适合粗分，0.8-1.0 标准，>1.0 细分。建议测试多个值。'},
        {'key': 'n_neighbors', 'label': '邻居数量', 'type': 'number', 'default': 15, 'help': '构建 KNN 图时的邻居数量。值越大聚类越平滑，越小越敏感。15 为默认值，小数据集可降至 10。'},
    ],
    'qc_reassess': [
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden',
         'help': '用于分组的聚类列名。通常为 leiden 或 leiden_0.8 等。'},
        {'key': 'doublet_threshold', 'label': 'Doublet 比例阈值', 'type': 'number', 'default': 0.3, 'step': 0.05,
         'help': 'doublet 比例高于此值的簇被标记为低质量。0.3 表示 30%。'},
        {'key': 'mt_threshold', 'label': 'MT 比例阈值', 'type': 'number', 'default': 15.0, 'step': 1.0,
         'help': '平均线粒体基因比例高于此值的簇被标记为低质量。'},
    ],
    'annotation': [
        {'key': 'method', 'label': '注释方法', 'type': 'select', 'options': ['auto_marker', 'manual'], 'default': 'auto_marker', 'help': '注释方法。auto_marker：使用内置 TME marker 基因自动打分。manual：手动指定 marker 基因。'},
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden', 'help': '用于分组的聚类列名。通常为 leiden 或 leiden_0.8 等。'},
        {'key': 'resolution', 'label': 'Leiden 分辨率', 'type': 'text', 'default': '0.8', 'help': '对应的 Leiden 分辨率，用于定位正确的聚类列。'},
        {'key': 'marker_set', 'label': 'Marker 基因集', 'type': 'select', 'options': ['TME', 'Immune', 'Blood'], 'default': 'TME',
         'help': 'TME：肿瘤微环境 marker（上皮、CAF、内皮、免疫细胞等）。Immune：免疫细胞 marker（T、B、NK、髓系等）。Blood：血液细胞 marker（HSC、红系、巨核、粒系等）。'},
        {'key': 'custom_markers', 'label': '自定义 Marker（可选）', 'type': 'textarea', 'default': '',
         'help': '每行一个细胞类型，格式：CellType:GENE1,GENE2。示例：\nT_cell:CD3D,CD3E,CD2\nB_cell:CD19,MS4A1,CD79A\nMacrophage:CD68,CD163,MSR1'},
    ],
    'deg': [
        {'key': 'groupby', 'label': '分组依据', 'type': 'text', 'default': '', 'help': '差异分析的分组依据列名。如 celltype、leiden、condition 等。留空则自动使用 leiden。'},
        {'key': 'reference', 'label': '参考组', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': 'rest', 'help': '参考组。rest 表示以所有其他组为对照。选择特定组则以该组为对照。'},
        {'key': 'method', 'label': '统计方法', 'type': 'select', 'options': ['wilcoxon', 't-test', 'logreg'], 'default': 'wilcoxon', 'help': '统计检验方法。Wilcoxon：非参数检验，最常用。t-test：参数检验。logreg：逻辑回归。'},
        {'key': 'n_genes', 'label': '显示 Top N 基因数', 'type': 'number', 'default': 20, 'help': '每个簇显示的 Top N 差异基因数。20 为标准值。'},
        {'key': 'show_dotplot', 'label': '生成 DEG Dotplot', 'type': 'checkbox', 'default': True,
         'help': '是否生成 Top 差异基因的 dotplot 可视化。'},
        {'key': 'plot_genes_umap', 'label': 'UMAP 展示基因（逗号分隔）', 'type': 'text', 'default': '',
         'help': '指定要在 UMAP 上展示表达分布的基因名，多个用逗号分隔。留空则不生成。'},
        {'key': 'custom_dotplot_genes', 'label': '自定义 Dotplot 基因（可选）', 'type': 'textarea', 'default': '',
         'help': '手动输入基因名，逗号或换行分隔。填写后 Dotplot 使用此列表而非自动 Top N DEG。'},
    ],
    'trajectory': [
        {'key': 'method', 'label': '轨迹方法', 'type': 'select', 'options': ['diffusion_map', 'slingshot'], 'default': 'diffusion_map', 'help': '轨迹推断方法。diffusion_map：基于扩散图的拟时序，适合连续过渡。slingshot：基于 MST 的轨迹，适合分支结构。'},
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden', 'help': '用于轨迹推断的聚类列名。'},
        {'key': 'plot_genes', 'label': '拟时序基因表达（可选）', 'type': 'textarea', 'default': '',
         'help': '手动输入基因名，逗号或换行分隔。生成这些基因沿拟时序的表达曲线图。最多 10 个基因。'},
    ],
    'proportion': [
        {'key': 'groupby', 'label': '分组依据', 'type': 'text', 'default': 'celltype', 'help': '统计比例的细胞类型列名。通常为 celltype 或 leiden。'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch', 'help': '用于比较的分组列名。如 batch、condition、treatment 等。'},
        {'key': 'compare_groups', 'label': '指定比较组（可选）', 'type': 'text', 'default': '',
         'help': '格式：GroupA-vs-GroupB，多个比较用分号分隔（如 A-vs-B;C-vs-D）。仅比较指定组的细胞比例差异。留空则比较所有组。'},
    ],
    'bulk_qc': [
        {'key': 'min_counts', 'label': '最小文库 reads 数', 'type': 'number', 'default': 100000, 'help': '最小文库 reads 数。低于此值的样本被过滤。人类/小鼠 RNA-seq 通常要求 ≥100000，小样本可降至 50000。'},
        {'key': 'min_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 5000, 'help': '每个样本检测到的最小基因数。低于此值的样本可能质量差。通常 5000-8000。'},
        {'key': 'max_mt_pct', 'label': '最大线粒体基因比例 (%)', 'type': 'number', 'default': 20.0, 'step': 0.1, 'help': '最大线粒体基因比例（%）。高于此值的样本可能降解严重。RNA-seq 通常 15-20%。'},
        {'key': 'max_ribo_pct', 'label': '最大核糖体基因比例 (%)', 'type': 'number', 'default': 40.0, 'step': 0.1, 'help': '最大核糖体基因比例（%）。RPL/RPS 基因比例过高提示 rRNA 污染。PolyA 建库通常 < 5-10%，rRNA 去除建库可至 40-50%。'},
        {'key': 'min_gini', 'label': '最小文库复杂度 (Gini)', 'type': 'number', 'default': 0, 'step': 0.01, 'help': '最小 Gini 系数（0 = 不过滤）。Gini > 0.8 提示文库复杂度低（PCR 过度扩增）。'},
        {'key': 'min_sample_expr', 'label': '基因最低表达样本数', 'type': 'number', 'default': 0, 'step': 1, 'help': '基因在至少 N 个样本中 CPM > 1 才保留。0 = 不过滤。建议设为样本总数的 10-20%。'},
        {'key': 'group_column', 'label': '分组列名（可选）', 'type': 'text', 'default': '', 'help': '样本分组列名（adata.obs 中的列）。留空则自动从样本名推断（取第一个分隔符前的前缀）。填写后启用组内/组间距离分析和分组着色图。'},
        {'key': 'detect_outliers', 'label': '检测离群样本', 'type': 'checkbox', 'default': True, 'help': '基于 PCA 马氏距离检测离群样本。仅在 summary 中告警，不自动剔除。'},
        {'key': 'filter_strategy', 'label': '过滤策略', 'type': 'select', 'options': ['conservative', 'standard', 'custom'], 'default': 'standard', 'help': 'conservative：宽松阈值（适合小样本）；standard：推荐阈值；custom：自定义所有阈值。'},
    ],
    'bulk_normalize': [
        {'key': 'method', 'label': '标准化方法', 'type': 'select', 'options': ['deseq2', 'tmm', 'cpm', 'vst', 'rlog', 'log2_quantile'], 'default': 'deseq2',
         'help': '标准化方法。差异分析：DESeq2（中位比率法，金标准）或 TMM（edgeR 方法，组成偏差大时更优）。可视化/高维：VST（方差稳定）或 rlog（小样本更稳定）。简单归一：CPM（每百万计数）或 log2 分位数。'},
        {'key': 'min_expr_value', 'label': '最小表达阈值 (CPM)', 'type': 'number', 'default': 1, 'step': 0.1,
         'help': '基因表达量需达到此 CPM 阈值才算有效表达。默认 1。'},
        {'key': 'min_expr_samples', 'label': '最小表达样本数', 'type': 'number', 'default': 3, 'step': 1,
         'help': '基因在至少 N 个样本中达到最小表达阈值才保留。0 = 不过滤。建议设为最小组的样本数。'},
        {'key': 'max_zero_pct', 'label': '最大零值比例 (%)', 'type': 'number', 'default': 0, 'step': 1,
         'help': '基因在超过此比例的样本中为零则被过滤。0 = 不过滤。建议 50-70%。'},
    ],
    'bulk_deg': [
        {'key': 'groupby', 'label': '分组列名', 'type': 'text', 'default': '', 'help': '分组列名。adata.obs 中用于区分实验组和对照组的列。如 condition、treatment、group。'},
        {'key': 'group1', 'label': '实验组', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': '', 'help': '实验组名称。将与对照组比较计算差异基因。'},
        {'key': 'group2', 'label': '对照组', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': '', 'help': '对照组名称。rest 表示以所有其他样本为对照。'},
        {'key': 'method', 'label': '统计方法', 'type': 'select', 'options': ['t-test', 'mann-whitney', 'deseq2', 'edger', 'limma'], 'default': 't-test', 'help': '统计方法。t-test：参数检验，适合正态分布数据，速度快。Mann-Whitney：非参数检验，不假设正态分布，更稳健。DESeq2：基于负二项分布的差异分析，RNA-seq 金标准，需要原始计数。edgeR：基于负二项分布模型和经验贝叶斯方法，适合多组比较和复杂实验设计。limma-voom：基于线性模型和经验贝叶斯收缩，适合复杂实验设计，稳健且灵敏。'},
        {'key': 'fc_threshold', 'label': 'Fold Change 阈值', 'type': 'number', 'default': 2.0, 'step': 0.1, 'help': 'Fold Change 阈值。log2FC > log2(fc) 为上调，< -log2(fc) 为下调。常用值：1.5（宽松）、2.0（标准）、4.0（严格）。'},
        {'key': 'pval_threshold', 'label': 'padj 显著性阈值', 'type': 'number', 'default': 0.05, 'step': 0.01, 'help': '调整后 p-value 显著性阈值。0.05 为标准，0.01 为严格，0.1 为宽松探索性分析。'},
        {'key': 'top_n', 'label': 'Top N 差异基因数', 'type': 'number', 'default': 20, 'help': '结果中展示的 Top N 差异基因数。用于火山图标注和 Top 基因 CSV 导出。'},
        {'key': 'base_mean_filter', 'label': '最低平均表达量', 'type': 'number', 'default': 1, 'step': 0.5, 'help': '过滤低表达基因。BaseMean 低于此值的基因不参与分析和绘图。建议 1-10。'},
        {'key': 'plot_genes', 'label': '额外展示基因（逗号分隔，可选）', 'type': 'text', 'default': '', 'help': '指定要额外绘制箱线图的基因名，多个用逗号分隔。Top 差异基因会自动生成箱线图，此字段用于补充其他感兴趣的基因。'},
        {'key': 'comparisons', 'label': '多组比较（可选）', 'type': 'text', 'default': '',
         'help': '多个比较用分号分隔，格式：A-vs-B;C-vs-D。填写后实验组/对照组参数被忽略。示例：DrugA-vs-Control;DrugB-vs-Control;DrugA-vs-DrugB'},
        {'key': 'custom_groups', 'label': '自定义合并组（可选）', 'type': 'textarea', 'default': '',
         'help': '每行一个定义，格式：新组名=原组1+原组2。示例：High=Treated_1h+Treated_3h。定义后可在多组比较中使用新组名。'},
        {'key': 'test_type', 'label': '检验类型', 'type': 'select', 'options': ['pairwise', 'lrt'], 'default': 'pairwise',
         'help': 'pairwise：两两比较（默认）。lrt：似然比检验（仅 edger），一次性检验所有组间是否有差异。'},
        {'key': 'auto_comparisons', 'label': '自动生成比较', 'type': 'select', 'options': ['manual', 'all_pairwise', 'vs_reference'], 'default': 'manual',
         'help': 'manual：手动输入比较。all_pairwise：自动生成所有两两配对。vs_reference：所有组 vs 参考组。'},
        {'key': 'reference_group', 'label': '参考组（可选）', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': '',
         'help': '指定参考组。确保 logFC 方向一致（正值=该组>参考组）。'},
        {'key': 'cooks_filter', 'label': "Cook's 距离过滤", 'type': 'checkbox', 'default': True,
         'help': "剔除 Cook's 距离过大的异常高表达基因。DESeq2 和 edgeR 支持。"},
        {'key': 'independent_filter', 'label': '独立过滤', 'type': 'checkbox', 'default': True,
         'help': '自动去除低表达基因，提升检测效力。DESeq2 支持。'},
        {'key': 'padj_method', 'label': 'p 值校正方法', 'type': 'select', 'options': ['fdr_bh', 'bonferroni', 'holm', 'fdr_by'], 'default': 'fdr_bh',
         'help': '多重检验校正方法。BH（Benjamini-Hochberg）：最常用。Bonferroni：最严格。Holm：逐步校正。BY：依赖性校正。'},
        {'key': 'regulation_filter', 'label': '差异方向', 'type': 'select', 'options': ['both', 'up', 'down'], 'default': 'both',
         'help': 'both：双向差异基因。up：仅输出上调基因。down：仅输出下调基因。'},
    ],
    'bulk_pca': [
        {'key': 'n_comps', 'label': 'PCA 主成分数量', 'type': 'number', 'default': 10, 'help': 'PCA 主成分数量。通常 5-10 即可。样本数少时自动降至 n_samples-1。'},
        {'key': 'color_by', 'label': '颜色分组列名（留空则不着色）', 'type': 'text', 'default': '', 'help': '用于着色的 obs 列名。留空则不着色。如 condition、batch、celltype 等。'},
        {'key': 'dimred_method', 'label': '降维方法', 'type': 'select', 'options': ['pca', 'umap', 'tsne'], 'default': 'pca', 'help': '降维可视化方法。PCA：线性降维，保留全局结构。UMAP：非线性降维，保留局部结构。t-SNE：非线性降维，适合发现聚类。'},
    ],
    'bulk_heatmap': [
        # === 基因选择 ===
        {'key': 'gene_import_source', 'label': '基因来源', 'type': 'select',
         'options': ['top_var', 'deg', 'expression_filter', 'manual'], 'default': 'top_var',
         'help': '基因来源。top_var：最高变异基因。deg：差异表达基因。expression_filter：表达式筛选。manual：手动输入。'},
        {'key': 'var_metric', 'label': '变异度量', 'type': 'select',
         'options': ['var', 'mad', 'cv', 'range'], 'default': 'var',
         'help': 'top_var 模式下的变异度量。var=方差，mad=中位绝对偏差，cv=变异系数，range=极差。'},
        {'key': 'deg_direction', 'label': 'DEG 方向过滤', 'type': 'select',
         'options': ['both', 'up', 'down'], 'default': 'both',
         'help': '仅显示上调/下调/全部差异基因。'},
        {'key': 'deg_sortby', 'label': 'DEG 排序方式', 'type': 'select',
         'options': ['padj', 'abs_logfc', 'logfc'], 'default': 'padj',
         'help': '差异基因在热图中的排列顺序。'},
        {'key': 'filter_expression', 'label': '表达式筛选', 'type': 'textarea', 'default': '',
         'help': '集合逻辑表达式（如 ALL:up），仅 gene_import_source=expression_filter 时生效。'},
        {'key': 'top_n', 'label': '显示基因数', 'type': 'number', 'default': 50, 'step': 5,
         'help': '热图中显示的基因数量。通常 30-100。'},
        {'key': 'custom_genes', 'label': '自定义基因列表', 'type': 'textarea', 'default': '',
         'help': '手动输入基因名，逗号或换行分隔。gene_import_source=manual 时生效。'},
        # === 数据变换 ===
        {'key': 'row_scaling', 'label': '行标准化', 'type': 'select',
         'options': ['zscore', 'center', 'none'], 'default': 'zscore',
         'help': 'zscore=Z 分数标准化（推荐），center=仅中心化，none=不变换。'},
        {'key': 'pseudocount', 'label': '伪计数', 'type': 'number', 'default': 1, 'step': 0.1,
         'help': 'log2 转换前加值，避免 log(0)。'},
        {'key': 'winsorize', 'label': '极端值截断', 'type': 'select',
         'options': ['none', '1pct', '5pct'], 'default': 'none',
         'help': '对极端值进行百分位截断，避免 outliers 主导颜色。'},
        {'key': 'clip_range', 'label': '颜色截断范围', 'type': 'text', 'default': '-3,3',
         'help': '标准化后的值截断范围，逗号分隔（如 -3,3）。留空不截断。'},
        {'key': 'missing_value', 'label': '缺失值处理', 'type': 'select',
         'options': ['ignore', 'mean_fill', 'zero_fill'], 'default': 'ignore',
         'help': '缺失值处理方式。ignore=忽略，mean_fill=行均值填充，zero_fill=填 0。'},
        # === 聚类 ===
        {'key': 'row_cluster', 'label': '基因聚类', 'type': 'select',
         'options': ['yes', 'no'], 'default': 'yes',
         'help': '是否对基因进行层次聚类。'},
        {'key': 'col_cluster', 'label': '样本聚类', 'type': 'select',
         'options': ['yes', 'no', 'group_order'], 'default': 'yes',
         'help': '样本聚类。group_order=按分组固定顺序排列。'},
        {'key': 'row_method', 'label': '基因聚类方法', 'type': 'select',
         'options': ['ward', 'complete', 'average', 'single', 'mcquitty'], 'default': 'ward',
         'help': '层次聚类算法。ward=最小方差（推荐），complete=最长距离，average=平均距离。'},
        {'key': 'col_method', 'label': '样本聚类方法', 'type': 'select',
         'options': ['ward', 'complete', 'average', 'single', 'mcquitty'], 'default': 'ward',
         'help': '样本层次聚类算法。'},
        {'key': 'row_metric', 'label': '基因距离度量', 'type': 'select',
         'options': ['euclidean', 'pearson', 'spearman', 'cosine'], 'default': 'euclidean',
         'help': '基因间距离计算方式。euclidean=欧氏距离，pearson=相关距离。'},
        {'key': 'col_metric', 'label': '样本距离度量', 'type': 'select',
         'options': ['euclidean', 'pearson', 'spearman', 'cosine'], 'default': 'euclidean',
         'help': '样本间距离计算方式。'},
        # === 视觉样式 ===
        {'key': 'colorscale', 'label': '颜色方案', 'type': 'select',
         'options': ['RdBu_r', 'RdYlBu', 'viridis', 'magma', 'Blues', 'PiYG'], 'default': 'RdBu_r',
         'help': '热图颜色方案。RdBu_r=红蓝发散（推荐），viridis=感知均匀。'},
        {'key': 'reverse_color', 'label': '反转颜色', 'type': 'checkbox', 'default': False,
         'help': '反转颜色映射方向。'},
        {'key': 'zmin', 'label': '最小值', 'type': 'text', 'default': 'auto',
         'help': '颜色映射最小值。auto=自动，或输入数字（如 -3）。'},
        {'key': 'zmax', 'label': '最大值', 'type': 'text', 'default': 'auto',
         'help': '颜色映射最大值。auto=自动，或输入数字（如 3）。'},
        {'key': 'show_gene_labels', 'label': '基因名显示', 'type': 'select',
         'options': ['all', 'top20', 'none'], 'default': 'all',
         'help': '基因名标签显示方式。'},
        {'key': 'show_sample_labels', 'label': '样本名显示', 'type': 'select',
         'options': ['all', 'none'], 'default': 'all',
         'help': '样本名标签显示方式。'},
        {'key': 'gene_font_size', 'label': '基因名字体大小', 'type': 'number', 'default': 8, 'step': 1,
         'help': '基因名标签字体大小。'},
        {'key': 'sample_font_size', 'label': '样本名字体大小', 'type': 'number', 'default': 9, 'step': 1,
         'help': '样本名标签字体大小。'},
        {'key': 'annotation_columns', 'label': '注释条列名', 'type': 'text', 'default': '',
         'help': '额外注释条列名，逗号分隔（如 group,batch）。groupby 列自动包含。'},
        {'key': 'groupby', 'label': '样本分组列名', 'type': 'text', 'default': '',
         'help': '主分组列名，用于默认注释条。留空则不添加。'},
        # === 第二阶段新增 ===
        {'key': 'deg_comparison', 'label': 'DEG 比较选择', 'type': 'text', 'default': '',
         'help': '指定使用哪组 DEG 结果（如 results_0），留空取最新。gene_import_source=deg 时生效。'},
        {'key': 'up_down_separate', 'label': '上下调分开排列', 'type': 'checkbox', 'default': False,
         'help': '上调基因在上、下调在下，中间留间隙。'},
        {'key': 'gene_annotation_columns', 'label': '基因注释列名', 'type': 'text', 'default': '',
         'help': '从 adata.var 中读取基因属性作为注释条，逗号分隔（如 gene_type,pathway）。'},
        {'key': 'log_transform', 'label': 'Log2 转换', 'type': 'select',
         'options': ['auto', 'yes', 'no'], 'default': 'auto',
         'help': 'auto=根据数据是否已标准化自动判断，yes=强制 log2(x+pseudocount)，no=不转换。'},
        {'key': 'corr_method', 'label': '相关性方法', 'type': 'select',
         'options': ['pearson', 'spearman'], 'default': 'pearson',
         'help': '样本相关性热图的计算方法。'},
        {'key': 'corr_colorscale', 'label': '相关性热图配色', 'type': 'select',
         'options': ['Blues', 'RdBu_r', 'viridis', 'YlOrRd'], 'default': 'Blues',
         'help': '相关性热图的颜色方案。'},
        {'key': 'annotation_palette', 'label': '注释条自定义配色', 'type': 'text', 'default': '',
         'help': '格式：ctrl=#1565c0,hmc3=#e53935。留空使用默认配色。'},
    ],
    'bulk_enrichment': [
        {'key': 'method', 'label': '富集方法', 'type': 'select', 'options': ['ORA', 'GSEA'], 'default': 'ORA',
         'help': 'ORA：超几何检验，输入 DE 基因列表，检验哪些通路过度代表。GSEA：秩检验，输入全基因按差异排序，检验通路在排序中的富集位置。'},
        {'key': 'database', 'label': '基因集数据库', 'type': 'select',
         'options': ['GO_BP', 'GO_MF', 'GO_CC', 'KEGG', 'WikiPathways', 'Reactome'], 'default': 'GO_BP',
         'help': 'GO_BP：生物过程。GO_MF：分子功能。GO_CC：细胞组分。KEGG：代谢和信号通路。WikiPathways：社区维护通路。Reactome：反应组数据库。'},
        {'key': 'organism', 'label': '物种', 'type': 'select', 'options': ['Human', 'Mouse'], 'default': 'Human',
         'help': 'Human：人类基因。Mouse：小鼠基因。基因 ID 需与所选物种匹配。'},
        {'key': 'pvalue_cutoff', 'label': '显著性阈值', 'type': 'number', 'default': 0.05, 'step': 0.01,
         'help': '调整后 p-value 截断值。0.05 为标准，0.01 为严格。'},
        {'key': 'top_n', 'label': '展示通路数', 'type': 'number', 'default': 20,
         'help': '可视化中显示的 Top N 显著通路数。'},
        {'key': 'input_source', 'label': 'DEG 结果文件路径', 'type': 'text', 'default': '',
         'help': '来自已完成的 DEG 分析的 CSV 结果文件路径。包含 gene 和 regulation/log2FC 列。'},
        {'key': 'split_direction', 'label': '分开分析上调/下调基因', 'type': 'checkbox', 'default': False, 'help': '开启后将 DEG 结果按 Up/Down 拆分，分别做 ORA 富集分析，生成独立的气泡图。'},
        {'key': 'custom_genes', 'label': '自定义基因列表（可选）', 'type': 'textarea', 'default': '', 'help': '手动输入基因名，逗号或换行分隔。填写后忽略 DEG 结果文件，直接用此列表做 ORA。'},
    ],
    'bulk_timecourse': [
        {'key': 'time_column', 'label': '时间列名', 'type': 'text', 'default': 'minute',
         'help': 'obs 中包含数值型时间点信息的列名，如 minute、hour、day。列中的值必须为数值。'},
        {'key': 'group_column', 'label': '分组列名（可选，用于交互检验）', 'type': 'text', 'default': '',
         'help': '用于分组的 obs 列名。填写后将执行 时间 x 分组 交互效应 F-test。留空则仅运行时序主效应分析。'},
        {'key': 'spline_df', 'label': 'Spline 自由度', 'type': 'number', 'default': 3,
         'help': '自然立方 spline 的自由度。3 为常用默认值，可拟合非线性时间趋势。增大可拟合更复杂的曲线，但消耗更多自由度。'},
        {'key': 'n_clusters', 'label': '聚类数', 'type': 'number', 'default': 6,
         'help': '模糊 c-means 聚类数。通常 4-8 可覆盖主要时间表达模式。需满足：显著时序基因数 >= 聚类数。'},
        {'key': 'fdr_threshold', 'label': 'FDR 阈值', 'type': 'number', 'default': 0.05, 'step': 0.01,
         'help': 'BH 校正后的 FDR 显著性阈值。0.05 为标准，0.01 为严格。'},
        {'key': 'pairwise_groups', 'label': '分组配对比较（可选）', 'type': 'text', 'default': '',
         'help': '格式：GroupA-vs-GroupB。在每个时间点对两组做 Welch t-test，生成时序差异热图。需同时填写分组列名。'},
    ],
    'bulk_deg_integration': [
        {'key': 'selected_comparisons', 'label': '选择比较（留空=全部）', 'type': 'multiselect',
         'api': '/api/projects/{pid}/deg-comparisons', 'default': '',
         'help': '勾选要参与整合分析的比较结果。不勾选则使用全部比较。'},
        {'key': 'min_comparisons', 'label': '最小比较数', 'type': 'number', 'default': 2, 'step': 1,
         'help': '基因至少在 N 个比较中显著才纳入一致性分析。建议 2-3。'},
        {'key': 'consistency_n', 'label': 'Top N 一致性基因', 'type': 'number', 'default': 50, 'step': 5,
         'help': '一致性评分最高的 Top N 基因用于热图和排名展示。'},
        {'key': 'fc_threshold', 'label': 'Fold Change 阈值', 'type': 'number', 'default': 2.0, 'step': 0.1,
         'help': '差异基因判定的 FC 阈值（与 bulk_deg 保持一致）。'},
        {'key': 'pval_threshold', 'label': 'padj 显著性阈值', 'type': 'number', 'default': 0.05, 'step': 0.01,
         'help': '差异基因判定的 padj 阈值（与 bulk_deg 保持一致）。'},
        {'key': 'filter_expression', 'label': '筛选表达式（可选）', 'type': 'textarea', 'default': '',
         'help': '通过集合逻辑表达式筛选目标基因集。支持 AND/OR/NOT/XOR 运算符和 ALL/ANY/ONLY 简写。'
                 '示例：ALL:up | hmc3-vs-ctrl:up AND rapa-vs-ctrl:down | ONLY[hmc3-vs-ctrl]:up'},
    ],
}


@analysis_bp.route('/<pid>/sc-analysis')
def sc_analysis_list(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    return render_template('sc_analysis.html', project=p, modules=SC_MODULE_LIST, tasks=tasks)


@analysis_bp.route('/<pid>/bulk-analysis')
def bulk_analysis_list(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    return render_template('bulk_analysis.html', project=p, modules=BULK_MODULE_LIST, tasks=tasks)


@analysis_bp.route('/<pid>/analyze/<module_name>', methods=['GET', 'POST'])
def analyze(pid, module_name):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    mod_info = next((m for m in MODULE_LIST if m['name'] == module_name), None)
    if not mod_info:
        flash(f'模块 "{module_name}" 不存在', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    schema = PARAM_SCHEMAS.get(module_name, [])
    tasks = AnalysisTask.get_by_project(pid)
    completed_tasks = [t for t in tasks if t.status == 'completed' and t.output_adata_path]

    is_bulk = module_name in BULK_MODULE_NAMES

    if request.method == 'POST':
        params = {}
        for param in schema:
            val = request.form.get(param['key'])
            if param['type'] == 'number':
                params[param['key']] = float(val) if val else param['default']
            elif param['type'] == 'checkbox':
                params[param['key']] = request.form.get(param['key']) == 'on'
            else:
                params[param['key']] = val or param['default']
        # 处理自动检测的分组映射
        auto_mapping = request.form.get('_auto_group_mapping')
        if auto_mapping:
            try:
                params['_auto_group_mapping'] = json.loads(auto_mapping)
            except Exception:
                pass
        input_path = request.form.get('input_path', '')
        if not input_path:
            flash('请选择输入数据', 'danger')
            return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
        # 注入 _visualization 和 _filters 到 params
        viz_json = request.form.get('_visualization', '')
        filters_json = request.form.get('_filters', '')
        if viz_json:
            try:
                params['_visualization'] = json.loads(viz_json)
            except json.JSONDecodeError:
                pass
        if filters_json:
            try:
                params['_filters'] = json.loads(filters_json)
            except json.JSONDecodeError:
                pass
        task = AnalysisTask(project_id=pid, module_name=module_name,
                           params_json=json.dumps(params))
        task.save()
        p.status = 'processing'
        p.save()
        from worker import submit_task
        submit_task(task.id, pid, module_name, params,
                   os.path.join(Config.DATA_DIR, 'projects', pid), input_path)
        return redirect(url_for('results.task_detail', pid=pid, task_id=task.id))

    uploads_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'uploads')
    uploaded_files = []
    if os.path.isdir(uploads_dir):
        for f in os.listdir(uploads_dir):
            fpath = os.path.join(uploads_dir, f)
            if is_bulk:
                if f.endswith(('.h5ad', '.csv', '.txt', '.xlsx', '.xls')):
                    uploaded_files.append({'name': f, 'path': fpath})
            else:
                if f.endswith('.h5ad'):
                    uploaded_files.append({'name': f, 'path': fpath})

    # For enrichment module, auto-populate input_source with DEG results
    if module_name == 'bulk_enrichment':
        results_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'results')
        deg_csv = os.path.join(results_dir, 'bulk_deg_results.csv')
        if os.path.exists(deg_csv):
            for param in schema:
                if param['key'] == 'input_source':
                    param['default'] = deg_csv

    sidebar_modules = BULK_MODULE_LIST if is_bulk else SC_MODULE_LIST
    return render_template('analysis_select.html', project=p, module=mod_info,
                          schema=schema, completed_tasks=completed_tasks,
                          uploaded_h5ad=uploaded_files, all_modules=sidebar_modules,
                          module_display_map=MODULE_DISPLAY_MAP)
