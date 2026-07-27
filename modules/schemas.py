import os
import json

SC_MODULE_LIST = [
    {'name': 'qc', 'display': '质控', 'desc': 'MT/ribo/hb 过滤 + 双细胞检测 + 细胞周期评分 + 复杂度过滤'},
    {'name': 'normalize', 'display': '标准化', 'desc': '数据标准化（log1p / Pearson 残差）'},
    {'name': 'hvg', 'display': '高变异基因', 'desc': '选择高变异基因，支持批次感知和基因过滤'},
    {'name': 'dimred', 'display': '降维分析', 'desc': 'PCA, UMAP'},
    {'name': 'batch_correct', 'display': '批次校正', 'desc': 'Harmony, ComBat, SysVI'},
    {'name': 'clustering', 'display': '聚类分析', 'desc': 'Leiden 聚类'},
    {'name': 'subcluster', 'display': '子簇精细分析', 'desc': '选定一个簇进行重聚类、差异表达、热图和通路富集'},
    {'name': 'qc_reassess', 'display': 'QC 重新评估', 'desc': '聚类后检查 doublet 和 QC 指标，标记低质量簇'},
    {'name': 'annotation', 'display': '细胞注释', 'desc': '基于 Marker 的细胞类型注释'},
    {'name': 'sc_timecourse', 'display': '单细胞时序动态', 'desc': '按真实时间点分析样本级细胞组成与伪 bulk 基因动态'},
    {'name': 'deg', 'display': '差异表达', 'desc': '差异表达基因分析'},
    {'name': 'trajectory', 'display': '轨迹分析', 'desc': '拟时序分析'},
    {'name': 'proportion', 'display': '比例分析', 'desc': '细胞比例分析'},
    {'name': 'cell_communication', 'display': '细胞通讯', 'desc': '基于 LIANA 的细胞间通讯分析'},
]

BULK_MODULE_LIST = [
    {'name': 'bulk_qc', 'display': 'Bulk RNA-seq 质控', 'desc': '文库大小、基因检测、离群值过滤'},
    {'name': 'bulk_normalize', 'display': '数据标准化', 'desc': 'DESeq2 / CPM / 分位数标准化'},
    {'name': 'bulk_pca', 'display': 'PCA / UMAP', 'desc': 'PCA 和 UMAP 降维可视化，检查样本分组和批次效应'},
    {'name': 'bulk_deg', 'display': '差异表达分析', 'desc': '组间差异基因检测（火山图、MA图）'},
    {'name': 'bulk_heatmap', 'display': 'Bulk 热图可视化', 'desc': 'Top 差异基因热图、样本相关性热图'},
    {'name': 'bulk_enrichment', 'display': '通路富集', 'desc': 'GO/KEGG/WikiPathways 通路富集分析（ORA / GSEA）'},
    {'name': 'bulk_timecourse', 'display': '时序分析（可选）', 'desc': '多时间点差异基因检测 + 轨迹聚类（仅时序实验需要）'},
    {'name': 'bulk_deg_integration', 'display': '多组差异整合分析', 'desc': '多组比较结果整合：Upset 图、一致性评分、logFC 矩阵分析（仅多组比较需要）'},
]

MODULE_LIST = SC_MODULE_LIST + BULK_MODULE_LIST

MODULE_DISPLAY_MAP = {m['name']: m['display'] for m in MODULE_LIST}
MODULE_DISPLAY_MAP['convert_10x'] = '单细胞数据导入'

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
        {'key': 'hb_perc', 'label': '最大血红蛋白比例 %（0 = 不过滤）', 'type': 'number', 'default': 0, 'step': 1.0, 'help': '过滤血红蛋白基因比例高于此值的细胞。高比例通常表示红细胞污染。设为 0 表示不过滤。建议范围 5-10%。'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch', 'help': 'adata.obs 中标识实验批次的列名。用于分批次运行 Scrublet 双细胞检测。'},
        {'key': 'batch_adaptive_qc', 'label': '批次自适应 QC', 'type': 'checkbox', 'default': False, 'help': '按批次独立计算 MAD 阈值过滤，适用于批次间质量差异大的数据。'},
        {'key': 'mad_multiplier', 'label': 'MAD 倍数', 'type': 'number', 'default': 3.0, 'step': 0.5, 'help': '批次自适应 QC 的 MAD 倍数。越大越宽松。默认 3.0（约对应 3σ）。'},
        {'key': 'save_counts_layer', 'label': '保存原始 counts 层', 'type': 'checkbox', 'default': True, 'help': '在 QC 过滤前将原始表达矩阵保存到 adata.layers["counts"]，供下游标准化使用。'},
        {'key': 'show_qc_filter_summary', 'label': '生成过滤前后 QC 对比图', 'type': 'checkbox', 'default': True, 'help': '输出 QC 过滤前后细胞数、基因数和核心 QC 指标的对比图，用于检查过滤强度是否合理。'},
        {'key': 'show_doublet_histogram', 'label': '生成 Doublet score 直方图', 'type': 'checkbox', 'default': True, 'help': '输出 Scrublet doublet score 分布图，用于判断双细胞阈值和残留双细胞风险。'},
    ],
    'normalize': [
        {'key': 'method', 'label': '标准化方法', 'type': 'select', 'options': ['log1p', 'pearson_residuals'], 'default': 'log1p', 'help': 'log1p：标准 log1p CPM（shiftlog），适合大多数分析。pearson_residuals：Pearson 残差标准化，对技术噪声更鲁棒。'},
        {'key': 'target_sum', 'label': '标准化目标总数', 'type': 'number', 'default': 10000, 'help': '每个细胞标准化后的总计数目标。10000 为 scanpy 默认值。'},
        {'key': 'clip_values', 'label': '裁剪 Pearson 残差', 'type': 'checkbox', 'default': True, 'help': '仅 pearson_residuals 模式生效。裁剪残差到 ±√n 范围，减少极端值影响。'},
        {'key': 'show_expression_distribution', 'label': '生成表达值分布图', 'type': 'checkbox', 'default': True, 'help': '展示标准化后表达值分布，并与 log1p 原始 counts 对比，用于发现标准化异常或极端值。'},
    ],
    'hvg': [
        {'key': 'n_top_genes', 'label': '高变异基因数量', 'type': 'number', 'default': 2000, 'help': '选择的高变异基因数量。2000 为标准值，适合大多数分析。基因数过少会丢失生物学信号，过多会引入噪声。'},
        {'key': 'batch_key', 'label': '批次列名（可选）', 'type': 'text', 'default': '', 'help': '用于批次感知 HVG 选择。留空则不进行批次校正。'},
        {'key': 'hvg_flavor', 'label': 'HVG 选择方法', 'type': 'select', 'options': ['seurat_v3', 'cell_ranger', 'seurat'], 'default': 'seurat_v3', 'help': 'seurat_v3：基于方差稳定的 HVG（推荐）。cell_ranger：Cell Ranger 方法。seurat：Seurat v1 方法（dispersion-based）。'},
        {'key': 'batch_hvg_strategy', 'label': '批次 HVG 合并策略', 'type': 'select', 'options': ['intersection', 'union'], 'default': 'intersection', 'help': '仅批次列名非空时生效。intersection：取各批次 HVG 的交集（更严格）。union：取并集（更宽松）。'},
        {'key': 'exclude_mt_genes', 'label': '排除线粒体基因', 'type': 'checkbox', 'default': False, 'help': '从 HVG 中排除线粒体基因（MT- 前缀）。'},
        {'key': 'exclude_cc_genes', 'label': '排除细胞周期基因', 'type': 'checkbox', 'default': False, 'help': '从 HVG 中排除 S 期和 G2M 期细胞周期基因。'},
        {'key': 'force_include_genes', 'label': '强制包含基因（可选）', 'type': 'textarea', 'default': '', 'help': '强制包含在 HVG 中的基因名，逗号或换行分隔。无论是否被选为 HVG 都会保留。'},
        {'key': 'cc_scoring', 'label': '细胞周期评分', 'type': 'checkbox', 'default': False, 'help': '计算 S 期和 G2M 期评分，存入 obs。'},
        {'key': 'regress_cc', 'label': '回归去除细胞周期', 'type': 'checkbox', 'default': False, 'help': '回归去除细胞周期效应（需先开启细胞周期评分）。用于消除细胞周期对下游分析的干扰。'},
        {'key': 'show_hvg_rank_plot', 'label': '生成 HVG Rank 图', 'type': 'checkbox', 'default': True, 'help': '按变异度排名展示基因，并标注 Top HVG，用于检查 HVG 选择是否合理。'},
    ],
    'dimred': [
        {'key': 'n_comps', 'label': 'PCA 主成分数量', 'type': 'number', 'default': 50, 'help': 'PCA 主成分数量。通常 30-50 即可捕获大部分方差。'},
        {'key': 'auto_n_comps', 'label': '自动选 PC', 'type': 'select', 'options': ['none', 'elbow', 'kneedle'], 'default': 'none', 'help': '自动选择主成分数。none：使用固定值。elbow：肘部法则。kneedle：Kneedle 算法（需安装 kneed 包）。'},
        {'key': 'umap_n_neighbors', 'label': 'UMAP 邻居数', 'type': 'number', 'default': 15, 'help': '构建 UMAP 邻居图的邻居数。值越大全局结构越清晰，局部结构越模糊。'},
        {'key': 'umap_min_dist', 'label': 'UMAP 最小距离', 'type': 'number', 'default': 0.5, 'step': 0.05, 'help': 'UMMAP 嵌入中点之间的最小距离。越小越紧密，越大越松散。范围 0.0-1.0。'},
        {'key': 'umap_metric', 'label': '距离度量', 'type': 'select', 'options': ['euclidean', 'cosine', 'correlation', 'manhattan'], 'default': 'euclidean', 'help': '计算邻居的距离度量。euclidean：欧氏距离（常用）。cosine：余弦距离（适合稀疏数据）。correlation：相关距离。'},
        {'key': 'umap_spread', 'label': 'UMAP spread', 'type': 'number', 'default': 1.0, 'step': 0.1, 'help': 'UMAP 嵌入的有效尺度。通常与 min_dist 配合调整。'},
        {'key': 'enable_tsne', 'label': '额外生成 t-SNE', 'type': 'checkbox', 'default': False, 'help': '除 UMAP 外额外计算 t-SNE 嵌入。t-SNE 更注重局部结构，但计算更慢。'},
        {'key': 'tsne_perplexity', 'label': 't-SNE 困惑度', 'type': 'number', 'default': 30, 'help': 't-SNE 困惑度，大致表示有效邻居数。通常 5-50，大数据集可增大。'},
        {'key': 'tsne_learning_rate', 'label': 't-SNE 学习率', 'type': 'number', 'default': 1000, 'help': 't-SNE 学习率。通常 100-1000，默认 1000。'},
        {'key': 'use_mde', 'label': '使用 MDE（加速 UMAP）', 'type': 'checkbox', 'default': False, 'help': '使用 Minimum Dystortion Embedding 替代标准 UMAP，速度更快但结果略有差异。'},
        {'key': 'show_pca_scatter', 'label': '生成 PCA Scatter 图', 'type': 'checkbox', 'default': True, 'help': '展示 PC1/PC2 散点图，用于在 UMAP 前检查离群细胞、批次结构或样本结构。'},
    ],
    'batch_correct': [
        {'key': 'method', 'label': '校正方法', 'type': 'select', 'options': ['harmony', 'combat', 'bbknn', 'scanorama', 'sysvi', 'scvi'], 'default': 'harmony', 'help': '批次整合方法。Harmony：PCA 空间校正；ComBat：Scanpy ComBat PCA；BBKNN：邻居图整合；Scanorama：MNN/全局整合；SysVI/scVI：深度生成模型，可 CPU/GPU 运行。'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch', 'help': 'adata.obs 中标识批次的列名。'},
        {'key': 'n_pcs', 'label': '主成分数量', 'type': 'number', 'default': 50, 'show_if': {'method': ['harmony', 'combat', 'scanorama']}, 'help': '用于 Harmony、ComBat 和 Scanorama 的主成分数量。BBKNN 使用已有 PCA 图；SysVI/scVI 使用潜在空间维度。'},
        {'key': 'max_epochs', 'label': '最大迭代轮数（仅 SysVI/scVI）', 'type': 'number', 'default': 60, 'show_if': {'method': ['sysvi', 'scvi']}, 'help': '深度学习方法的最大训练轮数。CPU 环境建议先用 40-80。'},
        {'key': 'harmony_theta', 'label': 'Harmony theta', 'type': 'number', 'default': 2.0, 'step': 0.5, 'show_if': {'method': 'harmony'}, 'help': 'Harmony 多样性惩罚。值越大强制批次混合越强。默认 2.0。'},
        {'key': 'harmony_lambda', 'label': 'Harmony lambda', 'type': 'number', 'default': 1.0, 'step': 0.1, 'show_if': {'method': 'harmony'}, 'help': 'Harmony 正则化强度。默认 1.0。'},
        {'key': 'harmony_max_iter', 'label': 'Harmony 最大迭代', 'type': 'number', 'default': 20, 'show_if': {'method': 'harmony'}, 'help': 'Harmony 最大校正迭代数。'},
        {'key': 'bbknn_neighbors_within_batch', 'label': 'BBKNN 每批次邻居数', 'type': 'number', 'default': 3, 'show_if': {'method': 'bbknn'}, 'help': 'BBKNN 在每个批次内搜索的邻居数。默认 3；如果某个批次细胞数更少，平台会自动下调到最小批次细胞数并在 summary 中记录。'},
        {'key': 'scvi_n_latent', 'label': 'scVI 潜在维度', 'type': 'number', 'default': 30, 'show_if': {'method': ['sysvi', 'scvi']}, 'help': 'scVI/SysVI 潜在空间维度。默认 30。'},
        {'key': 'scvi_n_hidden', 'label': 'scVI 隐藏层大小', 'type': 'number', 'default': 128, 'show_if': {'method': ['sysvi', 'scvi']}, 'help': 'scVI/SysVI 隐藏层神经元数。默认 128。'},
        {'key': 'scvi_n_layers', 'label': 'scVI 层数', 'type': 'number', 'default': 1, 'show_if': {'method': ['sysvi', 'scvi']}, 'help': 'scVI/SysVI 编码器/解码器层数。默认 1。'},
        {'key': 'scvi_dropout_rate', 'label': 'scVI dropout', 'type': 'number', 'default': 0.1, 'step': 0.05, 'show_if': {'method': ['sysvi', 'scvi']}, 'help': 'scVI/SysVI dropout 率。默认 0.1。'},
        {'key': 'scvi_learning_rate', 'label': 'scVI 学习率', 'type': 'number', 'default': 0.001, 'show_if': {'method': 'scvi'}, 'help': 'scVI 学习率。默认 0.001。'},
        {'key': 'sysvi_cycle_weight', 'label': 'SysVI cycle 权重', 'type': 'number', 'default': 5.0, 'step': 0.5, 'show_if': {'method': 'sysvi'}, 'help': 'SysVI latent cycle-consistency 权重。越高通常整合越强。'},
        {'key': 'sysvi_kl_weight', 'label': 'SysVI KL 权重', 'type': 'number', 'default': 1.0, 'step': 0.1, 'show_if': {'method': 'sysvi'}, 'help': 'SysVI KL loss 权重。降低可能保留更多生物差异。'},
        {'key': 'sysvi_prior', 'label': 'SysVI prior', 'type': 'select', 'options': ['vamp', 'standard_normal'], 'default': 'vamp', 'show_if': {'method': 'sysvi'}, 'help': 'SysVI prior 类型。默认 vamp。'},
        {'key': 'sysvi_n_prior_components', 'label': 'SysVI prior 组件数', 'type': 'number', 'default': 5, 'show_if': {'method': 'sysvi', 'sysvi_prior': 'vamp'}, 'help': 'VampPrior 组件数。'},
        {'key': 'evaluate_correction', 'label': '评估整合效果', 'type': 'checkbox', 'default': True, 'help': '计算批次整合关键指标并生成指标图表：batch ASW、cluster batch entropy、最大批次占比、邻居混合、图连通性等。'},
        {'key': 'evaluation_cluster_key', 'label': '评价分群列（可选）', 'type': 'text', 'default': '', 'show_if': {'evaluate_correction': True}, 'help': '用于计算 cluster 层面批次混合的 obs 列。留空时优先用 leiden；没有则自动生成 batch_eval_leiden。'},
        {'key': 'evaluation_resolution', 'label': '评价分群 resolution', 'type': 'number', 'default': 0.8, 'step': 0.1, 'show_if': {'evaluate_correction': True}, 'help': '自动生成 batch_eval_leiden 时使用的 resolution。'},
        {'key': 'evaluation_sample_size', 'label': 'ASW 抽样细胞数', 'type': 'number', 'default': 10000, 'show_if': {'evaluate_correction': True}, 'help': '计算 silhouette/ASW 时的抽样细胞数，避免大数据过慢。'},
        {'key': 'bio_label_key', 'label': '生物标签列（可选）', 'type': 'text', 'default': '', 'show_if': {'evaluate_correction': True}, 'help': '用于估计生物结构保留的标签列，如 celltype、annotation 或参考标签。留空时自动尝试 celltype/reference_celltype/leiden。'},
    ],
    'clustering': [
        {'key': 'resolutions', 'label': '聚类分辨率（逗号分隔）', 'type': 'text', 'default': '0.6,0.8,1.0', 'help': '聚类分辨率，多个值用逗号分隔。值越大聚类越细。'},
        {'key': 'n_neighbors', 'label': '邻居数量', 'type': 'number', 'default': 15, 'help': '构建 KNN 图时的邻居数量。'},
        {'key': 'clustering_method', 'label': '聚类算法', 'type': 'select', 'options': ['leiden', 'louvain'], 'default': 'leiden', 'help': '聚类算法。Leiden：社区检测算法，推荐。Louvain：经典社区检测。'},
        {'key': 'n_iterations', 'label': 'Leiden 迭代次数', 'type': 'number', 'default': 2, 'help': 'Leiden 算法迭代次数。-1 为运行至收敛。'},
        {'key': 'distance_metric', 'label': '距离度量', 'type': 'select', 'options': ['euclidean', 'cosine', 'correlation', 'manhattan'], 'default': 'euclidean', 'help': '邻居图的距离度量。'},
        {'key': 'use_corrected', 'label': '使用校正后表示', 'type': 'checkbox', 'default': True, 'help': '优先使用批次校正后的嵌入（如有）。'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch', 'help': '用于生成 cluster 批次组成图的 obs 列名，建议选择 batch/sample 等分类列；连续高基数列会自动跳过。数据中不存在该列时也会自动跳过。'},
        {'key': 'auto_select_resolution', 'label': '自动选择最优分辨率', 'type': 'checkbox', 'default': False, 'help': '使用聚类质量指标自动选择最优分辨率。'},
        {'key': 'resolution_metric', 'label': '评估指标', 'type': 'select', 'options': ['silhouette', 'calinski', 'davies_bouldin'], 'default': 'silhouette', 'help': '自动选择分辨率时的质量评估指标。'},
        {'key': 'primary_resolution', 'label': '主分辨率（可选）', 'type': 'text', 'default': '', 'help': '指定最终写入 leiden 的主分辨率，例如 0.8。留空则使用首个分辨率或自动选择结果。'},
        {'key': 'marker_selection_method', 'label': 'Marker 排名方法', 'type': 'select', 'options': ['wilcoxon', 'logreg'], 'default': 'wilcoxon', 'help': '按 cluster 独立排名 marker。Wilcoxon 默认提供校正 P 值；Logistic regression 主要按效应排名，小数据集无显著 P 值时会明确标记为探索性候选。'},
        {'key': 'marker_rank_genes', 'label': '每簇候选排名数', 'type': 'number', 'default': 200, 'step': 10, 'help': '每个 cluster 先保留多少个候选基因再按表达比例、logFC 和特异性筛选。'},
        {'key': 'marker_padj_cutoff', 'label': 'Marker 校正 P 值阈值', 'type': 'number', 'default': 0.05, 'step': 0.01, 'help': '数据驱动 marker 的 adjusted P-value 上限。'},
        {'key': 'marker_min_pct', 'label': 'Marker 最小表达比例', 'type': 'number', 'default': 0.1, 'step': 0.05, 'help': '基因在目标 cluster 中的最小检测比例，避免由极少数细胞驱动。'},
        {'key': 'marker_min_delta_pct', 'label': 'Marker 最小特异性差值', 'type': 'number', 'default': 0.05, 'step': 0.05, 'help': '目标 cluster 与其他 cluster 的检测比例差值下限。'},
        {'key': 'marker_min_per_cluster', 'label': '每簇最少 Marker 数', 'type': 'number', 'default': 2, 'step': 1, 'help': '每个 cluster 期望展示的最少 marker；小数据集不足时会保留实际可用数量。'},
        {'key': 'marker_max_per_cluster', 'label': '每簇最多 Marker 数', 'type': 'number', 'default': 5, 'step': 1, 'help': '每个 cluster 最多展示的 marker 数，避免 dotplot 过密；不同簇之间默认全局去重。'},
        {'key': 'show_labeled_umap', 'label': '生成带标签 Cluster UMAP', 'type': 'checkbox', 'default': True, 'help': '在主分辨率 UMAP 上显示 cluster 编号，便于人工复核分群是否符合预期。'},
        {'key': 'show_resolution_sankey', 'label': '生成分辨率流向图', 'type': 'checkbox', 'default': True, 'help': '用 Sankey 图展示不同 Leiden 分辨率之间的簇拆分关系，辅助选择合适分辨率。'},
        {'key': 'show_cluster_size_bar', 'label': '生成 Cluster 细胞数图', 'type': 'checkbox', 'default': True, 'help': '展示主分辨率每个 cluster 的细胞数量，用于识别过小簇、过度分裂或不均衡分群。'},
        {'key': 'show_cluster_batch_composition', 'label': '生成 Cluster 批次组成图', 'type': 'checkbox', 'default': True, 'help': '展示每个 cluster 的 batch/sample 组成比例，用于发现单一批次支配的分群。'},
    ],
    'subcluster': [
        {'key': 'source_cluster_key', 'label': '来源聚类列', 'type': 'text', 'default': 'leiden', 'help': '原始 AnnData.obs 中要精细拆分的聚类列，通常为 leiden。'},
        {'key': 'target_cluster', 'label': '目标簇编号', 'type': 'text', 'default': '', 'help': '要单独重聚类的簇编号，例如 3。必须与来源聚类列中的值完全一致。'},
        {'key': 'min_cells', 'label': '最小细胞数', 'type': 'number', 'default': 30, 'step': 1, 'help': '目标簇少于此数量时停止，避免对过少细胞产生不稳定子簇。'},
        {'key': 'n_neighbors', 'label': '子簇邻居数', 'type': 'number', 'default': 15, 'step': 1, 'help': '仅在目标簇细胞中重新构建 KNN 图；平台会自动限制为小于细胞数。'},
        {'key': 'resolution', 'label': '子簇聚类分辨率', 'type': 'number', 'default': 0.8, 'step': 0.1, 'help': '值越大拆分越细。建议从 0.4-1.2 多次比较后确定。'},
        {'key': 'clustering_method', 'label': '子簇聚类算法', 'type': 'select', 'options': ['leiden', 'louvain'], 'default': 'leiden', 'help': 'Leiden 为推荐的子簇社区发现方法。'},
        {'key': 'n_iterations', 'label': 'Leiden 迭代次数', 'type': 'number', 'default': 2, 'step': 1, 'show_if': {'clustering_method': 'leiden'}, 'help': 'Leiden 迭代次数；-1 表示收敛。'},
        {'key': 'distance_metric', 'label': '距离度量', 'type': 'select', 'options': ['euclidean', 'cosine', 'correlation'], 'default': 'euclidean', 'help': '子簇 KNN 图的距离度量。'},
        {'key': 'umap_min_dist', 'label': '子簇 UMAP 最小距离', 'type': 'number', 'default': 0.4, 'step': 0.05, 'help': '仅影响子簇 UMAP 展示，不改变子簇身份。'},
        {'key': 'deg_method', 'label': '子簇差异方法', 'type': 'select', 'options': ['wilcoxon', 't-test', 'logreg'], 'default': 'wilcoxon', 'help': '子簇 marker 差异表达方法。'},
        {'key': 'n_genes', 'label': '每子簇 DEG 数', 'type': 'number', 'default': 50, 'step': 10, 'help': '导出的每个子簇 marker/DEG 数量。提高此值可保留更多候选基因。'},
        {'key': 'marker_heatmap_top_n', 'label': '热图每子簇 Marker 数', 'type': 'number', 'default': 5, 'step': 1, 'help': '每个子簇纳入热图的 Top marker 数。'},
        {'key': 'run_enrichment', 'label': '运行子簇通路富集', 'type': 'checkbox', 'default': True, 'help': '对每个子簇的显著上调 marker 独立运行 Enrichr；网络不可用时不会影响其他结果。'},
        {'key': 'enrichment_database', 'label': '富集数据库', 'type': 'select', 'options': ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022'], 'default': 'GO_Biological_Process_2023', 'show_if': {'run_enrichment': True}, 'help': '每个子簇使用的 Enrichr 基因集库。'},
        {'key': 'organism', 'label': '物种', 'type': 'select', 'options': ['Human', 'Mouse'], 'default': 'Human', 'show_if': {'run_enrichment': True}, 'help': '富集服务所用物种。'},
        {'key': 'enrichment_pval_cutoff', 'label': '富集校正 P 值阈值', 'type': 'number', 'default': 0.05, 'step': 0.01, 'show_if': {'run_enrichment': True}, 'help': '仅使用达到此 Adjusted P-value 阈值的 marker 进行富集。'},
        {'key': 'enrichment_top_n', 'label': '每子簇展示 Top 通路数', 'type': 'number', 'default': 10, 'step': 1, 'show_if': {'run_enrichment': True}, 'help': '通路气泡图中每个子簇展示的最多通路数量。'},
    ],
    'qc_reassess': [
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden', 'help': '用于评估的分类聚类列名（通常为 leiden）；连续 QC 指标或高基数列会自动回退。'},
        {'key': 'doublet_threshold', 'label': 'Doublet 比例阈值', 'type': 'number', 'default': 0.3, 'step': 0.05, 'help': 'Doublet 比例高于此值的簇标记为低质量。'},
        {'key': 'mt_threshold', 'label': 'MT 比例阈值', 'type': 'number', 'default': 15.0, 'step': 1.0, 'help': '平均 MT 比例高于此值的簇标记为低质量。'},
        {'key': 'ribosomal_threshold', 'label': '核糖体比例阈值（0 = 不检查）', 'type': 'number', 'default': 0, 'step': 1.0, 'help': '平均核糖体比例高于此值的簇标记为低质量。0 表示不检查。'},
        {'key': 'min_cells_per_cluster', 'label': '最小细胞数', 'type': 'number', 'default': 10, 'help': '细胞数低于此值的簇标记为低质量。'},
        {'key': 'auto_remove', 'label': '自动移除低质量簇', 'type': 'checkbox', 'default': False, 'help': '自动从数据中移除标记为低质量的簇。'},
        {'key': 'show_qc_umap_panel', 'label': '生成 QC 指标 UMAP 面板', 'type': 'checkbox', 'default': True, 'help': '把 MT%、检测基因数、总 counts 和 doublet score 映射到 UMAP，用于定位低质量区域或疑似污染簇。'},
        {'key': 'show_cluster_qc_bar', 'label': '生成按簇 QC 汇总图', 'type': 'checkbox', 'default': True, 'help': '按 cluster 展示细胞数、检测基因数、MT% 和 doublet fraction，低质量簇用颜色标记。'},
    ],
    'annotation': [
        {'key': 'method', 'label': '注释方法', 'type': 'select', 'options': ['multi_evidence', 'auto_marker', 'manual'], 'default': 'multi_evidence', 'help': 'multi_evidence：Marker 规则、逐细胞 Top1/Top2 和 Cluster 一致性复核（推荐）；auto_marker：仅规则打分；manual：手工映射。CellTypist 作为可选参考证据，不直接替换 Marker 标签。'},
        {'key': 'annotation_version', 'label': '注释版本', 'type': 'text', 'default': 'v1', 'help': '用于区分自动注释、人工修订和不同 Marker 方案；不会覆盖输出中的旧 h5ad。'},
        {'key': 'annotation_comment', 'label': '注释备注（可选）', 'type': 'textarea', 'default': '', 'help': '记录本次注释的实验背景、人工判断或需要后续复核的事项。'},
        {'key': 'state_score_threshold', 'label': '细胞状态评分阈值', 'type': 'number', 'default': 0.35, 'step': 0.05, 'min': 0.0, 'max': 1.0, 'help': '仅用于标记 Cycling、IFN、Stress 等状态；不改变 cell_type。多个状态可同时进入 cell_state_flags。'},
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden', 'help': '用于分组的分类聚类列名（通常为 leiden）；连续 QC 指标不能作为 cluster。'},
        {'key': 'resolution', 'label': 'Leiden 分辨率', 'type': 'text', 'default': '0.8', 'help': '对应的 Leiden 分辨率，用于定位正确的聚类列。'},
        {'key': 'marker_set', 'label': 'Marker 基因集', 'type': 'select', 'options': ['Universal', 'Organoid', 'TME', 'Immune', 'Blood', 'PBMC'], 'default': 'Universal', 'help': 'Universal：未知组织的通用大谱系初注释（推荐起点）；Organoid：按类器官类型选择发育/组织谱系 marker；TME/Immune/Blood/PBMC：已知场景的细分 marker 集。初注释后可对子簇使用场景集或自定义 marker 精修。'},
        {'key': 'use_celltypist_reference', 'label': '启用 CellTypist 参考交叉验证', 'type': 'checkbox', 'default': False, 'show_if': {'method': ['multi_evidence', 'auto_marker']}, 'help': '使用本地已下载的人类 CellTypist 模型生成独立参考标签、置信度和冲突状态；只写入参考证据，不覆盖最终 Marker 注释。缺少依赖或模型时自动回退并记录警告。'},
        {'key': 'celltypist_model', 'label': 'CellTypist 人类模型', 'type': 'select', 'options': ['Immune_All_Low.pkl', 'Immune_All_High.pkl', 'Cells_Intestinal_Tract.pkl', 'Developing_Human_Organs.pkl', 'Developing_Human_Brain.pkl', 'Cells_Fetal_Lung.pkl', 'Cells_Lung_Airway.pkl', 'Human_Lung_Atlas.pkl', 'Nuclei_Lung_Airway.pkl', 'Healthy_Human_Liver.pkl', 'Adult_Human_PancreaticIslet.pkl', 'Fetal_Human_Pancreas.pkl', 'Healthy_Adult_Heart.pkl', 'Adult_Human_Skin.pkl', 'Fetal_Human_Skin.pkl', 'Adult_Human_Vascular.pkl', 'Cells_Adult_Breast.pkl', 'Pan_Fetal_Human.pkl'], 'default': 'Immune_All_Low.pkl', 'show_if': {'use_celltypist_reference': True}, 'help': '只从项目 data/references/celltypist 读取，不在分析期间自动下载；应选择与组织和物种匹配的模型。'},
        {'key': 'celltypist_mode', 'label': 'CellTypist 匹配模式', 'type': 'select', 'options': ['prob match', 'best match'], 'default': 'prob match', 'show_if': {'use_celltypist_reference': True}, 'help': 'prob match 可输出 Unknown/多标签，适合类器官参考复核；best match 总会选择一个最相近标签，需谨慎解释。'},
        {'key': 'celltypist_p_threshold', 'label': 'CellTypist 概率阈值', 'type': 'number', 'default': 0.5, 'step': 0.05, 'min': 0.0, 'max': 1.0, 'show_if': {'use_celltypist_reference': True}, 'help': 'prob match 的最低概率阈值，同时用于标记低置信度参考结果；低于阈值不参与一致性判断。'},
        {'key': 'celltypist_majority_voting', 'label': '启用 CellTypist 簇多数投票', 'type': 'checkbox', 'default': False, 'show_if': {'use_celltypist_reference': True}, 'help': '需要额外 over-clustering；默认关闭，避免在类器官过渡态中把混合簇强行平滑。'},
        {'key': 'organoid_type', 'label': '类器官类型', 'type': 'select', 'options': ['intestinal', 'cerebral', 'kidney', 'liver', 'lung', 'pancreatic', 'cardiac'], 'default': 'intestinal', 'show_if': {'marker_set': 'Organoid'}, 'help': '仅在 Marker 基因集选择 Organoid 时生效。面板用于第一轮候选注释，需结合实验阶段、物种和 marker 覆盖度人工复核；不要把类器官面板当作跨协议通用真值。'},
        {'key': 'maturity_time_key', 'label': '成熟度时间列（可选）', 'type': 'text', 'default': '', 'show_if': {'marker_set': 'Organoid'}, 'help': '填写 adata.obs 中的培养天数/时间点列名，如 culture_day、day、timepoint；留空时自动寻找这些列。不会把时间点直接当作成熟标签，而是与表达成熟度并列展示。'},
        {'key': 'custom_markers', 'label': '自定义 Marker（可选）', 'type': 'textarea', 'default': '', 'help': 'auto_marker 模式：CellType:GENE1,GENE2 格式。manual 模式：ClusterID:CellType 格式。'},
        {'key': 'negative_markers', 'label': '自定义负向 Marker（可选）', 'type': 'textarea', 'default': '', 'show_if': {'method': ['auto_marker', 'multi_evidence']}, 'help': '格式：CellType:GENE1,GENE2；这些基因表达会降低候选类型分数，但不会单独强制改成 Unknown。留空使用内置互斥 panel。'},
        {'key': 'negative_marker_weight', 'label': '负向 Marker 惩罚权重', 'type': 'number', 'default': 0.5, 'step': 0.1, 'min': 0.0, 'max': 2.0, 'show_if': {'method': ['auto_marker', 'multi_evidence']}, 'help': '建议 0.3-0.8。权重过高可能误伤真实过渡态，先查看负向证据和复核表再调整。'},
        {'key': 'doublet_score_threshold', 'label': '疑似 Doublet 阈值', 'type': 'number', 'default': 0.3, 'step': 0.05, 'min': 0.0, 'max': 1.0, 'show_if': {'method': ['auto_marker', 'multi_evidence']}, 'help': '两个独立 marker 模块分数接近时提高；仅标记 suspect_doublet/review，不自动删除细胞。'},
        {'key': 'ambient_score_threshold', 'label': '环境 RNA 信号阈值', 'type': 'number', 'default': 0.35, 'step': 0.05, 'min': 0.0, 'max': 1.0, 'show_if': {'method': ['auto_marker', 'multi_evidence']}, 'help': '基于高普遍性 marker 的启发式信号；没有 empty droplets 时只能用于提示复核，不能当作定量去污染结果。'},
        {'key': 'ambient_prevalence', 'label': '环境候选基因普遍性阈值', 'type': 'number', 'default': 0.5, 'step': 0.05, 'min': 0.0, 'max': 1.0, 'show_if': {'method': ['auto_marker', 'multi_evidence']}, 'help': '某 marker 在超过该比例细胞中被检测到时，纳入环境 RNA 启发式检查。'},
        {'key': 'marker_selection_method', 'label': 'Marker 排名方法', 'type': 'select', 'options': ['wilcoxon', 'logreg'], 'default': 'wilcoxon', 'help': '每个 cluster 独立排名 dotplot 基因。Wilcoxon 默认提供校正 P 值；Logistic regression 主要按效应排名，小数据集无显著 P 值时会标记为探索性候选。'},
        {'key': 'marker_rank_genes', 'label': '每簇候选排名数', 'type': 'number', 'default': 200, 'step': 10, 'help': '每个 cluster 先保留多少个候选基因再按表达比例、logFC 和特异性筛选。'},
        {'key': 'marker_padj_cutoff', 'label': 'Marker 校正 P 值阈值', 'type': 'number', 'default': 0.05, 'step': 0.01, 'help': '数据驱动 marker 的 adjusted P-value 上限。'},
        {'key': 'marker_min_pct', 'label': 'Marker 最小表达比例', 'type': 'number', 'default': 0.1, 'step': 0.05, 'help': '基因在目标 cluster 中的最小检测比例，避免由极少数细胞驱动。'},
        {'key': 'marker_min_delta_pct', 'label': 'Marker 最小特异性差值', 'type': 'number', 'default': 0.05, 'step': 0.05, 'help': '目标 cluster 与其他 cluster 的检测比例差值下限。'},
        {'key': 'marker_min_per_cluster', 'label': '每簇最少 Marker 数', 'type': 'number', 'default': 2, 'step': 1, 'help': '每个 cluster 期望展示的最少 marker；小数据集不足时会保留实际可用数量。'},
        {'key': 'marker_max_per_cluster', 'label': '每簇最多 Marker 数', 'type': 'number', 'default': 5, 'step': 1, 'help': '每个 cluster 最多展示的 marker 数，避免 dotplot 过密；不同簇之间默认全局去重。'},
        {'key': 'min_markers_per_type', 'label': '每类型最少可用 Marker 数', 'type': 'number', 'default': 2, 'step': 1, 'show_if': {'method': 'auto_marker'}, 'help': '当前数据中命中少于此数量的类型不参与打分，避免因基因面板缺失而误注释。'},
        {'key': 'min_annotation_score', 'label': '自动注释最低 Marker 得分', 'type': 'number', 'default': 0.0, 'step': 0.05, 'show_if': {'method': 'auto_marker'}, 'help': '大于 0 时，最高 marker score 低于阈值的细胞标为 Unknown。建议先查看 Marker 覆盖度和得分热图后调整。'},
        {'key': 'cluster_agreement_threshold', 'label': 'Cluster 最低标签一致率', 'type': 'number', 'default': 0.6, 'step': 0.05, 'show_if': {'method': 'multi_evidence'}, 'help': '同一 cluster 内 Marker 标签多数比例低于此值时标记为 Unknown，避免将混杂 cluster 强行命名。'},
        {'key': 'confidence_method', 'label': '置信度方法', 'type': 'select', 'options': ['none', 'entropy', 'score_margin'], 'default': 'none', 'help': '注释置信度评估方法。entropy：基于评分熵。score_margin：基于最高分与次高分差距。'},
        {'key': 'mark_unknown', 'label': '低置信度自动标 Unknown', 'type': 'checkbox', 'default': False, 'help': '默认关闭：低置信度只标记为 review。开启后，只有低置信度且 cluster 内标签不一致时才自动标为 Unknown，避免广泛 marker panel 的熵值偏低造成大面积误删。'},
        {'key': 'confidence_cutoff', 'label': '低置信度阈值', 'type': 'number', 'default': 0.2, 'min': 0.0, 'max': 1.0, 'step': 0.05, 'show_if': {'confidence_method': ['entropy', 'score_margin']}, 'help': '仅在开启自动标 Unknown 时生效；建议结合逐簇复核表调整，不建议单独用它否定 cluster 内一致的标签。'},
        {'key': 'merge_similar_threshold', 'label': '相似簇合并阈值（0 = 不合并）', 'type': 'number', 'default': 0, 'step': 0.05, 'help': '相似度高于此值的相邻簇合并为同一细胞类型。0 表示不合并。'},
        {'key': 'show_celltype_composition', 'label': '生成细胞类型组成图', 'type': 'checkbox', 'default': True, 'help': '展示每种注释细胞类型的数量和比例，用于检查注释组成和样本结构。'},
        {'key': 'show_marker_score_heatmap', 'label': '生成 Marker score 热图', 'type': 'checkbox', 'default': True, 'help': '按 cluster 展示各细胞类型 marker score 的相对强弱，用于解释自动注释依据。'},
        {'key': 'show_marker_expression_violin', 'label': '生成 Marker 表达验证图', 'type': 'checkbox', 'default': True, 'help': '按注释细胞类型展示核心 marker 表达分布，用于人工确认注释是否符合生物学预期。'},
        {'key': 'show_annotation_score_umap', 'label': '生成注释置信度 UMAP', 'type': 'checkbox', 'default': True, 'help': '把 annotation confidence 或 score margin 映射到 UMAP，用于定位低置信度区域和可能需要重分群的细胞。'},
    ],
    'sc_timecourse': [
        {'key': 'timepoint_key', 'label': '真实时间点列名', 'type': 'text', 'default': 'timepoint', 'help': '必填。adata.obs 中的真实采样时间列，如 timepoint、day、hour。不能用 DPT pseudotime 代替。至少 3 个时间点。'},
        {'key': 'time_order', 'label': '时间点顺序（可选）', 'type': 'text', 'default': '', 'help': '逗号分隔，如 D0,D3,D7。D0/D3/D7、0h/12h 等会自动排序；任意标签请显式指定，以免 D10 排在 D2 前。'},
        {'key': 'confirm_batch_is_biological_timepoint', 'label': '确认 batch 列确为真实采样时间', 'type': 'checkbox', 'default': False, 'help': '仅当 timepoint_key 填 batch 且该列实际记录采样时间时勾选。技术建库/测序 batch 不能产生时序 p 值/FDR。'},
        {'key': 'sample_key', 'label': '生物学重复 / 样本列名', 'type': 'text', 'default': 'sample_id', 'help': '强烈建议提供，如 sample_id 或 library_id。每个样本必须只属于一个时间点和条件；缺失、近似每细胞唯一或重复不足时仅输出描述性图，不报告 p 值/FDR。'},
        {'key': 'confirm_batch_is_biological_sample', 'label': '确认 batch 实为独立生物学样本', 'type': 'checkbox', 'default': False, 'help': '仅当 sample_key 填 batch 且该列确实是独立生物学样本时勾选。技术建库/测序 batch 不能作为统计重复。'},
        {'key': 'condition_key', 'label': '条件列名（可选）', 'type': 'text', 'default': '', 'help': '如 treatment、condition。填写后分别展示每个条件的时间曲线；当前版本不进行 time × condition 交互检验，不能据此声称处理改变时间轨迹。'},
        {'key': 'celltype_key', 'label': '细胞类型列名', 'type': 'text', 'default': 'celltype', 'help': '优先使用已注释 celltype；不存在时才回退 annotation 或 leiden。'},
        {'key': 'min_replicates_per_timepoint', 'label': '每时间点最少生物学重复', 'type': 'number', 'default': 2, 'min': 2, 'step': 1, 'help': '每个条件 × 时间点达到该样本数，才报告样本级 Kruskal 筛选 p 值与 BH-FDR。下限为 2，理想为 3 或以上。'},
        {'key': 'min_cells_per_celltype_sample', 'label': '每样本/细胞类型最少细胞数', 'type': 'number', 'default': 20, 'step': 1, 'help': '低于此数的 sample × celltype 不纳入伪 bulk 表达趋势，减少稀有细胞的高方差影响。'},
        {'key': 'max_genes', 'label': '候选动态基因数', 'type': 'number', 'default': 300, 'step': 50, 'help': '优先使用 HVG 并按方差选择的最大基因数。基因 FDR 只在该候选集合内校正，属于探索性筛选而非全转录组检验。'},
        {'key': 'top_dynamic_genes', 'label': '动态基因热图 Top N', 'type': 'number', 'default': 30, 'step': 5, 'help': '按 FDR 和趋势效应展示的细胞类型特异动态基因数。'},
        {'key': 'top_celltypes', 'label': '组成曲线展示细胞类型数', 'type': 'number', 'default': 8, 'step': 1, 'help': '按总体平均比例选择用于折线图和热图的细胞类型数，避免图例过密。'},
        {'key': 'enable_gene_trends', 'label': '运行伪 bulk 基因动态', 'type': 'checkbox', 'default': True, 'help': '优先从 layers["counts"] 聚合 sample × celltype 原始 counts 并进行样本级时间筛选；缺少 counts 层时仅给描述性均值曲线。'},
        {'key': 'show_timepoint_umap', 'label': '生成时间点 UMAP', 'type': 'checkbox', 'default': True, 'help': '按真实时间点着色已有 UMAP，用于检查时间状态是否与嵌入结构一致。'},
        {'key': 'show_composition_trajectory', 'label': '生成细胞组成时间曲线', 'type': 'checkbox', 'default': True, 'help': '展示每个时间点的样本平均细胞类型比例；CSV 同时保留每个样本的原始比例。'},
        {'key': 'show_composition_heatmap', 'label': '生成细胞组成热图', 'type': 'checkbox', 'default': True, 'help': '展示时间点 × 细胞类型组成，适合发现短暂或阶段特异的细胞群变化。'},
        {'key': 'show_gene_heatmap', 'label': '生成动态基因热图', 'type': 'checkbox', 'default': True, 'help': '按每个 celltype 内的时间变化展示 Top 动态基因。行标准化仅用于可视化，不改变统计结果。'},
    ],
    'deg': [
        {'key': 'groupby', 'label': '分组依据', 'type': 'text', 'default': '', 'help': '差异分析的分类分组列名。留空或无效时自动使用 leiden；连续 QC 指标不能作为分组。'},
        {'key': 'reference', 'label': '参考组', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': 'rest', 'help': '参考组。rest 表示以所有其他组为对照。'},
        {'key': 'method', 'label': '统计方法', 'type': 'select', 'options': ['wilcoxon', 't-test', 'logreg', 't-test_overestim_var'], 'default': 'wilcoxon', 'help': '统计检验方法。'},
        {'key': 'n_genes', 'label': '显示 Top N 基因数', 'type': 'number', 'default': 20, 'help': '每个簇显示的 Top N 差异基因数。'},
        {'key': 'show_dotplot', 'label': '生成 DEG Dotplot', 'type': 'checkbox', 'default': True, 'help': '是否生成 Top 差异基因的 dotplot。'},
        {'key': 'plot_genes_umap', 'label': 'UMAP 展示基因（逗号分隔）', 'type': 'text', 'default': '', 'help': '指定要在 UMAP 上展示的基因名。'},
        {'key': 'custom_dotplot_genes', 'label': '自定义 Dotplot 基因（可选）', 'type': 'textarea', 'default': '', 'help': '自定义 Dotplot 基因列表。'},
        {'key': 'pval_cutoff', 'label': 'padj 过滤阈值', 'type': 'number', 'default': 0.05, 'step': 0.01, 'help': '校正后 p 值过滤阈值。'},
        {'key': 'logfc_cutoff', 'label': 'logFC 过滤阈值', 'type': 'number', 'default': 1.0, 'step': 0.1, 'help': 'log2 fold change 过滤阈值。'},
        {'key': 'min_pct', 'label': '最小表达比例', 'type': 'number', 'default': 0.1, 'step': 0.05, 'help': '基因在任一组中的最小表达细胞比例。'},
        {'key': 'correction_method', 'label': '多重检验校正', 'type': 'select', 'options': ['benjamini_hochberg', 'bonferroni', 'BY'], 'default': 'benjamini_hochberg', 'help': '多重检验校正方法。'},
        {'key': 'volcano_top_n', 'label': '火山图标注基因数', 'type': 'number', 'default': 10, 'help': '火山图上自动标注的 Top N 基因数。'},
        {'key': 'volcano_genes', 'label': '火山图自定义标注基因', 'type': 'textarea', 'default': '', 'help': '火山图上自定义标注的基因名。'},
        {'key': 'show_deg_counts_bar', 'label': '生成显著 DEG 数量图', 'type': 'checkbox', 'default': True, 'help': '按分组统计显著上调/下调基因数量，用于快速判断各 cluster 差异信号强弱。'},
        {'key': 'show_top_marker_umap_panel', 'label': '生成 Top marker UMAP 面板', 'type': 'checkbox', 'default': True, 'help': '自动选择各簇 Top marker 并生成表达 UMAP 面板，用于验证 marker 空间分布。'},
        {'key': 'top_marker_umap_genes', 'label': 'Top marker UMAP 基因数', 'type': 'number', 'default': 6, 'help': '自动 marker UMAP 面板中最多展示的基因数。'},
        {'key': 'show_marker_heatmap', 'label': '生成 Cluster marker 热图', 'type': 'checkbox', 'default': True, 'help': '展示每个簇 Top marker 在各簇中的平均表达 z-score，便于检查分群和注释一致性。'},
        {'key': 'marker_heatmap_top_n', 'label': '热图每簇 Top marker 数', 'type': 'number', 'default': 3, 'help': '每个簇纳入 marker 热图的 Top 基因数量。数值越大热图越全面但也越拥挤。'},
    ],
    'trajectory': [
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden', 'help': '用于轨迹推断和可视化的分类聚类列名；连续 QC 指标会被忽略或回退。'},
        {'key': 'plot_genes', 'label': '拟时序基因表达（可选）', 'type': 'textarea', 'default': '', 'help': '手动输入基因名，逗号或换行分隔。最多 10 个基因。'},
        {'key': 'enable_paga', 'label': '启用 PAGA', 'type': 'checkbox', 'default': False, 'help': '生成 PAGA 轨迹图，展示簇间连接强度。'},
        {'key': 'show_pseudotime_distribution', 'label': '生成拟时序分布图', 'type': 'checkbox', 'default': True, 'help': '按 cluster 展示 DPT pseudotime 分布，用于判断轨迹方向和分支是否合理。'},
        {'key': 'paga_threshold', 'label': 'PAGA 连接阈值', 'type': 'number', 'default': 0.05, 'step': 0.01, 'help': 'PAGA 连接强度阈值，低于此值的连接不显示。'},
        {'key': 'n_diffcomps', 'label': '扩散图成分', 'type': 'number', 'default': 15, 'help': '扩散图计算的成分数量。'},
        {'key': 'start_cluster', 'label': '起始簇（留空=自动）', 'type': 'text', 'default': '', 'help': '伪时间计算的起始簇名。留空则自动选择。'},
        {'key': 'n_dcs', 'label': '伪时间扩散成分', 'type': 'number', 'default': 10, 'help': '用于伪时间计算的扩散成分数量。'},
        {'key': 'n_branchings', 'label': '分支点数量', 'type': 'number', 'default': 0, 'help': '允许的分支点数量。0 为无分支（线性轨迹）。'},
    ],
    'proportion': [
        {'key': 'groupby', 'label': '分组依据', 'type': 'text', 'default': 'celltype', 'help': '统计比例的分类细胞类型列名；连续或高基数列会拒绝。'},
        {'key': 'batch_key', 'label': '展示/比较列（兼容旧流程）', 'type': 'text', 'default': 'batch', 'help': '用于汇总柱状图和旧版细胞计数关联检验的分类列。新的样本级比较优先使用下方“条件列”。'},
        {'key': 'sample_key', 'label': '生物学样本列（推荐）', 'type': 'text', 'default': '', 'help': '如 sample_id、library_id。样本而非单个细胞是比例统计的独立重复；缺失、近乎一细胞一个 ID 或技术 batch 时只输出描述性结果。'},
        {'key': 'condition_key', 'label': '条件列（推荐）', 'type': 'text', 'default': '', 'help': '如 condition、treatment、group。与 sample_key 同时填写后，按样本计算每个细胞类型比例并输出 Mann-Whitney/Kruskal 与 BH-FDR。'},
        {'key': 'analysis_unit', 'label': '统计单位', 'type': 'select', 'options': ['auto', 'sample', 'cell'], 'default': 'auto', 'help': 'auto：有可靠 sample_key + condition_key 时进行样本级检验，否则仅描述性；sample：强制请求样本级结果但设计不合格时仍降级；cell：只保留细胞计数关联展示，不报告样本级 p 值。'},
        {'key': 'min_samples_per_condition', 'label': '每条件最少样本数', 'type': 'number', 'default': 2, 'min': 2, 'step': 1, 'help': '每个条件至少有该数量的独立生物学样本才报告样本级 p 值和 BH-FDR；理想为 3 或以上。'},
        {'key': 'min_cells_per_sample', 'label': '每样本最少细胞数', 'type': 'number', 'default': 10, 'min': 1, 'step': 1, 'help': '总细胞数低于此阈值的样本仍写入表格，但不进入样本级检验，避免低深度样本驱动差异。'},
        {'key': 'confirm_batch_is_biological_sample', 'label': '确认 batch 实为生物学样本', 'type': 'checkbox', 'default': False, 'help': '仅当 sample_key 填 batch 且该列实际代表独立生物学样本时启用。'},
        {'key': 'confirm_batch_is_biological_condition', 'label': '确认 batch 实为生物学条件', 'type': 'checkbox', 'default': False, 'help': '仅当 condition_key 填 batch 且该列实际记录生物学分组时启用。'},
        {'key': 'compare_groups', 'label': '指定比较组（可选）', 'type': 'text', 'default': '', 'help': '格式：GroupA-vs-GroupB，多个比较用分号分隔。'},
        {'key': 'stat_test', 'label': '细胞计数关联检验', 'type': 'select', 'options': ['chi_square', 'fisher_exact', 'permutation'], 'default': 'chi_square', 'help': '用于旧版 cell count 列联表关联展示，不将细胞视为生物学重复。样本级推断自动选择二组 Mann-Whitney 或多组 Kruskal，并对 celltype 做 BH-FDR。'},
        {'key': 'n_permutations', 'label': '置换检验次数', 'type': 'number', 'default': 1000, 'help': '置换检验的置换次数。仅 permutation 方法生效。'},
        {'key': 'min_cells_per_group', 'label': '最小细胞数', 'type': 'number', 'default': 10, 'help': '每组最小细胞数，低于此值的组不参与比较。'},
        {'key': 'show_proportion_heatmap', 'label': '生成比例热图', 'type': 'checkbox', 'default': True, 'help': '以 heatmap 展示每个样本/分组中的细胞类型比例，便于横向比较组成差异。'},
    ],
    'cell_communication': [
        {'key': 'cluster_key', 'label': '细胞类型列', 'type': 'text', 'default': 'celltype', 'help': '用于通讯分析的分类细胞类型列名。需先运行注释模块；连续/高基数列会拒绝。'},
        {'key': 'resource', 'label': '配体-受体数据库', 'type': 'select', 'options': ['consensus', 'cellcall', 'cellchatdb', 'omnipath'], 'default': 'consensus', 'help': 'consensus：综合多数据库（推荐）。cellcall/cellchatdb：特定数据库。omnipath：OmniPath 数据库。'},
        {'key': 'organism', 'label': '物种', 'type': 'select', 'options': ['human', 'mouse'], 'default': 'human', 'help': '物种选择，影响配体-受体对匹配。'},
        {'key': 'min_prop', 'label': '最小表达比例', 'type': 'number', 'default': 0.1, 'step': 0.05, 'help': '基因在细胞群中的最小表达比例，低于此值的不参与分析。'},
        {'key': 'top_n_interactions', 'label': '展示 Top N', 'type': 'number', 'default': 20, 'help': '展示 Top N 个最强相互作用。'},
        {'key': 'show_heatmap', 'label': '生成通讯热图', 'type': 'checkbox', 'default': True, 'help': '生成细胞类型间通讯数量热图。'},
        {'key': 'show_network', 'label': '生成通讯网络图', 'type': 'checkbox', 'default': True, 'help': '以网络图展示 Top source-target 通讯关系，节点大小代表通讯连接度。'},
    ],
    'bulk_qc': [
        {'key': 'min_counts', 'label': '最小文库 reads 数', 'type': 'number', 'default': 100000, 'help': '最小文库 reads 数。低于此值的样本被过滤。人类/小鼠 RNA-seq 通常要求 ≥100000，小样本可降至 50000。'},
        {'key': 'min_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 5000, 'help': '每个样本检测到的最小基因数。低于此值的样本可能质量差。通常 5000-8000。'},
        {'key': 'max_mt_pct', 'label': '最大线粒体基因比例 (%)', 'type': 'number', 'default': 20.0, 'step': 0.1, 'help': '最大线粒体基因比例（%）。高于此值的样本可能降解严重。RNA-seq 通常 15-20%。'},
        {'key': 'max_ribo_pct', 'label': '最大核糖体基因比例 (%)', 'type': 'number', 'default': 40.0, 'step': 0.1, 'help': '最大核糖体基因比例（%）。RPL/RPS 基因比例过高提示 rRNA 污染。PolyA 建库通常 < 5-10%，rRNA 去除建库可至 40-50%。'},
        {'key': 'min_gini', 'label': '最小文库复杂度 (Gini)', 'type': 'number', 'default': 0, 'step': 0.01, 'help': '最小 Gini 系数（0 = 不过滤）。Gini > 0.8 提示文库复杂度低（PCR 过度扩增）。'},
        {'key': 'min_sample_expr', 'label': '基因最低表达样本数', 'type': 'number', 'default': 0, 'step': 1, 'help': '基因在至少 N 个样本中 count >= 阈值才保留。0 = 不过滤。建议设为样本总数的 10-20%。'},
        {'key': 'min_count_threshold', 'label': '基因表达 count 阈值', 'type': 'number', 'default': 1, 'step': 1, 'help': '基因 count 达到此值才算有效表达（配合"最低表达样本数"使用）。默认 1（即 count > 0）。'},
        {'key': 'group_column', 'label': '分组列名（可选）', 'type': 'text', 'default': '', 'help': '样本分组列名（adata.obs 中的列）。留空则自动从样本名推断（取第一个分隔符前的前缀）。填写后启用组内/组间距离分析和分组着色图。'},
        {'key': 'detect_outliers', 'label': '检测离群样本', 'type': 'checkbox', 'default': True, 'help': '基于 PCA 马氏距离检测离群样本。仅在 summary 中告警，不自动剔除。'},
        {'key': 'filter_strategy', 'label': '过滤策略', 'type': 'select', 'options': ['standard', 'strict', 'custom'], 'default': 'standard', 'help': 'standard：推荐阈值；strict：严格阈值（适合大样本高质量数据）；custom：自定义所有阈值。'},
    ],
    'bulk_normalize': [
        {'key': 'method', 'label': '标准化方法', 'type': 'select', 'options': ['deseq2', 'tmm', 'cpm', 'vst', 'rlog', 'log2', 'log2_quantile'], 'default': 'deseq2',
         'help': '原始整数 count：DESeq2 中位比率（推荐）、TMM-CPM（组成偏差明显）、CPM、近似 VST/rlog（用于 PCA/热图，不是 DESeq2 原版变换）。FPKM/TPM：仅 log2(x+1)；log2_quantile 会强制样本分布一致，仅在该假设成立时使用。平台会拒绝将 FPKM/TPM 用于 count 方法。'},
        {'key': 'min_expr_value', 'label': '最小表达阈值 (CPM)', 'type': 'number', 'default': 1, 'step': 0.1,
         'help': '原始 count 时为 CPM 阈值；FPKM/TPM 的 log2 模式下为原始 FPKM/TPM 阈值。默认 1。'},
        {'key': 'min_expr_samples', 'label': '最小表达样本数', 'type': 'number', 'default': 3, 'step': 1,
         'help': '基因在至少 N 个样本中达到最小表达阈值才保留。0 = 不过滤。建议设为最小组的样本数。'},
        {'key': 'max_zero_pct', 'label': '最大零值比例 (%)', 'type': 'number', 'default': 0, 'step': 1,
         'help': '基因在超过此比例的样本中为零则被过滤。0 = 不过滤。建议 50-70%。'},
        {'key': 'cpm_target', 'label': 'CPM 缩放目标', 'type': 'number', 'default': 1000000, 'step': 100000,
         'help': 'CPM 标准化的缩放目标值。默认 1e6（标准 CPM）。'},
    ],
    'bulk_deg': [
        {'key': 'groupby', 'label': '分组列名', 'type': 'text', 'default': '', 'help': '分类分组列名。adata.obs 中用于区分实验组和对照组的列，如 condition、treatment、group；不能填连续 QC 指标。'},
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
    ],
    'bulk_pca': [
        {'key': 'n_comps', 'label': 'PCA 主成分数量', 'type': 'number', 'default': 10, 'help': 'PCA 主成分数量。通常 5-10 即可。样本数少时自动降至 n_samples-1。'},
        {'key': 'color_by', 'label': '颜色分组列名', 'type': 'text', 'default': '_auto_group_', 'help': '用于 PCA 着色。默认从样本名自动识别分组；也可从下方 obs 列选择真实 metadata（如 condition、batch）。无法识别时会明确提示，不会静默合并为同一组。'},
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
         'options': ['auto', 'all', 'none'], 'default': 'auto',
         'help': 'auto=样本较多时自动间隔显示（推荐）；all=显示全部；none=隐藏。'},
        {'key': 'gene_font_size', 'label': '基因名字体大小', 'type': 'number', 'default': 8, 'step': 1,
         'help': '基因名标签字体大小。'},
        {'key': 'sample_font_size', 'label': '样本名字体大小', 'type': 'number', 'default': 9, 'step': 1,
         'help': '样本名标签字体大小。'},
        {'key': 'annotation_columns', 'label': '注释条列名', 'type': 'text', 'default': '',
         'help': '额外注释条列名，逗号分隔（如 group,batch）。groupby 列自动包含。'},
        {'key': 'groupby', 'label': '样本分组列名', 'type': 'text', 'default': '',
         'help': '主分组列名，用于默认注释条。留空则不添加。'},
        {'key': 'sample_display_mode', 'label': '热图展示样本范围', 'type': 'select',
         'options': ['all', 'deg_groups', 'selected_groups', 'selected_samples'], 'default': 'all',
         'option_labels': {'all': '全部样本', 'deg_groups': '仅 DEG 两组',
                           'selected_groups': '仅勾选分组', 'selected_samples': '仅勾选样本'},
         'help': 'all=全部样本；deg_groups=只显示所选 DEG 比较的两组；selected_groups=仅显示勾选分组；selected_samples=仅显示指定样本。仅影响热图展示，不重新计算 DEG。'},
        {'key': 'sample_display_groups', 'label': '展示分组', 'type': 'dynamic_multiselect',
         'depends_on': 'groupby', 'default': '', 'show_if': {'sample_display_mode': 'selected_groups'},
         'help': '从样本分组列中勾选要展示的组。'},
        {'key': 'sample_display_names', 'label': '展示样本名', 'type': 'textarea', 'default': '',
         'show_if': {'sample_display_mode': 'selected_samples'},
         'help': '输入需要展示的样本名，逗号或换行分隔；样本名必须与输入矩阵完全一致。'},
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
        {'key': 'input_source', 'label': 'DEG 比较结果', 'type': 'text', 'default': '',
         'help': '选择一个已完成的 DEG 比较。ORA/GSEA 每次仅分析一个比较，结果表需包含 gene 和 regulation/log2FC 列。'},
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
        {'key': 'upset_strict', 'label': 'Upset 严格模式', 'type': 'checkbox', 'default': True,
         'help': '严格模式：仅显示"仅属于该组合"的基因。关闭则显示"至少属于该组合"的基因。'},
        {'key': 'exclude_mixed', 'label': '排除方向不一致基因', 'type': 'checkbox', 'default': False,
         'help': '一致性评分中排除方向不一致（Mixed）的基因。'},
        {'key': 'filter_expression', 'label': '筛选表达式（可选）', 'type': 'textarea', 'default': '',
         'help': '通过集合逻辑表达式筛选目标基因集。支持 AND/OR/NOT/XOR 运算符和 ALL/ANY/ONLY 简写。'
                 '示例：ALL:up | hmc3-vs-ctrl:up AND rapa-vs-ctrl:down | ONLY[hmc3-vs-ctrl]:up'},
        {'key': 'upset_top_n', 'label': 'Upset 图显示数', 'type': 'number', 'default': 20, 'step': 5,
         'help': 'Upset 图显示的交集模式数量。'},
        {'key': 'logfc_clip_range', 'label': 'logFC 截断范围', 'type': 'number', 'default': 5.0, 'step': 0.5,
         'help': 'logFC 热图的颜色截断范围（±值）。'},
        {'key': 'filter_show_n', 'label': '筛选热图显示基因数', 'type': 'number', 'default': 80, 'step': 10,
         'help': '筛选基因 logFC 热图最多显示的基因数。'},
    ],
}


def _condition_value_matches(actual, expected):
    if isinstance(expected, (list, tuple, set)):
        return any(_condition_value_matches(actual, item) for item in expected)
    if isinstance(expected, bool):
        if isinstance(actual, str):
            return (actual.lower() in {'1', 'true', 'yes', 'on'}) is expected
        return bool(actual) is expected
    return str(actual) == str(expected)


def param_is_active(param, values):
    """Return whether a schema parameter should be active for current form values."""
    conditions = param.get('show_if')
    if not conditions:
        return True
    return all(_condition_value_matches(values.get(key), expected)
               for key, expected in conditions.items())


def filter_active_params(schema, params):
    """Drop params hidden by show_if for the current parameter values."""
    values = {param['key']: param.get('default') for param in schema}
    values.update(params or {})
    active_keys = {
        param['key']
        for param in schema
        if param_is_active(param, values)
    }
    return {key: value for key, value in (params or {}).items() if key in active_keys}


def parse_form_params(schema, form):
    """从 Flask request.form 中解析参数，根据 schema 定义做类型转换。"""
    def form_values(key):
        if hasattr(form, 'getlist'):
            return form.getlist(key)
        value = form.get(key)
        return [] if value in (None, '') else [value]

    values = {}
    for param in schema:
        key = param['key']
        if param['type'] == 'checkbox':
            values[key] = form.get(key) == 'on'
        elif param['type'] == 'dynamic_multiselect':
            values[key] = ','.join(form_values(key))
        else:
            val = form.get(key)
            values[key] = val if val not in (None, '') else param['default']

    params = {}
    for param in schema:
        if not param_is_active(param, values):
            continue
        val = form.get(param['key'])
        if param['type'] == 'number':
            params[param['key']] = float(val) if val else param['default']
        elif param['type'] == 'checkbox':
            params[param['key']] = form.get(param['key']) == 'on'
        elif param['type'] == 'dynamic_multiselect':
            params[param['key']] = ','.join(form_values(param['key']))
        else:
            params[param['key']] = val or param['default']
    return filter_active_params(schema, params)


def list_upload_files(pid, is_bulk):
    """列出项目的上传文件。"""
    from config import Config
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
    return uploaded_files
