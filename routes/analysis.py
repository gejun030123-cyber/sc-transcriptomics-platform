import os
import json
from flask import Blueprint, render_template, request, redirect, url_for, flash
from models import Project, AnalysisTask
from config import Config

analysis_bp = Blueprint('analysis', __name__)

SC_MODULE_LIST = [
    {'name': 'qc', 'display': '质控', 'desc': '过滤细胞，去除双细胞'},
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
        {'key': 'marker_set', 'label': 'Marker 基因集', 'type': 'select', 'options': ['TME', 'Immune'], 'default': 'TME',
         'help': 'TME：肿瘤微环境 marker（上皮、CAF、内皮、免疫细胞等）。Immune：免疫细胞 marker（T、B、NK、髓系等）。'},
        {'key': 'custom_markers', 'label': '自定义 Marker（可选）', 'type': 'text', 'default': '',
         'help': '格式：CellType1:GENE1,GENE2;CellType2:GENE3,GENE4。留空则使用上方选择的预设基因集。'},
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
    ],
    'trajectory': [
        {'key': 'method', 'label': '轨迹方法', 'type': 'select', 'options': ['diffusion_map', 'slingshot'], 'default': 'diffusion_map', 'help': '轨迹推断方法。diffusion_map：基于扩散图的拟时序，适合连续过渡。slingshot：基于 MST 的轨迹，适合分支结构。'},
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden', 'help': '用于轨迹推断的聚类列名。'},
    ],
    'proportion': [
        {'key': 'groupby', 'label': '分组依据', 'type': 'text', 'default': 'celltype', 'help': '统计比例的细胞类型列名。通常为 celltype 或 leiden。'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch', 'help': '用于比较的分组列名。如 batch、condition、treatment 等。'},
    ],
    'bulk_qc': [
        {'key': 'min_counts', 'label': '最小文库 reads 数', 'type': 'number', 'default': 100000, 'help': '最小文库 reads 数。低于此值的样本被过滤。人类/小鼠 RNA-seq 通常要求 ≥100000，小样本可降至 50000。'},
        {'key': 'min_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 5000, 'help': '每个样本检测到的最小基因数。低于此值的样本可能质量差。通常 5000-8000。'},
        {'key': 'max_mt_pct', 'label': '最大线粒体基因比例 (%)', 'type': 'number', 'default': 20.0, 'step': 0.1, 'help': '最大线粒体基因比例（%）。高于此值的样本可能降解严重。RNA-seq 通常 15-20%。'},
    ],
    'bulk_normalize': [
        {'key': 'method', 'label': '标准化方法', 'type': 'select', 'options': ['deseq2', 'cpm', 'log2_quantile'], 'default': 'deseq2', 'help': '标准化方法。DESeq2：中位比率法，适用于差异分析前标准化，RNA-seq 金标准。CPM：每百万计数，简单但不考虑组成偏差。log2_quantile：分位数标准化，适合样本间可比性要求高的场景。'},
    ],
    'bulk_deg': [
        {'key': 'groupby', 'label': '分组列名', 'type': 'text', 'default': '', 'help': '分组列名。adata.obs 中用于区分实验组和对照组的列。如 condition、treatment、group。'},
        {'key': 'group1', 'label': '实验组', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': '', 'help': '实验组名称。将与对照组比较计算差异基因。'},
        {'key': 'group2', 'label': '对照组', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': '', 'help': '对照组名称。rest 表示以所有其他样本为对照。'},
        {'key': 'method', 'label': '统计方法', 'type': 'select', 'options': ['t-test', 'mann-whitney', 'deseq2'], 'default': 't-test', 'help': '统计方法。t-test：参数检验，适合正态分布数据，速度快。Mann-Whitney：非参数检验，不假设正态分布，更稳健。DESeq2：基于负二项分布的差异分析，RNA-seq 金标准，需要原始计数。'},
        {'key': 'fc_threshold', 'label': 'Fold Change 阈值', 'type': 'number', 'default': 2.0, 'step': 0.1, 'help': 'Fold Change 阈值。log2FC > log2(fc) 为上调，< -log2(fc) 为下调。常用值：1.5（宽松）、2.0（标准）、4.0（严格）。'},
        {'key': 'pval_threshold', 'label': 'padj 显著性阈值', 'type': 'number', 'default': 0.05, 'step': 0.01, 'help': '调整后 p-value 显著性阈值。0.05 为标准，0.01 为严格，0.1 为宽松探索性分析。'},
        {'key': 'top_n', 'label': 'Top N 差异基因数', 'type': 'number', 'default': 20, 'help': '结果中展示的 Top N 差异基因数。用于火山图标注和 Top 基因 CSV 导出。'},
        {'key': 'base_mean_filter', 'label': '最低平均表达量', 'type': 'number', 'default': 1, 'step': 0.5, 'help': '过滤低表达基因。BaseMean 低于此值的基因不参与分析和绘图。建议 1-10。'},
        {'key': 'plot_genes', 'label': '额外展示基因（逗号分隔，可选）', 'type': 'text', 'default': '', 'help': '指定要额外绘制箱线图的基因名，多个用逗号分隔。Top 差异基因会自动生成箱线图，此字段用于补充其他感兴趣的基因。'},
    ],
    'bulk_pca': [
        {'key': 'n_comps', 'label': 'PCA 主成分数量', 'type': 'number', 'default': 10, 'help': 'PCA 主成分数量。通常 5-10 即可。样本数少时自动降至 n_samples-1。'},
        {'key': 'color_by', 'label': '颜色分组列名（留空则不着色）', 'type': 'text', 'default': '', 'help': '用于着色的 obs 列名。留空则不着色。如 condition、batch、celltype 等。'},
        {'key': 'dimred_method', 'label': '降维方法', 'type': 'select', 'options': ['pca', 'umap', 'tsne'], 'default': 'pca', 'help': '降维可视化方法。PCA：线性降维，保留全局结构。UMAP：非线性降维，保留局部结构。t-SNE：非线性降维，适合发现聚类。'},
    ],
    'bulk_heatmap': [
        {'key': 'heatmap_type', 'label': '热图类型', 'type': 'select', 'options': ['top_var', 'deg'], 'default': 'top_var', 'help': '热图类型。top_var：显示最高变异的基因。deg：显示差异表达基因（需先运行 DEG 分析）。'},
        {'key': 'top_n', 'label': '显示基因数', 'type': 'number', 'default': 50, 'help': '热图中显示的基因数量。通常 30-100。过多会导致热图难以阅读。'},
        {'key': 'groupby', 'label': '样本分组列名（可选）', 'type': 'text', 'default': '', 'help': '样本分组列名，用于在热图旁添加分组注释条。留空则不添加。'},
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
