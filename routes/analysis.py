import os
import json
from flask import Blueprint, render_template, request, redirect, url_for, flash
from models import Project, AnalysisTask
from config import Config

analysis_bp = Blueprint('analysis', __name__)

MODULE_LIST = [
    {'name': 'qc', 'display': '质控', 'desc': '过滤细胞，去除双细胞'},
    {'name': 'preprocess', 'display': '预处理', 'desc': '标准化，选择高变异基因'},
    {'name': 'dimred', 'display': '降维分析', 'desc': 'PCA, UMAP'},
    {'name': 'batch_correct', 'display': '批次校正', 'desc': 'Harmony, ComBat, SysVI'},
    {'name': 'clustering', 'display': '聚类分析', 'desc': 'Leiden 聚类'},
    {'name': 'annotation', 'display': '细胞注释', 'desc': '基于 Marker 的细胞类型注释'},
    {'name': 'deg', 'display': '差异表达', 'desc': '差异表达基因分析'},
    {'name': 'trajectory', 'display': '轨迹分析', 'desc': '拟时序分析'},
    {'name': 'proportion', 'display': '比例分析', 'desc': '细胞比例分析'},
    {'name': 'bulk_qc', 'display': 'Bulk 质控', 'desc': '文库大小、基因检测、离群值过滤'},
    {'name': 'bulk_normalize', 'display': 'Bulk 标准化', 'desc': 'DESeq2 / CPM / 分位数标准化'},
    {'name': 'bulk_deg', 'display': 'Bulk 差异表达', 'desc': '组间差异基因检测（火山图、MA图）'},
    {'name': 'bulk_pca', 'display': 'Bulk PCA/UMAP', 'desc': 'PCA 和 UMAP 降维可视化'},
    {'name': 'bulk_heatmap', 'display': 'Bulk 热图', 'desc': 'Top 差异基因热图、样本相关性热图'},
]

PARAM_SCHEMAS = {
    'qc': [
        {'key': 'mito_perc', 'label': '最大线粒体比例', 'type': 'number', 'default': 0.2, 'step': 0.01},
        {'key': 'nUMIs', 'label': '最小 UMI 数', 'type': 'number', 'default': 500},
        {'key': 'detected_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 250},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch'},
    ],
    'preprocess': [
        {'key': 'n_top_genes', 'label': '高变异基因数量', 'type': 'number', 'default': 2000},
        {'key': 'target_sum', 'label': '标准化目标总数', 'type': 'number', 'default': 10000},
    ],
    'dimred': [
        {'key': 'n_comps', 'label': 'PCA 主成分数量', 'type': 'number', 'default': 50},
        {'key': 'use_mde', 'label': '使用 MDE（加速 UMAP）', 'type': 'checkbox', 'default': False},
    ],
    'batch_correct': [
        {'key': 'method', 'label': '校正方法', 'type': 'select', 'options': ['harmony', 'combat', 'sysvi'], 'default': 'harmony'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch'},
        {'key': 'n_pcs', 'label': '主成分数量', 'type': 'number', 'default': 50},
        {'key': 'max_epochs', 'label': '最大迭代轮数（仅 SysVI）', 'type': 'number', 'default': 200},
    ],
    'clustering': [
        {'key': 'resolutions', 'label': 'Leiden 分辨率（逗号分隔）', 'type': 'text', 'default': '0.6,0.8,1.0'},
        {'key': 'n_neighbors', 'label': '邻居数量', 'type': 'number', 'default': 15},
    ],
    'annotation': [
        {'key': 'method', 'label': '注释方法', 'type': 'select', 'options': ['auto_marker', 'manual'], 'default': 'auto_marker'},
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden'},
        {'key': 'resolution', 'label': 'Leiden 分辨率', 'type': 'text', 'default': '0.8'},
    ],
    'deg': [
        {'key': 'groupby', 'label': '分组依据', 'type': 'text', 'default': 'celltype'},
        {'key': 'reference', 'label': '参考组', 'type': 'text', 'default': 'rest'},
        {'key': 'method', 'label': '统计方法', 'type': 'select', 'options': ['wilcoxon', 't-test', 'logreg'], 'default': 'wilcoxon'},
        {'key': 'n_genes', 'label': '显示 Top N 基因数', 'type': 'number', 'default': 20},
    ],
    'trajectory': [
        {'key': 'method', 'label': '轨迹方法', 'type': 'select', 'options': ['diffusion_map', 'slingshot'], 'default': 'diffusion_map'},
        {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden'},
    ],
    'proportion': [
        {'key': 'groupby', 'label': '分组依据', 'type': 'text', 'default': 'celltype'},
        {'key': 'batch_key', 'label': '批次列名', 'type': 'text', 'default': 'batch'},
    ],
    'bulk_qc': [
        {'key': 'min_counts', 'label': '最小文库 reads 数', 'type': 'number', 'default': 100000},
        {'key': 'min_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 5000},
        {'key': 'max_mt_pct', 'label': '最大线粒体基因比例 (%)', 'type': 'number', 'default': 20.0, 'step': 0.1},
    ],
    'bulk_normalize': [
        {'key': 'method', 'label': '标准化方法', 'type': 'select', 'options': ['deseq2', 'cpm', 'log2_quantile'], 'default': 'deseq2'},
    ],
    'bulk_deg': [
        {'key': 'groupby', 'label': '分组列名', 'type': 'text', 'default': 'condition'},
        {'key': 'group1', 'label': '实验组名称', 'type': 'text', 'default': ''},
        {'key': 'group2', 'label': '对照组名称', 'type': 'text', 'default': ''},
        {'key': 'method', 'label': '统计方法', 'type': 'select', 'options': ['t-test', 'mann-whitney'], 'default': 't-test'},
        {'key': 'fc_threshold', 'label': 'Fold Change 阈值', 'type': 'number', 'default': 2.0, 'step': 0.1},
        {'key': 'pval_threshold', 'label': 'padj 显著性阈值', 'type': 'number', 'default': 0.05, 'step': 0.01},
        {'key': 'top_n', 'label': 'Top N 差异基因数', 'type': 'number', 'default': 20},
    ],
    'bulk_pca': [
        {'key': 'n_comps', 'label': 'PCA 主成分数量', 'type': 'number', 'default': 10},
        {'key': 'color_by', 'label': '颜色分组列名（留空则不着色）', 'type': 'text', 'default': ''},
        {'key': 'run_umap', 'label': '同时运行 UMAP', 'type': 'checkbox', 'default': True},
    ],
    'bulk_heatmap': [
        {'key': 'heatmap_type', 'label': '热图类型', 'type': 'select', 'options': ['top_var', 'deg'], 'default': 'top_var'},
        {'key': 'top_n', 'label': '显示基因数', 'type': 'number', 'default': 50},
        {'key': 'groupby', 'label': '样本分组列名（可选）', 'type': 'text', 'default': ''},
    ],
}

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
    uploaded_h5ad = []
    if os.path.isdir(uploads_dir):
        for f in os.listdir(uploads_dir):
            if f.endswith('.h5ad'):
                uploaded_h5ad.append({'name': f, 'path': os.path.join(uploads_dir, f)})
    return render_template('analysis_select.html', project=p, module=mod_info,
                          schema=schema, completed_tasks=completed_tasks,
                          uploaded_h5ad=uploaded_h5ad, all_modules=MODULE_LIST)
