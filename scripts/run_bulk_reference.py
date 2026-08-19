#!/usr/bin/env python3
"""Bulk RNA-seq Reference Runner

用指定参考表达矩阵跑通 Bulk RNA-seq 流水线：
  bulk_qc → bulk_normalize → bulk_pca → bulk_deg → bulk_heatmap
条件允许时（--run-optional）额外执行：bulk_enrichment、bulk_timecourse、bulk_deg_integration。

复用现有模块实例（MODULE_REGISTRY），不重写分析逻辑。

用法示例：
  # Dataset A
  python scripts/run_bulk_reference.py \
    --input /home/oelab/data/GJ/all.fpkm_anno.xls \
    --dataset-label dataset_a --project-name bulk_dataset_a \
    --groupby _auto_group

  # Dataset B
  python scripts/run_bulk_reference.py \
    --input /home/oelab/data/GJ/all.genes.expression.anno.xls \
    --dataset-label dataset_b --project-name bulk_dataset_b \
    --groupby _auto_group --run-optional
"""

import argparse
import glob
import json
import os
import shutil
import sys
import time
import traceback
from datetime import datetime

# 项目根目录加入 sys.path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from config import Config
from modules import MODULE_REGISTRY
from modules.schemas import PARAM_SCHEMAS

# 模块中文名映射
MODULE_DISPLAY = {
    'bulk_qc': 'Bulk RNA-seq 质控',
    'bulk_normalize': '数据标准化',
    'bulk_pca': 'PCA / UMAP 降维',
    'bulk_deg': '差异表达分析',
    'bulk_heatmap': '热图可视化',
    'bulk_enrichment': '通路富集分析',
    'bulk_timecourse': '时序分析',
    'bulk_deg_integration': '多组差异整合分析',
}

# 输入文件名 → 默认比较
DEFAULT_COMPARISONS = {
    'all.fpkm_anno.xls': 'hmc3-vs-ctrl;moclel-vs-ctrl;rapa-vs-ctrl',
    'all.genes.expression.anno.xls': 'NH4Cl-vs-Ctr;PEA-vs-Ctr;TMAO-vs-Ctr',
}


def build_default_params(module_name):
    """从 PARAM_SCHEMAS 读取模块默认参数。"""
    schema = PARAM_SCHEMAS.get(module_name, [])
    params = {}
    for field in schema:
        key = field['key']
        default = field.get('default')
        if field.get('type') == 'select' and 'options' in field:
            opts = field['options']
            if default not in opts and opts:
                default = opts[0]
        params[key] = default
    return params


def make_progress_callback(module_name):
    """生成带模块名前缀的进度回调。"""
    def cb(pct, msg):
        print(f"  [{module_name}] {pct:3d}% {msg}")
    return cb


def run_step(module_name, params, project_dir, input_path):
    """执行单个模块，返回 result dict。"""
    cls = MODULE_REGISTRY.get(module_name)
    if not cls:
        raise ValueError(f"未知模块: {module_name}")

    progress_cb = make_progress_callback(module_name)
    module = cls(project_dir=project_dir, params=params, progress_callback=progress_cb)
    return module.run(input_path)


def find_deg_csv(project_dir):
    """查找 bulk_deg 输出的 CSV 文件。"""
    results_dir = os.path.join(project_dir, 'results')
    patterns = ['bulk_deg_results*.csv', 'bulk_deg_*.csv', '*deg*.csv']
    for pat in patterns:
        matches = sorted(glob.glob(os.path.join(results_dir, pat)))
        if matches:
            return matches[0]
    return None


def generate_chinese_report(project_dir, report, input_file, comparisons):
    """生成中文运行报告。"""
    lines = []
    lines.append('# Bulk RNA-seq 参考运行报告\n')
    lines.append(f'**生成时间**: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
    lines.append(f'**输入文件**: `{input_file}`\n')
    lines.append(f'**项目目录**: `{project_dir}`\n')
    lines.append(f'**比较组**: `{comparisons}`\n')
    lines.append('')

    lines.append('## 模块执行结果\n')
    lines.append('| 模块 | 状态 | 输出 h5ad | CSV | Plotly JSON |')
    lines.append('|------|------|-----------|-----|-------------|')

    for step in report['steps']:
        name = step['module']
        display = MODULE_DISPLAY.get(name, name)
        status = '✅ 完成' if step['status'] == 'completed' else '❌ 失败'
        h5ad = step.get('output_adata_path', '—')
        if h5ad:
            h5ad = os.path.basename(h5ad)
        else:
            h5ad = '—'
        csv_count = sum(1 for rf in step.get('result_files', []) if rf.get('file_type') == 'csv')
        plotly_count = sum(1 for rf in step.get('result_files', []) if rf.get('file_type') == 'plotly_json')
        lines.append(f'| {display} | {status} | {h5ad} | {csv_count} | {plotly_count} |')

    lines.append('')

    # 失败模块详情
    failed_steps = [s for s in report['steps'] if s['status'] == 'failed']
    if failed_steps:
        lines.append('## 失败模块详情\n')
        for step in failed_steps:
            name = step['module']
            display = MODULE_DISPLAY.get(name, name)
            lines.append(f'### {display}\n')
            lines.append(f'**状态**: 失败\n')
            if step.get('error'):
                # 截取关键 traceback
                tb_lines = step['error'].strip().split('\n')
                key_lines = tb_lines[-5:] if len(tb_lines) > 5 else tb_lines
                lines.append('**关键错误**:\n')
                lines.append('```')
                lines.append('\n'.join(key_lines))
                lines.append('```\n')
            lines.append('')

    # 摘要统计
    lines.append('## 摘要统计\n')
    for step in report['steps']:
        if step.get('summary'):
            name = step['module']
            display = MODULE_DISPLAY.get(name, name)
            summary = step['summary']
            lines.append(f'### {display}\n')
            for key in ['n_cells', 'n_samples', 'n_genes', 'n_comparisons',
                        'n_degs', 'n_up', 'n_down', 'method', 'n_groups']:
                if key in summary:
                    lines.append(f'- **{key}**: {summary[key]}')
            lines.append('')

    # 下一步建议
    lines.append('## 下一步建议\n')
    if report['overall_status'] == 'completed':
        lines.append('- 核心流程全部完成，可检查各模块输出文件')
        lines.append('- 可使用 `--run-optional` 运行额外分析（enrichment、timecourse、deg_integration）')
    else:
        lines.append('- 修复失败模块后重新运行')
        for step in failed_steps:
            display = MODULE_DISPLAY.get(step['module'], step['module'])
            lines.append(f'- 检查 {display} 的输入数据和参数')
    lines.append('')

    report_path = os.path.join(project_dir, '中文运行报告.md')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    return report_path


def main():
    parser = argparse.ArgumentParser(
        description='Bulk RNA-seq Reference Runner — 跑通完整 Bulk 流水线',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--input', required=True, help='输入表达矩阵路径（CSV/TSV/h5ad/xls）')
    parser.add_argument('--dataset-label', default='', help='数据集标签（如 dataset_a）')
    parser.add_argument('--project-name', default='bulk_reference_run', help='项目目录名')
    parser.add_argument('--groupby', default='_auto_group', help='分组列名（默认 _auto_group）')
    parser.add_argument('--group1', default='', help='实验组名称')
    parser.add_argument('--group2', default='', help='对照组名称')
    parser.add_argument('--comparisons', default='', help='多组比较，格式：A-vs-B;C-vs-D')
    parser.add_argument('--time-column', default='', help='时间列名（有则跑 bulk_timecourse）')
    parser.add_argument(
        '--output-dir',
        default=os.environ.get(
            'BULK_REFERENCE_OUTPUT_DIR',
            os.path.join(Config.DATA_DIR, 'bulk_reference_output'),
        ),
        help='输出根目录（默认写入 DATA_DIR，可用 BULK_REFERENCE_OUTPUT_DIR 覆盖）',
    )
    parser.add_argument('--method', default='t-test', help='DEG 统计方法（默认 t-test）')
    parser.add_argument('--normalize-method', default='deseq2', help='归一化方法（默认 deseq2）')
    parser.add_argument('--fc-threshold', type=float, default=2.0, help='DEG Fold Change 阈值（默认 2.0）')
    parser.add_argument('--padj-threshold', type=float, default=0.05, help='DEG FDR/padj 阈值（默认 0.05）')
    parser.add_argument('--top-n', type=int, default=50, help='DEG/热图展示基因数（默认 50）')
    parser.add_argument('--min-expr-value', type=float, default=1.0, help='标准化前最小表达阈值（默认 1）')
    parser.add_argument('--min-expr-samples', type=int, default=3, help='基因至少在 N 个样本中达到表达阈值（默认 3）')
    parser.add_argument('--organism', default='Human', help='物种（默认 Human）')
    parser.add_argument('--run-optional', action='store_true', help='运行可选步骤（enrichment/timecourse/deg_integration）')
    parser.add_argument('--heatmap-source', default='top_var', choices=['top_var', 'deg'],
                        help='热图基因来源（默认 top_var，设为 deg 使用差异基因）')

    args = parser.parse_args()

    # --- 校验输入文件 ---
    if not os.path.isfile(args.input):
        print(f"❌ 输入文件不存在: {args.input}", file=sys.stderr)
        sys.exit(1)

    # --- 自动推断比较 ---
    comparisons_str = args.comparisons
    if not comparisons_str:
        input_basename = os.path.basename(args.input)
        if input_basename in DEFAULT_COMPARISONS:
            comparisons_str = DEFAULT_COMPARISONS[input_basename]
            print(f"📋 自动推断比较: {comparisons_str}")
        elif args.group1 and args.group2:
            comparisons_str = f"{args.group1}-vs-{args.group2}"
        else:
            print("❌ 未提供 --comparisons 且无法从文件名推断，请手动指定", file=sys.stderr)
            sys.exit(1)

    comparisons_list = [c.strip() for c in comparisons_str.split(';') if c.strip()]
    if args.group1 and args.group2:
        pair = f"{args.group1}-vs-{args.group2}"
        if pair not in comparisons_list:
            comparisons_list.append(pair)
            comparisons_str = ';'.join(comparisons_list)

    # --- 创建项目目录 ---
    project_dir = os.path.join(args.output_dir, args.project_name)
    for subdir in ['uploads', 'intermediate', 'results', 'plots']:
        os.makedirs(os.path.join(project_dir, subdir), exist_ok=True)

    # 复制输入文件到 uploads/
    input_filename = os.path.basename(args.input)
    upload_path = os.path.join(project_dir, 'uploads', input_filename)
    if os.path.abspath(args.input) != os.path.abspath(upload_path):
        shutil.copy2(args.input, upload_path)
        print(f"✅ 输入文件已复制到: {upload_path}")
    else:
        print(f"✅ 输入文件已在 uploads 目录: {upload_path}")

    # --- 初始化报告 ---
    report = {
        'project_name': args.project_name,
        'dataset_label': args.dataset_label or args.project_name,
        'input_file': args.input,
        'groupby': args.groupby,
        'comparisons': comparisons_str,
        'time_column': args.time_column,
        'started_at': datetime.now().isoformat(),
        'steps': [],
    }

    current_input = upload_path
    core_failed = False

    # --- 定义流水线步骤 ---
    core_steps = ['bulk_qc', 'bulk_normalize', 'bulk_pca', 'bulk_deg', 'bulk_heatmap']
    optional_steps = []
    if args.run_optional:
        optional_steps.append('bulk_deg_integration')
        optional_steps.append('bulk_enrichment')
        if args.time_column:
            optional_steps.append('bulk_timecourse')

    all_steps = core_steps + optional_steps

    for step_name in all_steps:
        print(f"\n{'='*60}")
        print(f"▶ 步骤: {step_name} ({MODULE_DISPLAY.get(step_name, step_name)})")
        print(f"  输入: {current_input}")
        print(f"{'='*60}")

        # 构建参数
        params = build_default_params(step_name)
        params['groupby'] = args.groupby

        if step_name == 'bulk_qc':
            params['group_column'] = args.groupby

        elif step_name == 'bulk_deg':
            params['comparisons'] = comparisons_str
            if args.group1 and args.group2:
                params['group1'] = args.group1
                params['group2'] = args.group2
            params['method'] = args.method
            params['fc_threshold'] = args.fc_threshold
            params['pval_threshold'] = args.padj_threshold
            params['top_n'] = args.top_n

        elif step_name == 'bulk_normalize':
            params['method'] = args.normalize_method
            params['min_expr_value'] = args.min_expr_value
            params['min_expr_samples'] = args.min_expr_samples

        elif step_name == 'bulk_pca':
            params['color_by'] = args.groupby

        elif step_name == 'bulk_heatmap':
            params['groupby'] = args.groupby
            params['gene_import_source'] = args.heatmap_source
            params['heatmap_type'] = args.heatmap_source
            params['top_n'] = args.top_n

        elif step_name == 'bulk_enrichment':
            params['organism'] = args.organism
            params['pvalue_cutoff'] = args.padj_threshold
            params['top_n'] = min(args.top_n, 30)
            params['split_direction'] = True
            # 自动查找 DEG CSV 作为 input_source
            deg_csv = find_deg_csv(project_dir)
            if deg_csv:
                params['input_source'] = deg_csv
                print(f"  📂 使用 DEG CSV: {deg_csv}")
            else:
                print(f"  ⚠️ 未找到 DEG CSV，跳过 enrichment")
                report['steps'].append({
                    'module': step_name,
                    'display': MODULE_DISPLAY.get(step_name, step_name),
                    'status': 'skipped',
                    'started_at': datetime.now().isoformat(),
                    'finished_at': datetime.now().isoformat(),
                    'summary': {},
                    'result_files': [],
                    'output_adata_path': None,
                    'error': '跳过: 未找到 DEG CSV',
                })
                continue

        elif step_name == 'bulk_timecourse':
            params['time_column'] = args.time_column
            params['group_column'] = args.groupby

        elif step_name == 'bulk_deg_integration':
            # 传字符串而非 list
            # 留空让模块自动使用所有 DEG CSV（selected_comparisons 期望文件 key，非比较名）
            params['selected_comparisons'] = ''
            params['fc_threshold'] = args.fc_threshold
            params['pval_threshold'] = args.padj_threshold
            params['consistency_n'] = args.top_n

        step_record = {
            'module': step_name,
            'display': MODULE_DISPLAY.get(step_name, step_name),
            'status': 'pending',
            'started_at': None,
            'finished_at': None,
            'summary': None,
            'params': params,
            'result_files': [],
            'output_adata_path': None,
            'error': None,
        }

        t0 = time.time()
        step_record['started_at'] = datetime.now().isoformat()

        try:
            result = run_step(step_name, params, project_dir, current_input)
            elapsed = time.time() - t0

            step_record['status'] = 'completed'
            step_record['finished_at'] = datetime.now().isoformat()
            step_record['summary'] = result.get('summary', {})
            step_record['result_files'] = [
                {
                    'file_path': rf.get('file_path', ''),
                    'file_type': rf.get('file_type', ''),
                    'category': rf.get('category', ''),
                    'label': rf.get('label', ''),
                }
                for rf in (result.get('result_files') or [])
            ]
            step_record['output_adata_path'] = result.get('output_adata')

            if result.get('output_adata'):
                current_input = result['output_adata']

            if step_name == 'bulk_deg':
                comparison_names = (result.get('summary') or {}).get('comparisons', [])
                if comparison_names:
                    label_path = os.path.join(project_dir, 'results', 'bulk_deg_comparison_labels.json')
                    label_map = {
                        f'bulk_deg_results_{idx}.csv': name
                        for idx, name in enumerate(comparison_names)
                    }
                    with open(label_path, 'w', encoding='utf-8') as f:
                        json.dump(label_map, f, ensure_ascii=False, indent=2)

            print(f"  ✅ {step_name} 完成 ({elapsed:.1f}s)")
            summary = result.get('summary', {})
            if summary:
                for key in ['n_cells', 'n_samples', 'n_genes', 'n_comparisons',
                            'n_degs', 'n_up', 'n_down', 'method']:
                    if key in summary:
                        print(f"     {key}: {summary[key]}")

        except Exception as e:
            elapsed = time.time() - t0
            tb = traceback.format_exc()
            step_record['status'] = 'failed'
            step_record['finished_at'] = datetime.now().isoformat()
            step_record['error'] = tb

            print(f"  ❌ {step_name} 失败 ({elapsed:.1f}s)")
            print(f"  错误:\n{tb}")

            if step_name in core_steps:
                core_failed = True
            else:
                print(f"  ⚠️ 可选步骤 {step_name} 失败，继续执行后续步骤")

        report['steps'].append(step_record)

    # --- 写入 JSON 报告 ---
    report['finished_at'] = datetime.now().isoformat()
    report['overall_status'] = 'failed' if core_failed else 'completed'

    report_path = os.path.join(project_dir, 'bulk_reference_report.json')
    with open(report_path, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    # --- 生成中文报告 ---
    chinese_report_path = generate_chinese_report(project_dir, report, args.input, comparisons_str)

    # --- 打印总结 ---
    print(f"\n{'='*60}")
    print("📋 流水线执行总结")
    print(f"{'='*60}")
    for step in report['steps']:
        icon = '✅' if step['status'] == 'completed' else ('⏭️' if step['status'] == 'skipped' else '❌')
        print(f"  {icon} {step['display']}: {step['status']}")
    print(f"\nJSON 报告: {report_path}")
    print(f"中文报告: {chinese_report_path}")
    print(f"项目目录: {project_dir}")

    if core_failed:
        print("\n❌ 核心流程因步骤失败而中止")
        sys.exit(1)
    else:
        print("\n✅ 核心流程全部完成")


if __name__ == '__main__':
    main()
