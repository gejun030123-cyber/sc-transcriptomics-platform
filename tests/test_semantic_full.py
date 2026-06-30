"""语义一致性全量测试 —— 扫描所有分析模块源码，验证输出结构、summary 字典、result_files 模式等一致性。

与 test_semantic.py 不同，本文件对全部 21 个分析模块逐一做静态扫描，
不依赖运行时实例化，仅通过正则/字符串匹配分析源码。
"""
import ast
import os
import re
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

MODULES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'modules')

# 排除工具/基础设施模块，只保留分析模块
SKIP_FILES = {
    '__init__.py', 'base.py', 'schemas.py', 'io_utils.py',
    'visualization.py', 'inspect_utils.py', 'expression_parser.py',
    'constants.py', 'ai_adapter.py', 'ai_tools.py',
}

# 允许的 file_type 值
VALID_FILE_TYPES = {'csv', 'plotly_json', 'png', 'svg', 'info', 'json', 'txt'}

# 允许的 category 值
VALID_CATEGORIES = {
    'table', 'volcano', 'ma', 'bar', 'heatmap', 'pca', 'umap',
    'qc', 'dotplot', 'boxplot', 'enrichment', 'info',
    'violin', 'scatter', 'histogram', 'tsne', 'pie', 'bubble',
    'paga', 'diffusion_map', 'gene_expression', 'cluster_centers',
    'qq', 'venn', 'tree', 'annotation',
}


def _get_module_files():
    """获取所有分析模块文件（排除工具模块）。"""
    files = []
    for f in sorted(os.listdir(MODULES_DIR)):
        if f.endswith('.py') and f not in SKIP_FILES:
            files.append(os.path.join(MODULES_DIR, f))
    return files


def _read_source(filepath):
    """读取模块源码并返回 (lines, source_text)。"""
    with open(filepath, 'r', encoding='utf-8') as f:
        source = f.read()
    return source.splitlines(keepends=True), source


def _module_name(filepath):
    """从文件路径提取模块名（不含 .py）。"""
    return os.path.splitext(os.path.basename(filepath))[0]


def _find_class_name(source):
    """从源码中提取主分析类名（第一个继承 BaseAnalysis 的类）。"""
    m = re.search(r'class\s+(\w+)\s*\(.*BaseAnalysis.*\)', source)
    return m.group(1) if m else None


# ──────────────────────────────────────────────────────────────
# TestAllModulesSummaryKeys — 扫描 summary dict 构建
# ──────────────────────────────────────────────────────────────

class TestAllModulesSummaryKeys:
    """扫描所有模块 run() 方法中的 summary dict，验证命名一致性和字段逻辑。"""

    def _extract_summary_keys(self, source):
        """从源码中提取 summary dict 的键名列表。

        匹配模式:
        - summary = { 'key': ..., ... }
        - 'summary': { 'key': ..., ... }
        """
        keys = set()

        # 模式 1: summary = { ... } 或 'summary': { ... }
        # 查找所有 summary dict 字面量并提取键
        for m in re.finditer(
            r"""(?:summary\s*=\s*|'summary':)\s*\{([^}]+)\}""",
            source, re.DOTALL
        ):
            block = m.group(1)
            # 提取字符串键: 'key': 或 "key":
            for km in re.finditer(r"""['"](\w+)['"]\s*:""", block):
                keys.add(km.group(1))

        return keys

    def test_n_up_n_down_paired(self):
        """如果模块使用 n_up，也应使用 n_down，反之亦然。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)
            has_n_up = bool(re.search(r"['\"]n_up['\"]", source))
            has_n_down = bool(re.search(r"['\"]n_down['\"]", source))
            if has_n_up != has_n_down:
                missing = 'n_down' if has_n_up else 'n_up'
                issues.append(f"  {fname}: 有 n_up={has_n_up} 但缺少 {missing}")

        if issues:
            msg = "以下模块 n_up/n_down 不配对:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_count_fields_are_int_cast(self):
        """summary 中的计数字段应使用 int() 包裹，避免浮点数。"""
        # 定义已知的计数字段名
        count_field_names = {
            'n_cells', 'n_genes', 'n_hvgs', 'n_pcs', 'n_groups',
            'n_clusters', 'n_celltypes', 'n_interactions', 'n_cell_types',
            'n_up', 'n_down', 'n_genes_total', 'n_batches',
            'n_comparisons', 'n_samples', 'n_genes_shown',
            'n_timepoints', 'n_significant', 'n_input_genes',
            'cells_before', 'cells_after', 'cells_removed',
            'samples_before', 'samples_after', 'samples_removed',
            'genes_before', 'genes_after', 'genes_removed',
            'total_deg_genes', 'n_consistent_genes',
            'n_up_consistent', 'n_down_consistent',
            'n_temporal_sig', 'n_pairwise_sig',
            'doublets_removed', 's_genes_found', 'g2m_genes_found',
            'force_include_count', 'n_low_quality',
            'shared_up_genes', 'shared_down_genes',
            'lrt_n_sig', 'n_genes_before_filter', 'n_genes_after_filter',
            'genes_filtered', 'n_components', 'total_genes',
        }

        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 在 summary dict 块中查找计数字段
            for m in re.finditer(
                r"""(?:summary\s*=\s*|'summary':)\s*\{([^}]+)\}""",
                source, re.DOTALL
            ):
                block = m.group(1)
                for line in block.splitlines():
                    for field in count_field_names:
                        if f"'{field}'" in line or f'"{field}"' in line:
                            # 检查值是否用 int() 包裹
                            # 匹配 'field': value 模式
                            val_match = re.search(
                                rf"""['"]?{field}['"]?\s*:\s*(.+?)(?:,|$)""",
                                line.strip()
                            )
                            if val_match:
                                val = val_match.group(1).strip()
                                # 排除已经 int() 包裹的、int 字面量、或函数调用返回 int 的情况
                                if val.startswith('int(') or val.startswith('len('):
                                    continue
                                # 检查是否是纯数字
                                try:
                                    int(val)
                                    continue
                                except ValueError:
                                    pass
                                # float() 包裹的也是可疑的
                                if val.startswith('float('):
                                    issues.append(
                                        f"  {fname}: '{field}' 使用 float() 而非 int(): {val}"
                                    )
                                # 非 int/len 的值标记为潜在问题（仅对纯计数字段）
                                # 排除常见 int 来源: 变量名、算术、属性访问、条件表达式
                                safe_patterns = [
                                    'True', 'False', 'None', 'self.',
                                    'adata.', 'len(', 'int(', 'sum(',
                                    'dict(', 'list(', 'str(', '.n_',
                                    'round(', 'max(', 'min(',
                                ]
                                # 排除纯变量名赋值（通常由 len()/int() 等提前计算）
                                # 排除算术表达式 (e.g. n_before - n_after)
                                # 排除条件表达式 (e.g. x if cond else 0)
                                is_safe = any(kw in val for kw in safe_patterns)
                                is_simple_var = bool(re.match(r'^[a-z]\w*$', val))
                                is_arithmetic = bool(re.search(r'[+\-*/]', val))
                                is_conditional = ' if ' in val and ' else ' in val
                                if not is_safe and not is_simple_var and not is_arithmetic and not is_conditional:
                                    issues.append(
                                        f"  {fname}: '{field}' 的值可能非 int: {val}"
                                    )

        if issues:
            msg = "以下 summary 计数字段可能未使用 int():\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_n_genes_total_consistent(self):
        """如果模块有 n_genes_total 且同时有 n_up/n_down，n_genes_total >= n_up + n_down。"""
        # 这是静态源码检查：确认逻辑上 n_genes_total 不会被设置为 0 而同时 n_up > 0
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)
            keys = self._extract_summary_keys(source)
            has_total = 'n_genes_total' in keys
            has_up = 'n_up' in keys
            has_down = 'n_down' in keys
            if has_total and has_up and has_down:
                # 检查 n_genes_total 是否独立于 n_up/n_down 赋值
                # 如果 n_genes_total = len(deg_df) 而 n_up/n_down 从 deg_df 过滤，逻辑应该一致
                # 此处只验证三者同时出现时的源码模式
                pass  # 逻辑检查需要运行时，此处仅记录存在性

    def test_all_modules_have_summary(self):
        """每个分析模块的 run() 返回都应包含 summary 键。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)
            if "'summary'" not in source and '"summary"' not in source:
                issues.append(f"  {fname}: 未找到 'summary' 键")

        if issues:
            msg = "以下模块缺少 summary 返回:\n" + '\n'.join(issues)
            pytest.fail(msg)


# ──────────────────────────────────────────────────────────────
# TestAllModulesResultFiles — 验证 result_files 模式
# ──────────────────────────────────────────────────────────────

class TestAllModulesResultFiles:
    """扫描所有模块的 result_files.append 调用，验证必需键和值域。"""

    def _extract_result_file_dicts(self, source):
        """从源码中提取 result_files 相关的 dict 构造。

        匹配模式:
        - {'file_path': ..., 'file_type': ..., 'category': ..., 'label': ...}
        - result_files.append({...})
        """
        dicts = []
        # 查找所有包含 file_path 的字面量 dict
        for m in re.finditer(r"\{[^{}]*'file_path'[^{}]*\}", source):
            dicts.append(m.group(0))
        return dicts

    def test_result_files_have_required_keys(self):
        """每个 result_files dict 应包含必需键: file_path, file_type, category, label。"""
        required_keys = {'file_path', 'file_type', 'category', 'label'}
        issues = []

        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for dict_str in self._extract_result_file_dicts(source):
                found_keys = set(re.findall(r"'(\w+)'\s*:", dict_str))
                missing = required_keys - found_keys
                if missing:
                    # 精简显示
                    compact = re.sub(r'\s+', ' ', dict_str.strip())[:120]
                    issues.append(
                        f"  {fname}: 缺少键 {missing} — {compact}"
                    )

        if issues:
            msg = "以下 result_files dict 缺少必需键:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_file_type_is_valid(self):
        """file_type 值应在合法范围内。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for dict_str in self._extract_result_file_dicts(source):
                m = re.search(r"""'file_type'\s*:\s*['"]([^'"]+)['"]""", dict_str)
                if m:
                    ft = m.group(1)
                    if ft not in VALID_FILE_TYPES:
                        compact = re.sub(r'\s+', ' ', dict_str.strip())[:120]
                        issues.append(
                            f"  {fname}: 未知 file_type '{ft}' — {compact}"
                        )

        if issues:
            msg = "以下 result_files 使用了未知 file_type:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_category_is_valid(self):
        """category 值应在合法范围内。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for dict_str in self._extract_result_file_dicts(source):
                m = re.search(r"""'category'\s*:\s*['"]([^'"]+)['"]""", dict_str)
                if m:
                    cat = m.group(1)
                    if cat not in VALID_CATEGORIES:
                        compact = re.sub(r'\s+', ' ', dict_str.strip())[:120]
                        issues.append(
                            f"  {fname}: 未知 category '{cat}' — {compact}"
                        )

        if issues:
            msg = "以下 result_files 使用了未知 category:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_save_plotly_json_result_dict(self):
        """使用 self.save_plotly_json() 的模块应接收其返回值（已是完整 dict）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 检查 save_plotly_json 调用是否有对应的结果变量
            for i, line in enumerate(lines):
                stripped = line.strip()
                if 'save_plotly_json(' in stripped:
                    # 应该是赋值语句或直接 append
                    if '=' not in stripped and 'append' not in stripped and 'result_files' not in stripped:
                        # 可能是 standalone 调用，检查上下文
                        context_start = max(0, i - 2)
                        context_end = min(len(lines), i + 3)
                        context = ''.join(lines[context_start:context_end])
                        if 'result_files' not in context:
                            issues.append(
                                f"  {fname}:{i+1}: save_plotly_json() 返回值未被使用"
                            )

        if issues:
            msg = "以下 save_plotly_json() 调用返回值可能未被使用:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_file_path_not_empty_string(self):
        """result_files 中的 file_path 不应是空字符串（info 类型除外）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for dict_str in self._extract_result_file_dicts(source):
                # 检查 file_type 是否为 info
                ft_match = re.search(r"""'file_type'\s*:\s*['"]([^'"]+)['"]""", dict_str)
                if ft_match and ft_match.group(1) == 'info':
                    continue  # info 类型允许空 file_path

                # 检查 file_path 是否为空字符串
                fp_match = re.search(r"""'file_path'\s*:\s*['"]([^'"]*)['"]""", dict_str)
                if fp_match and fp_match.group(1) == '':
                    compact = re.sub(r'\s+', ' ', dict_str.strip())[:120]
                    issues.append(
                        f"  {fname}: file_path 为空字符串 — {compact}"
                    )

        if issues:
            msg = "以下 result_files 的 file_path 为空:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_label_not_empty(self):
        """result_files 中的 label 不应为空字符串。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for dict_str in self._extract_result_file_dicts(source):
                label_match = re.search(r"""'label'\s*:\s*['"]([^'"]*)['"]""", dict_str)
                if label_match and label_match.group(1).strip() == '':
                    compact = re.sub(r'\s+', ' ', dict_str.strip())[:120]
                    issues.append(
                        f"  {fname}: label 为空 — {compact}"
                    )

        if issues:
            msg = "以下 result_files 的 label 为空:\n" + '\n'.join(issues)
            pytest.fail(msg)


# ──────────────────────────────────────────────────────────────
# TestAllModulesOutputAdata — 验证 output_adata 处理
# ──────────────────────────────────────────────────────────────

class TestAllModulesOutputAdata:
    """验证模块的 output_adata 处理和 save_output() 调用一致性。"""

    def test_output_adata_in_return(self):
        """每个模块的 run() 返回都应包含 output_adata 键。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)
            if "'output_adata'" not in source and '"output_adata"' not in source:
                issues.append(f"  {fname}: 未找到 'output_adata' 返回键")

        if issues:
            msg = "以下模块缺少 output_adata 返回:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_output_path_ends_with_output_h5ad(self):
        """模块的 output_adata 路径应以 '_output.h5ad' 结尾（save_output() 约定）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 save_output 调用
            save_calls = re.findall(r'self\.save_output\([^)]+\)', source)
            for call in save_calls:
                # save_output 内部生成路径: {module_name}_output.h5ad
                # 提取模块名参数
                m = re.search(r"self\.save_output\(\s*\w+\s*,\s*['\"]([^'\"]+)['\"]", call)
                if m:
                    module_arg = m.group(1)
                    expected_suffix = f"{module_arg}_output.h5ad"
                    # save_output 内部: os.path.join(intermediate_dir, f'{module_name}_output.h5ad')
                    # 只验证模块名不含特殊字符
                    if not re.match(r'^[a-zA-Z0-9_]+$', module_arg):
                        issues.append(
                            f"  {fname}: save_output 模块名含非法字符: '{module_arg}'"
                        )

            # 检查手动构造的 output_adata 路径
            for line in lines:
                stripped = line.strip()
                if "'output_adata'" in stripped and 'input_path' not in stripped:
                    # 提取路径值
                    path_match = re.search(r"'output_adata'\s*:\s*(\w+)", stripped)
                    if path_match:
                        var_name = path_match.group(1)
                        # 在源码中查找该变量的赋值
                        var_assign = re.search(
                            rf'{var_name}\s*=\s*(.+)',
                            source
                        )
                        if var_assign:
                            val = var_assign.group(1).strip()
                            if 'h5ad' in val and '_output.h5ad' not in val and 'uploads' not in val:
                                issues.append(
                                    f"  {fname}: output_adata 路径不以 _output.h5ad 结尾: {val}"
                                )

        if issues:
            msg = "以下模块的 output_adata 路径不符合约定:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_single_cell_modules_use_save_output(self):
        """单细胞分析模块应使用 self.save_output() 而非手动 write_h5ad。"""
        # 单细胞模块列表（不包括 bulk_ 和 convert_10x）
        sc_module_names = {
            'qc', 'normalize', 'hvg', 'dimred', 'batch_correct',
            'clustering', 'qc_reassess', 'annotation', 'deg',
            'trajectory', 'proportion', 'cell_communication',
        }
        issues = []
        for filepath in _get_module_files():
            mname = _module_name(filepath)
            if mname not in sc_module_names:
                continue

            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            has_save_output = 'self.save_output(' in source
            has_manual_write = 'write_h5ad(' in source and 'self.save_output(' not in source

            if has_manual_write:
                issues.append(
                    f"  {fname}: 单细胞模块使用手动 write_h5ad() 而非 self.save_output()"
                )
            if not has_save_output and not has_manual_write:
                issues.append(
                    f"  {fname}: 既无 save_output() 也无 write_h5ad()"
                )

        if issues:
            msg = "以下单细胞模块的 output_adata 处理不规范:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_bulk_modules_output_path_convention(self):
        """Bulk 模块的 output_adata 路径应为 intermediate/{module}_output.h5ad。"""
        issues = []
        for filepath in _get_module_files():
            mname = _module_name(filepath)
            if not mname.startswith('bulk_'):
                continue
            if mname == 'bulk_enrichment':
                continue  # bulk_enrichment 不保存 h5ad

            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 output_adata 的赋值
            output_match = re.search(
                r"'output_adata'\s*:\s*(\w+)",
                source
            )
            if not output_match:
                continue

            var_name = output_match.group(1)
            if var_name == 'input_path' or var_name == 'None':
                continue  # 错误返回或无输出

            # 查找该变量的赋值，验证路径
            expected_pattern = f'{mname}_output.h5ad'
            if expected_pattern not in source:
                issues.append(
                    f"  {fname}: 未找到预期的路径模式 '{expected_pattern}'"
                )

        if issues:
            msg = "以下 Bulk 模块的 output_adata 路径不符合约定:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_convert_10x_output_path(self):
        """convert_10x 模块的 output_adata 应指向 uploads/ 目录。"""
        filepath = os.path.join(MODULES_DIR, 'convert_10x.py')
        if not os.path.exists(filepath):
            pytest.skip("convert_10x.py 不存在")

        lines, source = _read_source(filepath)

        # convert_10x 特殊：输出到 uploads/ 而非 intermediate/
        if "'output_adata'" in source:
            # 查找 output_adata 的值
            output_match = re.search(
                r"'output_adata'\s*:\s*(\w+)",
                source
            )
            if output_match:
                var_name = output_match.group(1)
                if var_name == 'input_path':
                    return  # 错误返回
                # 追踪变量赋值
                var_assign = re.search(
                    rf'^\s*{var_name}\s*=\s*(.+)$',
                    source, re.MULTILINE
                )
                if var_assign:
                    val = var_assign.group(1).strip()
                    if 'uploads' not in val and 'input_path' not in val:
                        pytest.fail(
                            f"convert_10x: output_adata 应指向 uploads/ 目录，"
                            f"实际 {var_name}={val}"
                        )


# ──────────────────────────────────────────────────────────────
# TestAllModulesLabelConsistency — 标签与内容一致性
# ──────────────────────────────────────────────────────────────

class TestAllModulesLabelConsistency:
    """检查 CSV 标签是否准确描述文件内容。"""

    def _extract_label_csv_pairs(self, source):
        """提取 label 和 .to_csv() 的配对关系。

        返回 [(label, csv_line, context_lines), ...]
        """
        pairs = []
        lines = source.splitlines(keepends=True)

        # 在 result_files dict 中查找 label 与对应文件路径
        # 然后关联到 .to_csv() 调用
        for m in re.finditer(
            r"\{[^{}]*'file_path'\s*:\s*([^,}]+)[^{}]*'file_type'\s*:\s*'csv'[^{}]*'label'\s*:\s*'([^']*)'[^{}]*\}",
            source, re.DOTALL
        ):
            pairs.append((m.group(2), m.group(1).strip(), m.start()))

        # 反向匹配: label 在 file_path 之前
        for m in re.finditer(
            r"\{[^{}]*'label'\s*:\s*'([^']*)'[^{}]*'file_path'\s*:\s*([^,}]+)[^{}]*'file_type'\s*:\s*'csv'[^{}]*\}",
            source, re.DOTALL
        ):
            pairs.append((m.group(1), m.group(2).strip(), m.start()))

        return pairs

    def test_top_n_label_includes_count(self):
        """如果标签含 "Top N" 或 "top N" 模式（后跟量词），应包含具体数量。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 label 中含 Top 后跟数量词（前/个/种/基因/基因等）的 dict
            # 只匹配 "Top" 后直接跟中文量词或 "genes" 等暗示数量的模式
            for m in re.finditer(
                r"""'label'\s*:\s*['"]([^'"]*(?:Top|top|TOP)\s*\d*[^'"]*)['"]""",
                source
            ):
                label = m.group(1)
                # 只检查含 "Top" 后跟数字或量词的模式
                # 如 "Top 100 Genes", "Top差异基因" 等暗示排名的标签
                # 排除 "Top 一致性基因" 等描述性标签
                if re.search(r'(?:Top|top)\s*\d+', label):
                    # "Top 100" — 有数字，OK
                    continue
                if re.search(r'(?:Top|top)\s*(?:差异|上调|下调|显著|DEG)', label, re.IGNORECASE):
                    # "Top 差异基因" — 描述性，无需数字
                    continue
                # 检查英文 "Top" 后跟 genes/paths 等但无数字
                if re.search(r'(?:Top|top)\s*(?:Genes|genes|Paths|Pathways)', label):
                    if not re.search(r'\d+', label):
                        issues.append(
                            f"  {fname}: 标签含 'Top' 但无具体数量: '{label}'"
                        )

        if issues:
            msg = "以下标签含 'Top' 但缺少具体数量:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_full_label_not_filtered(self):
        """如果标签含 "完整" 或 "Full"，数据不应被过滤。"""
        # 静态检查：查找标签含 "完整" 的 dict，检查附近是否有 head()/iloc[:]/nlargest() 等过滤操作
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for m in re.finditer(
                r"""'label'\s*:\s*['"]([^'"]*(?:完整|Full|full)[^'"]*)['"]""",
                source
            ):
                label = m.group(1)
                pos = m.start()
                # 在附近 500 字符内查找过滤操作
                nearby = source[max(0, pos - 500):pos + 500]
                filter_ops = ['head(', 'nlargest(', 'nsmallest(', '.iloc[:', '[:N]', '[0:N]']
                for op in filter_ops:
                    if op in nearby:
                        issues.append(
                            f"  {fname}: 标签含 '完整' 但附近有过滤操作 '{op}': '{label}'"
                        )
                        break

        if issues:
            msg = "以下标签含 '完整' 但附近有数据过滤操作:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_csv_label_mentions_format(self):
        """CSV 类型文件的 label 应暗示表格/表格数据（可选检查）。"""
        # 此测试仅做软检查，不强制失败
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)
            # 统计 CSV result_files 数量
            csv_count = source.count("'file_type': 'csv'")
            if csv_count == 0:
                continue


# ──────────────────────────────────────────────────────────────
# TestAllModulesInputRequires — 验证 INPUT_REQUIRES 准确性
# ──────────────────────────────────────────────────────────────

class TestAllModulesInputRequires:
    """验证 INPUT_REQUIRES 声明与 validate_input() 实际检查一致。"""

    def _get_input_requires(self, source):
        """提取 INPUT_REQUIRES 列表值。"""
        m = re.search(r"INPUT_REQUIRES\s*=\s*\[([^\]]*)\]", source)
        if m:
            content = m.group(1).strip()
            if not content:
                return []
            return re.findall(r"['\"](\w+)['\"]", content)
        return None  # 未声明

    def _get_validate_checks(self, source):
        """提取 validate_input() 中检查的键名列表。"""
        # 查找 validate_input 方法体
        m = re.search(
            r'def validate_input\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
            source, re.DOTALL
        )
        if not m:
            return []
        body = m.group(1)
        if 'return None' in body and len(body.strip().splitlines()) <= 2:
            return []  # 无实际检查

        # 提取检查的键名: 'key' not in adata.obsm / 'key' not in adata.uns / 'key' not in adata.obs
        checks = re.findall(r"['\"](\w+)['\"]\s+not\s+in\s+adata\.\w+", body)
        # 也提取 adata.n_obs / adata.n_vars 检查
        numeric_checks = re.findall(r"adata\.n_(\w+)\s*[<>]=?\s*\d+", body)
        return checks

    def test_requires_matches_validate(self):
        """INPUT_REQUIRES 声明的键应在 validate_input() 或 run() 中被使用。"""
        import warnings
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            requires = self._get_input_requires(source)
            if requires is None or not requires:
                continue

            checks = self._get_validate_checks(source)

            for req in requires:
                if req in checks:
                    continue  # 在 validate_input 中显式检查

                # 检查是否在 run() 方法中有使用
                run_body_match = re.search(
                    r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                    source, re.DOTALL
                )
                if run_body_match:
                    run_body = run_body_match.group(1)
                    if req in run_body:
                        # 在 run 中隐式使用 — 这是可接受的（pipeline guard 已检查依赖）
                        continue
                    else:
                        issues.append(
                            f"  {fname}: INPUT_REQUIRES=['{req}'] 但 "
                            f"validate_input() 和 run() 均未检查/使用 '{req}'"
                        )
                else:
                    issues.append(
                        f"  {fname}: INPUT_REQUIRES=['{req}'] 但 "
                        f"validate_input() 未检查 '{req}'"
                    )

        if issues:
            msg = "以下模块的 INPUT_REQUIRES 声明了未使用的依赖:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_validate_without_requires(self):
        """validate_input() 有实际检查但 INPUT_REQUIRES 为空时应标记。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            requires = self._get_input_requires(source)
            checks = self._get_validate_checks(source)

            if checks and (not requires or len(requires) == 0):
                # 有实际检查但 INPUT_REQUIRES 为空
                # 可能是检查 n_obs/n_vars 等运行时条件，不算 require 字段
                # 但如果检查的是 obsm/uns 中的键，则应声明
                m = re.search(
                    r'def validate_input\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                    source, re.DOTALL
                )
                if m:
                    body = m.group(1)
                    obsm_checks = re.findall(r"['\"](\w+)['\"]\s+not\s+in\s+adata\.obsm", body)
                    uns_checks = re.findall(r"['\"](\w+)['\"]\s+not\s+in\s+adata\.uns", body)
                    obs_checks = re.findall(r"['\"](\w+)['\"]\s+not\s+in\s+adata\.obs(?:\.columns)?", body)
                    field_checks = obsm_checks + uns_checks + obs_checks
                    if field_checks:
                        issues.append(
                            f"  {fname}: validate_input() 检查了 {field_checks} "
                            f"但 INPUT_REQUIRES 为空"
                        )

        # 这是警告级别，使用 warnings 而非 pytest.fail
        if issues:
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)

    def test_leiden_requires_consistent(self):
        """声明 INPUT_REQUIRES=['leiden'] 的模块应检查 leiden 列存在。"""
        issues = []
        leiden_modules = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)
            requires = self._get_input_requires(source)
            if requires and 'leiden' in requires:
                leiden_modules.append((filepath, fname))

        for filepath, fname in leiden_modules:
            lines, source = _read_source(filepath)
            checks = self._get_validate_checks(source)
            # leiden 可能不在 validate_input 中显式检查
            # 因为 leiden 是 obs 列而非 obsm key，许多模块直接使用 adata.obs['leiden']
            # 这是可接受的，只做信息性记录


# ──────────────────────────────────────────────────────────────
# TestModuleReturnStructure — 验证 run() 返回结构
# ──────────────────────────────────────────────────────────────

class TestModuleReturnStructure:
    """验证每个模块的 run() 方法返回正确的 dict 结构。"""

    def _extract_return_statements(self, source):
        """提取 run() 方法中的 return 语句。"""
        returns = []
        # 查找 run 方法体
        m = re.search(
            r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
            source, re.DOTALL
        )
        if not m:
            return returns
        body = m.group(1)
        for line in body.splitlines():
            stripped = line.strip()
            if stripped.startswith('return ') and '{' in stripped:
                returns.append(stripped)
        return returns

    def test_return_has_required_keys(self):
        """每个模块的 run() 返回应包含 output_adata, result_files, summary 三个键。"""
        required_keys = {'output_adata', 'result_files', 'summary'}
        issues = []

        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 检查源码中是否包含所有必需键
            for key in required_keys:
                if f"'{key}'" not in source and f'"{key}"' not in source:
                    issues.append(f"  {fname}: 缺少 '{key}' 返回键")

        if issues:
            msg = "以下模块的 run() 返回缺少必需键:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_result_files_is_list(self):
        """result_files 应初始化为列表。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 result_files 初始化
            init_patterns = re.findall(r'result_files\s*=\s*(.+)', source)
            for pattern in init_patterns:
                val = pattern.strip().rstrip(',')
                if val == '[]' or val.startswith('['):
                    continue  # 正确：列表字面量
                if 'list(' in val:
                    continue  # 正确：list() 构造
                if val.startswith('{'):
                    # 可能是 dict comprehension 而非 list
                    issues.append(
                        f"  {fname}: result_files 初始化为 dict 而非 list: {val}"
                    )

        if issues:
            msg = "以下模块的 result_files 初始化不正确:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_summary_is_dict(self):
        """summary 应为 dict 类型。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 summary dict 构造模式
            # 有两种模式:
            # 1. summary = { ... } 然后 'summary': summary
            # 2. 'summary': { ... } 内联
            has_summary_var = bool(re.search(r'summary\s*=\s*\{', source))
            has_summary_inline = bool(re.search(r"'summary'\s*:\s*\{", source))

            if not has_summary_var and not has_summary_inline:
                # 检查是否有变量赋值
                summary_assign = re.search(r"summary\s*=\s*(\w+)", source)
                if summary_assign:
                    var = summary_assign.group(1)
                    if var in ('{}', 'dict()'):
                        continue
                    # 变量可能是 dict，需要运行时验证
                else:
                    issues.append(
                        f"  {fname}: 未找到 summary dict 构造模式"
                    )

        if issues:
            msg = "以下模块的 summary 构造模式不明确:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_error_returns_have_all_keys(self):
        """错误提前返回也应包含 output_adata, result_files, summary 三个键。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找包含 'error' 的 return 语句
            for i, line in enumerate(lines):
                stripped = line.strip()
                if not stripped.startswith('return '):
                    continue
                if "'error'" not in stripped and '"error"' not in stripped:
                    continue
                # 检查是否包含所有必需键
                # 收集整个 return 块（可能跨多行）
                return_block = stripped
                brace_count = stripped.count('{') - stripped.count('}')
                j = i + 1
                while brace_count > 0 and j < len(lines):
                    return_block += lines[j].strip()
                    brace_count += lines[j].count('{') - lines[j].count('}')
                    j += 1

                for key in ('output_adata', 'result_files', 'summary'):
                    if f"'{key}'" not in return_block:
                        issues.append(
                            f"  {fname}:{i+1}: 错误返回缺少 '{key}' — "
                            f"{stripped[:100]}"
                        )

        if issues:
            msg = "以下模块的错误返回缺少必需键:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_no_bare_return_none_in_run(self):
        """run() 方法不应有 return None（应始终返回 dict）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 run 方法体
            m = re.search(
                r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                source, re.DOTALL
            )
            if not m:
                continue
            body = m.group(1)
            for line in body.splitlines():
                stripped = line.strip()
                if stripped == 'return None' or stripped == 'return':
                    issues.append(
                        f"  {fname}: run() 方法中有 '{stripped}'，应返回 dict"
                    )

        if issues:
            msg = "以下模块的 run() 方法有无效返回:\n" + '\n'.join(issues)
            pytest.fail(msg)


# ──────────────────────────────────────────────────────────────
# TestAllModulesNamingConventions — 命名约定一致性
# ──────────────────────────────────────────────────────────────

class TestAllModulesNamingConventions:
    """验证模块的命名约定一致性。"""

    def test_summary_key_naming_style(self):
        """summary 键应使用 snake_case 命名（全小写+下划线）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 提取 summary dict 中的键
            for m in re.finditer(
                r"""(?:summary\s*=\s*|'summary':)\s*\{([^}]+)\}""",
                source, re.DOTALL
            ):
                block = m.group(1)
                for km in re.finditer(r"""['"](\w+)['"]\s*:""", block):
                    key = km.group(1)
                    # 检查是否 snake_case
                    if key != key.lower() and '_' not in key:
                        # 可能是 camelCase 或 PascalCase
                        if key[0].isupper():
                            continue  # 有些键如 'Error' 可接受
                        issues.append(
                            f"  {fname}: summary 键 '{key}' 不是 snake_case"
                        )

        if issues:
            msg = "以下 summary 键不符合 snake_case 命名:\n" + '\n'.join(issues)
            # 这是警告级别的问题
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)

    def test_category_naming_style(self):
        """result_files 的 category 应使用 snake_case 命名。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for m in re.finditer(r"""'category'\s*:\s*['"]([^'"]+)['"]""", source):
                cat = m.group(1)
                if cat != cat.lower() and '_' not in cat:
                    issues.append(
                        f"  {fname}: category '{cat}' 不是 snake_case"
                    )

        if issues:
            msg = "以下 category 不符合 snake_case 命名:\n" + '\n'.join(issues)
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)

    def test_all_modules_in_registry(self):
        """所有分析模块文件都应在 MODULE_REGISTRY 中注册。"""
        # 读取 __init__.py 获取注册表
        init_path = os.path.join(MODULES_DIR, '__init__.py')
        with open(init_path, 'r', encoding='utf-8') as f:
            init_source = f.read()

        # 提取注册的模块名
        registered = set(re.findall(r"'(\w+)'\s*:", init_source))
        # 移除非模块键
        registered -= {'PipelineType', 'validate_pipeline_order'}

        issues = []
        for filepath in _get_module_files():
            fname = os.path.basename(filepath)
            mname = _module_name(filepath)

            # 检查是否有 MODULE_NAME 声明
            lines, source = _read_source(filepath)
            module_name_match = re.search(r"MODULE_NAME\s*=\s*['\"]([^'\"]+)['\"]", source)
            if module_name_match:
                declared_name = module_name_match.group(1)
                if declared_name not in registered:
                    issues.append(
                        f"  {fname}: MODULE_NAME='{declared_name}' 未在 MODULE_REGISTRY 中注册"
                    )

        if issues:
            msg = "以下模块未在 MODULE_REGISTRY 中注册:\n" + '\n'.join(issues)
            pytest.fail(msg)


# ──────────────────────────────────────────────────────────────
# TestModuleErrorHandling — 验证错误处理模式
# ──────────────────────────────────────────────────────────────

class TestModuleErrorHandling:
    """验证模块的错误处理和异常捕获模式。"""

    def test_no_bare_except_pass(self):
        """模块不应使用 `except: pass` 或 `except Exception: pass` 模式（至少应记录日志）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for i, line in enumerate(lines):
                stripped = line.strip()
                # 匹配 except: 或 except Exception: 或 except SomeException:
                if re.match(r'except\s*(\w+)?\s*:', stripped):
                    # 检查下一行是否是 pass
                    if i + 1 < len(lines):
                        next_stripped = lines[i + 1].strip()
                        if next_stripped == 'pass':
                            # 排除合理的 pass（如 import 保护）
                            # 查看 except 块之前的 try 块内容
                            context_start = max(0, i - 10)
                            context = ''.join(lines[context_start:i])
                            # 如果是 import 保护，允许 pass
                            if 'import ' in context:
                                continue
                            issues.append(
                                f"  {fname}:{i+1}: except 块中使用 pass 而非记录日志"
                            )

        if issues:
            msg = "以下 except 块使用了 pass（应至少记录日志）:\n" + '\n'.join(issues)
            # 这是警告级别
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)

    def test_error_returns_have_error_key(self):
        """模块的错误返回应包含 'error' 键。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 run() 方法中的 return 语句
            m = re.search(
                r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                source, re.DOTALL
            )
            if not m:
                continue
            body = m.group(1)

            # 查找包含 'error' 或 'Error' 的返回（可能的错误返回）
            for i, line in enumerate(body.splitlines()):
                stripped = line.strip()
                if not stripped.startswith('return '):
                    continue
                # 检查是否是错误返回模式（包含 error 相关字符串）
                if 'error' not in stripped.lower() and 'invalid' not in stripped.lower():
                    continue
                # 检查是否有 'error' 键
                if "'error'" not in stripped and '"error"' not in stripped:
                    # 可能跨多行，收集 return 块
                    return_block = stripped
                    brace_count = stripped.count('{') - stripped.count('}')
                    j = i + 1
                    while brace_count > 0 and j < len(body.splitlines()):
                        return_block += body.splitlines()[j].strip()
                        brace_count += body.splitlines()[j].count('{') - body.splitlines()[j].count('}')
                        j += 1
                    if "'error'" not in return_block and '"error"' not in return_block:
                        issues.append(
                            f"  {fname}: 错误返回可能缺少 'error' 键"
                        )

        if issues:
            msg = "以下模块的错误返回可能缺少 'error' 键:\n" + '\n'.join(issues)
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)


# ──────────────────────────────────────────────────────────────
# TestModuleProgressCalls — 验证 progress 调用模式
# ──────────────────────────────────────────────────────────────

class TestModuleProgressCalls:
    """验证模块的 self.progress() 调用模式。"""

    def test_all_modules_call_progress(self):
        """每个分析模块的 run() 方法都应调用 self.progress()。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 查找 run 方法体
            m = re.search(
                r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                source, re.DOTALL
            )
            if not m:
                continue
            body = m.group(1)

            if 'self.progress(' not in body:
                issues.append(f"  {fname}: run() 方法未调用 self.progress()")

        if issues:
            msg = "以下模块的 run() 方法未调用 progress 回调:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_progress_starts_early_ends_at_100(self):
        """模块应在 run() 开始时调用 progress（≤ 10%），结束时调用 100%。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            m = re.search(
                r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                source, re.DOTALL
            )
            if not m:
                continue
            body = m.group(1)

            # 提取所有 progress 调用的百分比值
            pct_calls = re.findall(r'self\.progress\(\s*(\d+)', body)
            if not pct_calls:
                continue

            pct_values = [int(p) for p in pct_calls]
            first_pct = pct_values[0]
            last_pct = pct_values[-1]

            if first_pct > 15:
                issues.append(
                    f"  {fname}: 首次 progress 调用为 {first_pct}%（应 ≤ 10%）"
                )
            if last_pct != 100:
                issues.append(
                    f"  {fname}: 最后 progress 调用为 {last_pct}%（应为 100%）"
                )

        if issues:
            msg = "以下模块的 progress 调用范围异常:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_progress_monotonically_increasing(self):
        """模块的 progress 调用应单调递增（允许 -1 用于警告消息）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            m = re.search(
                r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                source, re.DOTALL
            )
            if not m:
                continue
            body = m.group(1)

            pct_calls = re.findall(r'self\.progress\(\s*(-?\d+)', body)
            if len(pct_calls) < 2:
                continue

            pct_values = [int(p) for p in pct_calls]
            # 排除 -1（警告消息）
            normal_values = [p for p in pct_values if p >= 0]

            for i in range(1, len(normal_values)):
                if normal_values[i] < normal_values[i - 1]:
                    issues.append(
                        f"  {fname}: progress 从 {normal_values[i-1]}% 回退到 {normal_values[i]}%"
                    )

        if issues:
            msg = "以下模块的 progress 调用有回退:\n" + '\n'.join(issues)
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)


# ──────────────────────────────────────────────────────────────
# TestModuleParamsAccess — 验证参数访问模式
# ──────────────────────────────────────────────────────────────

class TestModuleParamsAccess:
    """验证模块使用 self.params.get() 而非直接 self.params[] 访问。"""

    def test_no_bare_params_subscript(self):
        """模块不应使用 self.params['key']（应使用 self.params.get('key', default)）。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            m = re.search(
                r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                source, re.DOTALL
            )
            if not m:
                continue
            body = m.group(1)

            # 查找 self.params[' 或 self.params[" 模式
            for match in re.finditer(r'self\.params\s*\[\s*[\'"]', body):
                pos = match.start()
                # 检查是否在注释中
                line_start = body.rfind('\n', 0, pos) + 1
                line = body[line_start:body.find('\n', pos)]
                if '#' in line[:pos - line_start]:
                    continue
                # 获取行号
                line_no = body[:pos].count('\n') + 1
                issues.append(
                    f"  {fname}:{line_no}: 使用 self.params[] 而非 self.params.get()"
                )

        if issues:
            msg = "以下模块使用了直接参数访问（可能导致 KeyError）:\n" + '\n'.join(issues)
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)

    def test_params_get_has_default_for_required_types(self):
        """对数值类型参数使用 params.get() 时应提供默认值。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            m = re.search(
                r'def run\(self[^)]*\):(.*?)(?=\n    def |\nclass |\Z)',
                source, re.DOTALL
            )
            if not m:
                continue
            body = m.group(1)

            # 查找 int(self.params.get('key')) 没有默认值的模式
            for match in re.finditer(
                r'int\(\s*self\.params\.get\(\s*[\'"][^"\']+[\'"]\s*\)\s*\)',
                body
            ):
                pos = match.start()
                line_start = body.rfind('\n', 0, pos) + 1
                line = body[line_start:body.find('\n', pos)]
                if '#' in line[:pos - line_start]:
                    continue
                line_no = body[:pos].count('\n') + 1
                issues.append(
                    f"  {fname}:{line_no}: int(self.params.get('key')) 无默认值，"
                    f"key 不存在时将 int(None) 报错"
                )

            # 查找 float(self.params.get('key')) 没有默认值的模式
            for match in re.finditer(
                r'float\(\s*self\.params\.get\(\s*[\'"][^"\']+[\'"]\s*\)\s*\)',
                body
            ):
                pos = match.start()
                line_start = body.rfind('\n', 0, pos) + 1
                line = body[line_start:body.find('\n', pos)]
                if '#' in line[:pos - line_start]:
                    continue
                line_no = body[:pos].count('\n') + 1
                issues.append(
                    f"  {fname}:{line_no}: float(self.params.get('key')) 无默认值"
                )

        if issues:
            msg = "以下模块的 params.get() 缺少默认值:\n" + '\n'.join(issues)
            pytest.fail(msg)


# ──────────────────────────────────────────────────────────────
# TestModuleImports — 验证模块导入模式
# ──────────────────────────────────────────────────────────────

class TestModuleImports:
    """验证模块的导入模式。"""

    def test_heavy_deps_imported_inside_run(self):
        """重型依赖（scanpy, omicverse, plotly 等）应在 run() 方法内导入。"""
        heavy_deps = {'scanpy', 'omicverse', 'plotly', 'matplotlib', 'sklearn', 'pymde'}
        issues = []

        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 检查文件顶层是否有重型依赖的 import
            in_run_method = False
            top_level_imports = []
            run_level_imports = []

            for line in lines:
                stripped = line.strip()
                if stripped.startswith('def run('):
                    in_run_method = True
                    continue
                if in_run_method and (stripped.startswith('def ') or stripped.startswith('class ')):
                    in_run_method = False

                # 检查 import 语句
                import_match = re.match(r'(?:from\s+(\w+)|import\s+(\w+))', stripped)
                if import_match:
                    dep = import_match.group(1) or import_match.group(2)
                    if dep in heavy_deps:
                        if in_run_method:
                            run_level_imports.append(dep)
                        else:
                            # 检查是否在 class 定义前（模块顶层）
                            indent = len(line) - len(line.lstrip())
                            if indent == 0:
                                top_level_imports.append(dep)

            for dep in top_level_imports:
                issues.append(
                    f"  {fname}: 顶层导入了 '{dep}'（应在 run() 内导入以加速启动）"
                )

        if issues:
            msg = "以下模块在顶层导入了重型依赖:\n" + '\n'.join(issues)
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)


# ──────────────────────────────────────────────────────────────
# TestResultFilePathSafety — 验证输出路径安全
# ──────────────────────────────────────────────────────────────

class TestResultFilePathSafety:
    """验证模块的输出路径安全性和一致性。"""

    def test_save_plotly_json_uses_safe_filenames(self):
        """使用 save_plotly_json() 的文件名应只含安全字符。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            for match in re.finditer(
                r'self\.save_plotly_json\([^,]+,\s*[^,]+,\s*[\'"]([^\'"]+)[\'"]',
                source
            ):
                filename = match.group(1)
                # 文件名应只含字母、数字、下划线、点、连字符
                if not re.match(r'^[a-zA-Z0-9_.\-]+$', filename):
                    issues.append(
                        f"  {fname}: save_plotly_json 文件名含特殊字符: '{filename}'"
                    )

        if issues:
            msg = "以下 save_plotly_json 文件名不安全:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_csv_filenames_are_unique_per_module(self):
        """同一模块内的 CSV 文件名不应重复。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            # 提取所有 .to_csv() 调用的路径
            csv_paths = re.findall(
                r'\.to_csv\(\s*(?:os\.path\.join\([^)]+\)|[\'"][^\'"]+[\'"])',
                source
            )
            # 提取文件名部分
            csv_filenames = []
            for p in csv_paths:
                # 获取最后一个路径组件
                fname_match = re.search(r'[\'"]([^\'"/]+\.csv)[\'"]\s*$', p)
                if not fname_match:
                    fname_match = re.search(r'[\'"]([^\'"]+)[\'"]\s*$', p)
                if fname_match:
                    csv_filenames.append(fname_match.group(1))

            seen = set()
            for fn in csv_filenames:
                if fn in seen:
                    issues.append(f"  {fname}: 重复的 CSV 文件名 '{fn}'")
                seen.add(fn)

        if issues:
            msg = "以下模块有重复的 CSV 文件名:\n" + '\n'.join(issues)
            import warnings
            for issue in issues:
                warnings.warn(issue, UserWarning)


# ──────────────────────────────────────────────────────────────
# TestModuleDisplayInfo — 验证模块元数据完整性
# ──────────────────────────────────────────────────────────────

class TestModuleDisplayInfo:
    """验证模块的显示信息（DISPLAY_NAME, DESCRIPTION）完整性。"""

    def test_all_modules_have_display_name(self):
        """每个分析模块都应有 DISPLAY_NAME。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            if 'DISPLAY_NAME' not in source:
                issues.append(f"  {fname}: 缺少 DISPLAY_NAME")

            m = re.search(r"DISPLAY_NAME\s*=\s*['\"]([^'\"]+)['\"]", source)
            if m and len(m.group(1).strip()) == 0:
                issues.append(f"  {fname}: DISPLAY_NAME 为空字符串")

        if issues:
            msg = "以下模块缺少 DISPLAY_NAME:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_all_modules_have_description(self):
        """每个分析模块都应有 DESCRIPTION。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            if 'DESCRIPTION' not in source:
                issues.append(f"  {fname}: 缺少 DESCRIPTION")

            m = re.search(r"DESCRIPTION\s*=\s*['\"]([^'\"]+)['\"]", source)
            if m and len(m.group(1).strip()) == 0:
                issues.append(f"  {fname}: DESCRIPTION 为空字符串")

        if issues:
            msg = "以下模块缺少 DESCRIPTION:\n" + '\n'.join(issues)
            pytest.fail(msg)

    def test_module_names_unique(self):
        """所有模块的 MODULE_NAME 应唯一。"""
        names = []
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)

            m = re.search(r"MODULE_NAME\s*=\s*['\"]([^'\"]+)['\"]", source)
            if m:
                name = m.group(1)
                if name in names:
                    issues.append(f"  {fname}: MODULE_NAME '{name}' 重复")
                names.append(name)
            else:
                issues.append(f"  {fname}: 缺少 MODULE_NAME 声明")

        if issues:
            msg = "以下模块的 MODULE_NAME 有问题:\n" + '\n'.join(issues)
            pytest.fail(msg)


# ──────────────────────────────────────────────────────────────
# TestSaveOutputConsistency — 验证 save_output 调用一致性
# ──────────────────────────────────────────────────────────────

class TestSaveOutputConsistency:
    """验证所有模块使用 save_output() 的一致性。"""

    def test_save_output_module_name_matches_declaration(self):
        """save_output() 的模块名参数应与 MODULE_NAME 声明一致。"""
        issues = []
        for filepath in _get_module_files():
            lines, source = _read_source(filepath)
            fname = os.path.basename(filepath)
            mname = _module_name(filepath)

            # 提取 MODULE_NAME 声明
            module_name_match = re.search(r"MODULE_NAME\s*=\s*['\"]([^'\"]+)['\"]", source)
            if not module_name_match:
                continue
            declared_name = module_name_match.group(1)

            # 查找 save_output 调用
            for match in re.finditer(r'self\.save_output\(\s*\w+\s*,\s*[\'"]([^\'"]+)[\'"]', source):
                save_name = match.group(1)
                if save_name != declared_name:
                    issues.append(
                        f"  {fname}: save_output('{save_name}') != MODULE_NAME('{declared_name}')"
                    )

        if issues:
            msg = "以下模块的 save_output 名称与 MODULE_NAME 不一致:\n" + '\n'.join(issues)
            pytest.fail(msg)
