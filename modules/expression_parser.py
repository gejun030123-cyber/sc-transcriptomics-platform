# modules/expression_parser.py
"""
多组差异分析表达式解析器。

支持语法：
    atom    := comparison:direction | comparison:direction(thresholds)
              | ALL:dir | ANY:dir | SHARED:dir | ONLY[comp]:dir
    expr    := xor_expr (XOR xor_expr)*
    xor_expr := and_expr ((AND) and_expr)*
    and_expr := NOT and_expr | atom | (expr)

运算符别名：
    AND / ∩ / &    交集
    OR  / ∪ / |    并集
    NOT / - / \\    差集
    XOR / △        对称差集

⚠️ 运算符优先级（从高到低）：
    一元 NOT > OR > AND > XOR > 二元 NOT（差集）
注意这与常见布尔优先级（AND 高于 OR）相反：
    ``A:up OR B:up AND C:up`` 解析为 ``(A OR B) AND C``。
请显式加括号表达复杂逻辑，避免歧义。
"""
import re


class ParseError(Exception):
    pass


# --- AST 节点 ---

class Atom:
    __slots__ = ('comparison', 'direction', 'padj_override', 'logfc_override')

    def __init__(self, comparison, direction, padj_override=None, logfc_override=None):
        self.comparison = comparison
        self.direction = direction  # 'up', 'down', 'both'
        self.padj_override = padj_override
        self.logfc_override = logfc_override

    def to_dict(self):
        d = {'type': 'atom', 'comparison': self.comparison, 'direction': self.direction}
        if self.padj_override is not None:
            d['padj'] = self.padj_override
        if self.logfc_override is not None:
            d['logfc'] = self.logfc_override
        return d


class BinOp:
    __slots__ = ('op', 'left', 'right')

    def __init__(self, op, left, right):
        self.op = op  # 'AND', 'OR', 'XOR'
        self.left = left
        self.right = right

    def to_dict(self):
        return {'type': 'binop', 'op': self.op,
                'left': self.left.to_dict(), 'right': self.right.to_dict()}


class UnaryNot:
    __slots__ = ('operand',)

    def __init__(self, operand):
        self.operand = operand

    def to_dict(self):
        return {'type': 'not', 'operand': self.operand.to_dict()}


class Shorthand:
    __slots__ = ('kind', 'direction', 'target')

    def __init__(self, kind, direction, target=None):
        self.kind = kind  # 'ALL', 'ANY', 'SHARED', 'ONLY'
        self.direction = direction
        self.target = target  # only for ONLY[comp]

    def to_dict(self):
        d = {'type': 'shorthand', 'kind': self.kind, 'direction': self.direction}
        if self.target:
            d['target'] = self.target
        return d


# --- Tokenizer ---

OP_ALIASES = {
    'AND': 'AND', '∩': 'AND', '&': 'AND',
    'OR': 'OR', '∪': 'OR', '|': 'OR',
    'NOT': 'NOT', '-': 'NOT', '\\': 'NOT',
    'XOR': 'XOR', '△': 'XOR',
}

DIRECTION_SET = {'up', 'down', 'both'}


def _parse_thresholds(thresh_str):
    """Parse 'padj<<0.01|logFC>2' → (padj_val, logfc_val)"""
    padj_val = None
    logfc_val = None
    for part in thresh_str.split('|'):
        part = part.strip()
        if '<<' in part:
            padj_val = float(part.split('<<')[1].strip())
        elif '>' in part:
            logfc_val = float(part.split('>')[1].strip())
    return padj_val, logfc_val


def tokenize(expr):
    """Tokenize expression string into token list."""
    tokens = []
    i = 0
    s = expr.strip()
    while i < len(s):
        # Skip whitespace
        if s[i].isspace():
            i += 1
            continue

        # Parentheses
        if s[i] == '(':
            tokens.append(('LPAREN', '('))
            i += 1
            continue
        if s[i] == ')':
            tokens.append(('RPAREN', ')'))
            i += 1
            continue

        # Multi-char operators: AND, OR, NOT, XOR (check BEFORE atoms)
        m = re.match(r'(AND|OR|NOT|XOR)\b', s[i:], re.IGNORECASE)
        if m:
            op = m.group(1).upper()
            tokens.append(('OP', op))
            i += m.end()
            continue

        # Single-char operators: ∩ ∪ & | \ - △
        if s[i] in '∩∪&|\\-△':
            op = OP_ALIASES.get(s[i], None)
            if op:
                tokens.append(('OP', op))
                i += 1
                continue

        # ONLY[comp]:dir shorthand
        m = re.match(r'ONLY\[([^\]]+)\]\s*:\s*(up|down|both)', s[i:], re.IGNORECASE)
        if m:
            tokens.append(('SHORTHAND', ('ONLY', m.group(2).lower(), m.group(1).strip())))
            i += m.end()
            continue

        # ALL/ANY/SHARED:dir shorthand
        m = re.match(r'(ALL|ANY|SHARED)\s*:\s*(up|down|both)', s[i:], re.IGNORECASE)
        if m:
            tokens.append(('SHORTHAND', (m.group(1).upper(), m.group(2).lower(), None)))
            i += m.end()
            continue

        # Atom: comparison:direction with optional (thresholds)
        # Match comparison name (any non-whitespace, non-colon, non-paren chars)
        # followed by :direction and optional (thresholds)
        m = re.match(
            r'([^\s:()]+?)\s*:\s*(up|down|both)'
            r'(?:\s*\(([^)]+)\))?',
            s[i:], re.IGNORECASE
        )
        if m:
            comp = m.group(1).strip()
            direction = m.group(2).lower()
            thresh_str = m.group(3)
            padj_ov = logfc_ov = None
            if thresh_str:
                padj_ov, logfc_ov = _parse_thresholds(thresh_str)
            tokens.append(('ATOM', Atom(comp, direction, padj_ov, logfc_ov)))
            i += m.end()
            continue

        # Unknown character
        raise ParseError(f"无法解析字符 '{s[i]}'（位置 {i}）")

    return tokens


# --- Parser (recursive descent) ---

class _Parser:
    def __init__(self, tokens):
        self.tokens = tokens
        self.pos = 0

    def peek(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def consume(self):
        tok = self.tokens[self.pos]
        self.pos += 1
        return tok

    def parse(self):
        if not self.tokens:
            raise ParseError("表达式为空")
        node = self._xor_expr()
        if self.pos < len(self.tokens):
            tok = self.tokens[self.pos]
            raise ParseError(f"表达式末尾有多余内容: {tok}")
        return node

    def _xor_expr(self):
        left = self._and_expr()
        while self.peek() and self.peek()[0] == 'OP' and self.peek()[1] == 'XOR':
            self.consume()
            right = self._and_expr()
            left = BinOp('XOR', left, right)
        # Binary NOT (set-difference): lower precedence than AND
        if self.peek() and self.peek()[0] == 'OP' and self.peek()[1] == 'NOT':
            self.consume()
            right = self._xor_expr()
            return BinOp('NOT', left, right)
        return left

    def _and_expr(self):
        # AND 的优先级低于 OR（文档约定：OR > AND > XOR > 二元 NOT）。
        left = self._or_expr()
        while self.peek() and self.peek()[0] == 'OP' and self.peek()[1] == 'AND':
            self.consume()
            right = self._or_expr()
            left = BinOp('AND', left, right)
        return left

    def _unary_not(self):
        # 一元 NOT 只绑定单个 primary（原子/括号/简写），因此
        # NOT A:up OR B:up 解析为 (NOT A) OR B，与文档优先级一致；
        # 需要整体取反时请使用括号 NOT (A:up OR B:up)。
        if self.peek() and self.peek()[0] == 'OP' and self.peek()[1] == 'NOT':
            self.consume()
            operand = self._unary_not()
            return UnaryNot(operand)
        return self._primary()

    def _or_expr(self):
        left = self._unary_not()
        while self.peek() and self.peek()[0] == 'OP' and self.peek()[1] == 'OR':
            self.consume()
            right = self._unary_not()
            left = BinOp('OR', left, right)
        return left

    def _primary(self):
        tok = self.peek()
        if tok is None:
            raise ParseError("表达式意外结束")

        if tok[0] == 'LPAREN':
            self.consume()
            node = self._xor_expr()
            if not self.peek() or self.peek()[0] != 'RPAREN':
                raise ParseError("缺少右括号 ')'")
            self.consume()
            return node

        if tok[0] == 'SHORTHAND':
            self.consume()
            kind, direction, target = tok[1]
            return Shorthand(kind, direction, target)

        if tok[0] == 'ATOM':
            self.consume()
            return tok[1]

        raise ParseError(f"意外的 token: {tok}")


def parse(expr):
    """Parse expression string into AST."""
    tokens = tokenize(expr)
    return _Parser(tokens).parse()


# --- Evaluator ---

def _get_atom_genes(atom, gene_sets, padj_matrix, logfc_matrix,
                    default_pval, default_log2fc):
    """Get gene set for a single Atom node."""
    comp = atom.comparison
    if comp not in gene_sets:
        raise ParseError(f"比较 '{comp}' 不存在。可用比较: {', '.join(sorted(gene_sets.keys()))}")

    pval = atom.padj_override if atom.padj_override is not None else default_pval
    log2fc_thresh = atom.logfc_override if atom.logfc_override is not None else default_log2fc

    pvals = padj_matrix[comp]
    abs_fcs = abs(logfc_matrix[comp])
    sig_mask = (pvals < pval) & (abs_fcs >= log2fc_thresh)

    if atom.direction == 'up':
        return set(logfc_matrix.index[sig_mask & (logfc_matrix[comp] > 0)])
    elif atom.direction == 'down':
        return set(logfc_matrix.index[sig_mask & (logfc_matrix[comp] < 0)])
    else:  # both
        return set(logfc_matrix.index[sig_mask])


def _expand_shorthand(sh, all_comparisons, gene_sets, padj_matrix, logfc_matrix,
                      default_pval, default_log2fc):
    """Expand shorthand into equivalent set operations."""
    comps = sorted(all_comparisons)

    if sh.kind == 'ALL' or sh.kind == 'SHARED':
        # Intersection across all comparisons
        result = None
        for c in comps:
            atom = Atom(c, sh.direction)
            gs = _get_atom_genes(atom, gene_sets, padj_matrix, logfc_matrix,
                                 default_pval, default_log2fc)
            result = gs if result is None else result & gs
        return result or set()

    if sh.kind == 'ANY':
        result = set()
        for c in comps:
            atom = Atom(c, sh.direction)
            gs = _get_atom_genes(atom, gene_sets, padj_matrix, logfc_matrix,
                                 default_pval, default_log2fc)
            result |= gs
        return result

    if sh.kind == 'ONLY':
        target = sh.target
        if target not in gene_sets:
            raise ParseError(f"ONLY 目标比较 '{target}' 不存在")
        target_genes = _get_atom_genes(
            Atom(target, sh.direction), gene_sets, padj_matrix, logfc_matrix,
            default_pval, default_log2fc)
        others = set()
        for c in comps:
            if c != target:
                others |= _get_atom_genes(
                    Atom(c, sh.direction), gene_sets, padj_matrix, logfc_matrix,
                    default_pval, default_log2fc)
        return target_genes - others

    raise ParseError(f"未知简写类型: {sh.kind}")


def evaluate(ast, all_comparisons, gene_sets, padj_matrix, logfc_matrix,
             default_pval=0.05, default_log2fc=1.0):
    """Evaluate AST to produce a set of gene names."""
    if isinstance(ast, Atom):
        return _get_atom_genes(ast, gene_sets, padj_matrix, logfc_matrix,
                               default_pval, default_log2fc)

    if isinstance(ast, Shorthand):
        return _expand_shorthand(ast, all_comparisons, gene_sets,
                                 padj_matrix, logfc_matrix,
                                 default_pval, default_log2fc)

    if isinstance(ast, UnaryNot):
        all_genes = set(padj_matrix.index)
        inner = evaluate(ast.operand, all_comparisons, gene_sets,
                         padj_matrix, logfc_matrix, default_pval, default_log2fc)
        return all_genes - inner

    if isinstance(ast, BinOp):
        left = evaluate(ast.left, all_comparisons, gene_sets,
                        padj_matrix, logfc_matrix, default_pval, default_log2fc)
        right = evaluate(ast.right, all_comparisons, gene_sets,
                         padj_matrix, logfc_matrix, default_pval, default_log2fc)
        if ast.op == 'AND':
            return left & right
        if ast.op == 'OR':
            return left | right
        if ast.op == 'XOR':
            return left ^ right
        if ast.op == 'NOT':
            return left - right

    raise ParseError(f"未知 AST 节点类型: {type(ast)}")


def _collect_comparisons_from_ast(ast):
    """Walk AST and collect all comparison names referenced."""
    comps = set()
    if isinstance(ast, Atom):
        comps.add(ast.comparison)
    elif isinstance(ast, Shorthand):
        if ast.target:
            comps.add(ast.target)
    elif isinstance(ast, UnaryNot):
        comps |= _collect_comparisons_from_ast(ast.operand)
    elif isinstance(ast, BinOp):
        comps |= _collect_comparisons_from_ast(ast.left)
        comps |= _collect_comparisons_from_ast(ast.right)
    return comps


def validate(expr, available_comparisons):
    """
    Validate expression against available comparisons.

    Returns: (ast_or_None, error_msg_or_None, warnings_list)
    """
    if not expr or not expr.strip():
        return None, "表达式为空", []

    try:
        ast = parse(expr.strip())
    except ParseError as e:
        return None, str(e), []

    warnings = []
    used = _collect_comparisons_from_ast(ast)
    avail = set(available_comparisons)
    missing = used - avail
    if missing:
        return None, f"比较不存在: {', '.join(sorted(missing))}", warnings

    unused = avail - used
    if unused:
        warnings.append(f"未使用的比较: {', '.join(sorted(unused))}")

    return ast, None, warnings
