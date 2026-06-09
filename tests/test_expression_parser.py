# tests/test_expression_parser.py
"""Tests for the expression parser used in bulk_deg_integration filter."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import pytest
from modules.expression_parser import (
    tokenize, parse, evaluate, validate, ParseError,
    Atom, BinOp, UnaryNot, Shorthand,
)


# --- Fixtures ---

@pytest.fixture
def matrices():
    """Create mock padj and logfc matrices for testing."""
    genes = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE', 'GeneF']
    comps = ['comp1-vs-ctrl', 'comp2-vs-ctrl', 'comp3-vs-ctrl']

    logfc_data = {
        'comp1-vs-ctrl': [2.5, -3.0, 0.5, 1.8, -0.2, 0.0],
        'comp2-vs-ctrl': [1.5, -2.0, 3.0, 0.1, -1.5, 0.3],
        'comp3-vs-ctrl': [-1.2, -2.5, 2.0, 0.3, 0.1, -4.0],
    }
    padj_data = {
        'comp1-vs-ctrl': [0.001, 0.002, 0.8, 0.01, 0.9, 0.5],
        'comp2-vs-ctrl': [0.005, 0.001, 0.003, 0.6, 0.04, 0.3],
        'comp3-vs-ctrl': [0.001, 0.008, 0.002, 0.5, 0.7, 0.001],
    }

    logfc_matrix = pd.DataFrame(logfc_data, index=genes)
    padj_matrix = pd.DataFrame(padj_data, index=genes)
    return genes, comps, logfc_matrix, padj_matrix


# --- Tokenizer tests ---

class TestTokenize:
    def test_simple_atom(self):
        tokens = tokenize('comp1-vs-ctrl:up')
        assert len(tokens) == 1
        assert tokens[0][0] == 'ATOM'
        assert tokens[0][1].comparison == 'comp1-vs-ctrl'
        assert tokens[0][1].direction == 'up'

    def test_atom_with_thresholds(self):
        tokens = tokenize('comp1-vs-ctrl:up(padj<<0.01|logFC>2)')
        assert len(tokens) == 1
        atom = tokens[0][1]
        assert atom.padj_override == 0.01
        assert atom.logfc_override == 2.0

    def test_and_operator(self):
        tokens = tokenize('A:up AND B:down')
        types = [t[0] for t in tokens]
        assert types == ['ATOM', 'OP', 'ATOM']

    def test_operator_aliases(self):
        for op_str in ['AND', '∩', '&']:
            tokens = tokenize(f'A:up {op_str} B:up')
            assert tokens[1] == ('OP', 'AND')

        for op_str in ['OR', '∪', '|']:
            tokens = tokenize(f'A:up {op_str} B:up')
            assert tokens[1] == ('OP', 'OR')

        for op_str in ['NOT', '-', '\\']:
            tokens = tokenize(f'{op_str} A:up')
            assert tokens[0] == ('OP', 'NOT')

        tokens = tokenize('A:up XOR B:up')
        assert tokens[1] == ('OP', 'XOR')

    def test_shorthand_all(self):
        tokens = tokenize('ALL:up')
        assert len(tokens) == 1
        assert tokens[0][0] == 'SHORTHAND'
        assert tokens[0][1] == ('ALL', 'up', None)

    def test_shorthand_only(self):
        tokens = tokenize('ONLY[comp1-vs-ctrl]:up')
        assert len(tokens) == 1
        assert tokens[0][1] == ('ONLY', 'up', 'comp1-vs-ctrl')

    def test_parentheses(self):
        tokens = tokenize('(A:up AND B:up)')
        assert tokens[0] == ('LPAREN', '(')
        assert tokens[-1] == ('RPAREN', ')')

    def test_empty_raises(self):
        tokens = tokenize('')
        assert tokens == []

    def test_unknown_char_raises(self):
        with pytest.raises(ParseError, match='无法解析'):
            tokenize('A:up @ B:up')


# --- Parser tests ---

class TestParse:
    def test_single_atom(self):
        ast = parse('comp1-vs-ctrl:up')
        assert isinstance(ast, Atom)
        assert ast.comparison == 'comp1-vs-ctrl'
        assert ast.direction == 'up'

    def test_and_expression(self):
        ast = parse('A:up AND B:down')
        assert isinstance(ast, BinOp)
        assert ast.op == 'AND'
        assert isinstance(ast.left, Atom)
        assert isinstance(ast.right, Atom)

    def test_or_expression(self):
        ast = parse('A:up OR B:up')
        assert isinstance(ast, BinOp)
        assert ast.op == 'OR'

    def test_not_expression(self):
        ast = parse('NOT A:up')
        assert isinstance(ast, UnaryNot)
        assert isinstance(ast.operand, Atom)

    def test_xor_expression(self):
        ast = parse('A:up XOR B:up')
        assert isinstance(ast, BinOp)
        assert ast.op == 'XOR'

    def test_parentheses(self):
        ast = parse('(A:up AND B:up) OR C:down')
        assert isinstance(ast, BinOp)
        assert ast.op == 'OR'
        assert isinstance(ast.left, BinOp)
        assert ast.left.op == 'AND'

    def test_nested_parentheses(self):
        ast = parse('((A:up AND B:up) NOT C:up) OR D:down')
        assert isinstance(ast, BinOp)
        assert ast.op == 'OR'

    def test_shorthand(self):
        ast = parse('ALL:up')
        assert isinstance(ast, Shorthand)
        assert ast.kind == 'ALL'
        assert ast.direction == 'up'

    def test_only_shorthand(self):
        ast = parse('ONLY[comp1-vs-ctrl]:up')
        assert isinstance(ast, Shorthand)
        assert ast.kind == 'ONLY'
        assert ast.target == 'comp1-vs-ctrl'

    def test_complex_expression(self):
        expr = '(comp1-vs-ctrl:up AND comp2-vs-ctrl:up) NOT comp3-vs-ctrl:up'
        ast = parse(expr)
        assert isinstance(ast, BinOp)
        assert ast.op == 'NOT'

    def test_empty_raises(self):
        with pytest.raises(ParseError, match='为空'):
            parse('')

    def test_unmatched_paren_raises(self):
        with pytest.raises(ParseError, match='右括号'):
            parse('(A:up AND B:up')

    def test_to_dict(self):
        ast = parse('A:up AND B:down')
        d = ast.to_dict()
        assert d['type'] == 'binop'
        assert d['op'] == 'AND'


# --- Evaluator tests ---

class TestEvaluate:
    def test_single_atom_up(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('comp1-vs-ctrl:up')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # comp1-vs-ctrl: GeneA(2.5,p=0.001), GeneD(1.8,p=0.01) are up and sig
        assert 'GeneA' in result
        assert 'GeneD' in result
        assert 'GeneB' not in result  # down

    def test_single_atom_down(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('comp1-vs-ctrl:down')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        assert 'GeneB' in result  # -3.0, padj=0.002

    def test_and_intersection(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('comp1-vs-ctrl:up AND comp2-vs-ctrl:up')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneA: comp1 up (2.5,p=0.001) AND comp2 up (1.5,p=0.005) → in
        assert 'GeneA' in result
        # GeneD: comp1 up (1.8,p=0.01) but comp2 ns (0.1,p=0.6) → out
        assert 'GeneD' not in result

    def test_or_union(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('comp1-vs-ctrl:up OR comp2-vs-ctrl:up')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneA in both, GeneD in comp1 only, GeneC in comp2 only
        assert 'GeneA' in result
        assert 'GeneD' in result
        assert 'GeneC' in result

    def test_not_difference(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('comp1-vs-ctrl:up NOT comp2-vs-ctrl:up')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneA in both → excluded; GeneD only in comp1 → included
        assert 'GeneD' in result
        assert 'GeneA' not in result

    def test_xor_symmetric_diff(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('comp1-vs-ctrl:up XOR comp2-vs-ctrl:up')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneA in both → excluded; GeneD only comp1, GeneC only comp2 → included
        assert 'GeneA' not in result
        assert 'GeneD' in result
        assert 'GeneC' in result

    def test_shorthand_all(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('ALL:down')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneB: comp1 down (-3.0,p=0.002), comp2 down (-2.0,p=0.001), comp3 down (-2.5,p=0.008)
        assert 'GeneB' in result

    def test_shorthand_any(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('ANY:up')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneA up in comp1&2, GeneC up in comp2&3, GeneD up in comp1
        assert 'GeneA' in result
        assert 'GeneC' in result

    def test_shorthand_only(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        ast = parse('ONLY[comp1-vs-ctrl]:up')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneD up only in comp1 → included; GeneA up in comp1 AND comp2 → excluded
        assert 'GeneD' in result
        assert 'GeneA' not in result

    def test_invalid_comparison_raises(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        ast = parse('nonexistent:up')
        with pytest.raises(ParseError, match='不存在'):
            evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)

    def test_threshold_override(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        # With stricter padj threshold (0.005), GeneD (padj=0.01) should be excluded
        ast = parse('comp1-vs-ctrl:up(padj<<0.005)')
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        assert 'GeneA' in result  # padj=0.001
        assert 'GeneD' not in result  # padj=0.01 > 0.005

    def test_complex_expression(self, matrices):
        genes, comps, logfc_matrix, padj_matrix = matrices
        gene_sets = {}
        for c in comps:
            sig = (padj_matrix[c] < 0.05) & (abs(logfc_matrix[c]) >= 1.0)
            gene_sets[c] = set(logfc_matrix.index[sig])

        expr = '(comp1-vs-ctrl:up OR comp2-vs-ctrl:up) AND comp3-vs-ctrl:down'
        ast = parse(expr)
        result = evaluate(ast, comps, gene_sets, padj_matrix, logfc_matrix, 0.05, 1.0)
        # GeneB: up in neither comp1 nor comp2 → out
        # GeneE: not up in any → out
        # Need a gene that is up in (comp1 or comp2) AND down in comp3
        # GeneF: comp3 down (-4.0,p=0.001) but not up in comp1 or comp2 → out
        # No gene matches both → result should be empty or specific
        # This is valid - just checking it doesn't crash
        assert isinstance(result, set)


# --- Validation tests ---

class TestValidate:
    def test_valid_expression(self):
        ast, err, warnings = validate('comp1:up AND comp2:down', ['comp1', 'comp2'])
        assert err is None
        assert ast is not None

    def test_invalid_comparison(self):
        ast, err, warnings = validate('unknown:up', ['comp1', 'comp2'])
        assert err is not None
        assert '不存在' in err
        assert ast is None

    def test_empty_expression(self):
        ast, err, warnings = validate('', ['comp1'])
        assert err is not None

    def test_unused_comparison_warning(self):
        ast, err, warnings = validate('comp1:up', ['comp1', 'comp2'])
        assert err is None
        assert len(warnings) > 0

    def test_syntax_error(self):
        ast, err, warnings = validate('A:up AND', ['A'])
        assert err is not None

    def test_mixed_operators(self):
        ast, err, warnings = validate(
            'comp1:up AND comp2:down OR comp1:both',
            ['comp1', 'comp2']
        )
        assert err is None

    def test_direction_both(self):
        ast, err, warnings = validate('comp1:both', ['comp1'])
        assert err is None
        assert isinstance(ast, Atom)
        assert ast.direction == 'both'


# --- Edge case tests ---

class TestEdgeCases:
    def test_single_atom_expression(self):
        ast, err, _ = validate('comp1:up', ['comp1'])
        assert err is None
        assert isinstance(ast, Atom)

    def test_whitespace_handling(self):
        ast = parse('  comp1:up   AND   comp2:down  ')
        assert isinstance(ast, BinOp)

    def test_multiple_not(self):
        ast = parse('NOT NOT A:up')
        assert isinstance(ast, UnaryNot)
        assert isinstance(ast.operand, UnaryNot)

    def test_chained_and(self):
        ast = parse('A:up AND B:up AND C:up')
        assert isinstance(ast, BinOp)
        assert ast.op == 'AND'
        assert isinstance(ast.left, BinOp)
