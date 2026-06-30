# tests/test_schemas.py
"""Tests for modules/schemas.py — parse_form_params, module metadata."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from modules.schemas import (
    parse_form_params, SC_MODULE_LIST, BULK_MODULE_LIST, MODULE_LIST,
    MODULE_DISPLAY_MAP, STATUS_MAP, PARAM_SCHEMAS,
)


class TestModuleMetadata:
    """测试模块元数据定义。"""

    def test_sc_module_count(self):
        """单细胞模块有 12 个。"""
        assert len(SC_MODULE_LIST) == 12

    def test_bulk_module_count(self):
        """Bulk 模块有 8 个。"""
        assert len(BULK_MODULE_LIST) == 8

    def test_module_list_is_sc_plus_bulk(self):
        """MODULE_LIST = SC + Bulk。"""
        assert MODULE_LIST == SC_MODULE_LIST + BULK_MODULE_LIST

    def test_display_map_complete(self):
        """每个模块都有中文显示名。"""
        for m in MODULE_LIST:
            assert m['name'] in MODULE_DISPLAY_MAP
            assert MODULE_DISPLAY_MAP[m['name']] == m['display']

    def test_status_map_keys(self):
        """状态映射包含所有状态。"""
        assert set(STATUS_MAP.keys()) == {'pending', 'running', 'completed', 'failed', 'processing'}

    def test_param_schemas_cover_all_modules(self):
        """PARAM_SCHEMAS 覆盖所有模块。"""
        for m in MODULE_LIST:
            assert m['name'] in PARAM_SCHEMAS, f"PARAM_SCHEMAS 缺少 {m['name']}"

    def test_param_schema_structure(self):
        """每个参数定义都有 key, label, type, default。"""
        for module_name, params in PARAM_SCHEMAS.items():
            for p in params:
                assert 'key' in p, f"{module_name}: 参数缺少 key"
                assert 'label' in p, f"{module_name}: 参数缺少 label"
                assert 'type' in p, f"{module_name}: 参数 {p.get('key')} 缺少 type"
                assert 'default' in p, f"{module_name}: 参数 {p.get('key')} 缺少 default"

    def test_param_types_valid(self):
        """参数类型都是合法值。"""
        valid_types = {'number', 'select', 'checkbox', 'text', 'textarea', 'dynamic_select', 'multiselect'}
        for module_name, params in PARAM_SCHEMAS.items():
            for p in params:
                assert p['type'] in valid_types, f"{module_name}.{p['key']}: 无效类型 {p['type']}"


class TestParseFormParams:
    """测试 parse_form_params 函数。"""

    def test_number_type(self):
        """number 类型：字符串转 float。"""
        schema = [{'key': 'resolution', 'type': 'number', 'default': 1.0}]
        form = {'resolution': '0.5'}
        result = parse_form_params(schema, form)
        assert result['resolution'] == 0.5
        assert isinstance(result['resolution'], float)

    def test_number_default(self):
        """number 类型：空值使用默认值。"""
        schema = [{'key': 'resolution', 'type': 'number', 'default': 1.0}]
        form = {}
        result = parse_form_params(schema, form)
        assert result['resolution'] == 1.0

    def test_checkbox_checked(self):
        """checkbox 类型：'on' → True。"""
        schema = [{'key': 'batch_aware', 'type': 'checkbox', 'default': False}]
        form = {'batch_aware': 'on'}
        result = parse_form_params(schema, form)
        assert result['batch_aware'] is True

    def test_checkbox_unchecked(self):
        """checkbox 类型：不在 form 中 → False。"""
        schema = [{'key': 'batch_aware', 'type': 'checkbox', 'default': False}]
        form = {}
        result = parse_form_params(schema, form)
        assert result['batch_aware'] is False

    def test_select_type(self):
        """select 类型：保留字符串值。"""
        schema = [{'key': 'method', 'type': 'select', 'default': 'leiden'}]
        form = {'method': 'louvain'}
        result = parse_form_params(schema, form)
        assert result['method'] == 'louvain'

    def test_select_default(self):
        """select 类型：空值使用默认值。"""
        schema = [{'key': 'method', 'type': 'select', 'default': 'leiden'}]
        form = {}
        result = parse_form_params(schema, form)
        assert result['method'] == 'leiden'

    def test_multiple_params(self):
        """多个参数混合类型。"""
        schema = [
            {'key': 'method', 'type': 'select', 'default': 'leiden'},
            {'key': 'resolution', 'type': 'number', 'default': 1.0},
            {'key': 'batch_aware', 'type': 'checkbox', 'default': False},
        ]
        form = {'method': 'louvain', 'resolution': '0.8', 'batch_aware': 'on'}
        result = parse_form_params(schema, form)
        assert result['method'] == 'louvain'
        assert result['resolution'] == 0.8
        assert result['batch_aware'] is True

    def test_text_type(self):
        """text 类型：保留字符串。"""
        schema = [{'key': 'batch_key', 'type': 'text', 'default': 'batch'}]
        form = {'batch_key': 'sample_id'}
        result = parse_form_params(schema, form)
        assert result['batch_key'] == 'sample_id'

    def test_empty_schema(self):
        """空 schema → 空 dict。"""
        result = parse_form_params([], {'anything': 'value'})
        assert result == {}
