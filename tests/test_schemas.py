# tests/test_schemas.py
"""Tests for modules/schemas.py — parse_form_params, module metadata."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from modules.schemas import (
    parse_form_params, SC_MODULE_LIST, BULK_MODULE_LIST, MODULE_LIST,
    MODULE_DISPLAY_MAP, STATUS_MAP, PARAM_SCHEMAS, filter_active_params,
)


class TestModuleMetadata:
    """测试模块元数据定义。"""

    def test_sc_module_count(self):
        """主单细胞流程不再展示兼容用的批量 CSV 导出。"""
        assert len(SC_MODULE_LIST) == 18
        assert all(item['name'] != 'sc_csv_export' for item in SC_MODULE_LIST)

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
        valid_types = {'number', 'select', 'checkbox', 'text', 'textarea', 'dynamic_select', 'dynamic_multiselect', 'multiselect'}
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

    def test_dynamic_multiselect_joins_checked_values(self):
        from werkzeug.datastructures import MultiDict

        schema = [{'key': 'sample_display_groups', 'type': 'dynamic_multiselect', 'default': ''}]
        result = parse_form_params(schema, MultiDict([
            ('sample_display_groups', 'ctrl'), ('sample_display_groups', 'treat'),
        ]))
        assert result == {'sample_display_groups': 'ctrl,treat'}

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

    def test_show_if_skips_inactive_params(self):
        """show_if 不匹配的参数不应被表单解析保存。"""
        schema = [
            {'key': 'method', 'type': 'select', 'default': 'harmony'},
            {
                'key': 'harmony_theta',
                'type': 'number',
                'default': 2.0,
                'show_if': {'method': 'harmony'},
            },
            {
                'key': 'bbknn_neighbors_within_batch',
                'type': 'number',
                'default': 3,
                'show_if': {'method': 'bbknn'},
            },
        ]
        form = {
            'method': 'bbknn',
            'harmony_theta': '9',
            'bbknn_neighbors_within_batch': '4',
        }
        result = parse_form_params(schema, form)
        assert result == {
            'method': 'bbknn',
            'bbknn_neighbors_within_batch': 4.0,
        }

    def test_show_if_supports_checkbox_conditions(self):
        """show_if 应正确识别 checkbox 的开关状态。"""
        schema = [
            {'key': 'evaluate_correction', 'type': 'checkbox', 'default': True},
            {
                'key': 'evaluation_sample_size',
                'type': 'number',
                'default': 10000,
                'show_if': {'evaluate_correction': True},
            },
        ]
        assert parse_form_params(schema, {}) == {'evaluate_correction': False}
        assert parse_form_params(schema, {
            'evaluate_correction': 'on',
            'evaluation_sample_size': '500',
        }) == {
            'evaluate_correction': True,
            'evaluation_sample_size': 500.0,
        }

    def test_filter_active_params_reuses_schema_conditions(self):
        """非表单入口也应复用 show_if 过滤无关参数。"""
        schema = [
            {'key': 'method', 'type': 'select', 'default': 'harmony'},
            {'key': 'harmony_theta', 'type': 'number', 'default': 2.0, 'show_if': {'method': 'harmony'}},
            {'key': 'bbknn_neighbors_within_batch', 'type': 'number', 'default': 3, 'show_if': {'method': 'bbknn'}},
        ]
        result = filter_active_params(schema, {
            'method': 'bbknn',
            'harmony_theta': 9.0,
            'bbknn_neighbors_within_batch': 4.0,
        })
        assert result == {
            'method': 'bbknn',
            'bbknn_neighbors_within_batch': 4.0,
        }

    def test_bulk_heatmap_defaults_to_auto_scope_and_exposes_comparison_for_top_var(self):
        schema = PARAM_SCHEMAS['bulk_heatmap']
        by_key = {item['key']: item for item in schema}

        assert by_key['sample_display_mode']['default'] == 'auto'
        assert by_key['deg_comparison']['show_if'] == {
            'gene_import_source': ['top_var', 'deg'],
        }
        params = parse_form_params(schema, {
            'gene_import_source': 'top_var',
            'deg_comparison': '/project/results/bulk_deg_results_0.csv',
        })
        assert params['sample_display_mode'] == 'auto'
        assert params['deg_comparison'].endswith('bulk_deg_results_0.csv')

    def test_organoid_type_is_active_only_for_organoid_panel(self):
        """类器官类型只在选择 Organoid 面板时传入后端。"""
        schema = PARAM_SCHEMAS['annotation']
        result = filter_active_params(schema, {
            'marker_set': 'Universal',
            'organoid_type': 'kidney',
        })
        assert 'organoid_type' not in result

        result = filter_active_params(schema, {
            'marker_set': 'Organoid',
            'organoid_type': 'kidney',
        })
        assert result['organoid_type'] == 'kidney'

    def test_organoid_maturity_and_quality_controls_are_active(self):
        """成熟度时间元数据、负向 marker 与质量证据参数可由 Schema 传递。"""
        schema = PARAM_SCHEMAS['annotation']
        values = filter_active_params(schema, {
            'method': 'multi_evidence',
            'marker_set': 'Organoid',
            'maturity_time_key': 'culture_day',
            'annotation_version': 'v2',
            'annotation_comment': '人工确认肾单位边界',
            'state_score_threshold': 0.4,
            'negative_marker_weight': 0.7,
            'lineage_mixture_threshold': 0.25,
            'ambient_score_threshold': 0.4,
        })
        assert values['maturity_time_key'] == 'culture_day'
        assert values['annotation_version'] == 'v2'
        assert values['annotation_comment'] == '人工确认肾单位边界'
        assert values['state_score_threshold'] == 0.4
        assert values['negative_marker_weight'] == 0.7
        assert values['lineage_mixture_threshold'] == 0.25
        assert values['ambient_score_threshold'] == 0.4

        colorectal = filter_active_params(schema, {
            'method': 'multi_evidence',
            'marker_set': 'Colorectal',
            'maturity_time_key': 'culture_day',
        })
        assert colorectal['maturity_time_key'] == 'culture_day'

    def test_celltypist_reference_controls_are_active(self):
        """CellTypist is an explicit local reference option, not a default method."""
        schema = PARAM_SCHEMAS['annotation']
        defaults = {item['key']: item.get('default') for item in schema}
        assert defaults['use_celltypist_reference'] is False
        values = filter_active_params(schema, {
            'method': 'multi_evidence',
            'use_celltypist_reference': True,
            'celltypist_model': 'Cells_Intestinal_Tract.pkl',
            'celltypist_mode': 'prob match',
            'celltypist_p_threshold': 0.6,
            'celltypist_majority_voting': False,
        })
        assert values['use_celltypist_reference'] is True
        assert values['celltypist_model'] == 'Cells_Intestinal_Tract.pkl'
        assert values['celltypist_mode'] == 'prob match'
        assert values['celltypist_p_threshold'] == 0.6
