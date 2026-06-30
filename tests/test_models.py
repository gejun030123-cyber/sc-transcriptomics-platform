# tests/test_models.py
"""Tests for models.py — gen_id, model constructors, to_dict."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest


class TestGenId:
    """测试 gen_id 函数。"""

    def test_returns_string(self):
        """返回字符串。"""
        from models import gen_id
        result = gen_id()
        assert isinstance(result, str)

    def test_length_12(self):
        """长度为 12。"""
        from models import gen_id
        assert len(gen_id()) == 12

    def test_unique(self):
        """多次调用返回不同值。"""
        from models import gen_id
        ids = {gen_id() for _ in range(100)}
        assert len(ids) == 100  # 全部唯一

    def test_alphanumeric(self):
        """只包含字母、数字、连字符。"""
        import re
        from models import gen_id
        for _ in range(50):
            assert re.match(r'^[a-f0-9-]+$', gen_id())


class TestProjectModel:
    """测试 Project 模型。"""

    def test_default_values(self):
        """默认值正确。"""
        from models import Project
        p = Project()
        assert p.name == ''
        assert p.description == ''
        assert p.status == 'empty'
        assert p.metadata_json == '{}'

    def test_custom_values(self):
        """自定义值。"""
        from models import Project
        p = Project(id='test-123', name='My Project', description='desc', status='data_ready')
        assert p.id == 'test-123'
        assert p.name == 'My Project'
        assert p.status == 'data_ready'

    def test_to_dict(self):
        """to_dict 包含所有字段。"""
        from models import Project
        p = Project(id='test-123', name='Test')
        d = p.to_dict()
        assert d['id'] == 'test-123'
        assert d['name'] == 'Test'
        assert 'created_at' in d
        assert 'updated_at' in d
        assert 'status' in d

    def test_auto_id(self):
        """未指定 id 时自动生成。"""
        from models import Project
        p = Project()
        assert p.id is not None
        assert len(p.id) == 12


class TestAnalysisTaskModel:
    """测试 AnalysisTask 模型。"""

    def test_default_values(self):
        """默认值正确。"""
        from models import AnalysisTask
        t = AnalysisTask(project_id='proj-1', module_name='qc')
        assert t.project_id == 'proj-1'
        assert t.module_name == 'qc'
        assert t.status == 'pending'
        assert t.progress == 0

    def test_to_dict(self):
        """to_dict 包含所有字段。"""
        from models import AnalysisTask
        t = AnalysisTask(project_id='proj-1', module_name='qc')
        d = t.to_dict()
        assert d['project_id'] == 'proj-1'
        assert d['module_name'] == 'qc'
        assert d['status'] == 'pending'
        assert 'id' in d


class TestResultFileModel:
    """测试 ResultFile 模型。"""

    def test_default_values(self):
        """默认值正确。"""
        from models import ResultFile
        rf = ResultFile(task_id='task-1', project_id='proj-1', file_type='plotly_json', file_path='/tmp/test.json')
        assert rf.task_id == 'task-1'
        assert rf.file_type == 'plotly_json'

    def test_to_dict(self):
        """to_dict 包含所有字段。"""
        from models import ResultFile
        rf = ResultFile(task_id='task-1', project_id='proj-1', file_type='csv', file_path='/tmp/test.csv')
        d = rf.to_dict()
        assert d['file_type'] == 'csv'
        assert d['file_path'] == '/tmp/test.csv'
