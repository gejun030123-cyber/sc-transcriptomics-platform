import json

import pytest


def test_llm_annotation_sends_cluster_only_summary(monkeypatch):
    from modules import llm_annotation

    monkeypatch.setattr(llm_annotation, 'get_effective_ai_config', lambda: {
        'api_url': 'https://llm.example/v1',
        'api_key': 'test-key',
        'model': 'test-model',
        'provider': 'openai',
    })
    captured = {}

    class FakeResponse:
        status_code = 200

        def json(self):
            return {
                'choices': [{
                    'message': {'content': json.dumps({
                        'annotations': [
                            {
                                'cluster': '0', 'cell_type': 'T cells',
                                'confidence': 'high', 'rationale': 'CD3D and TRAC.',
                                'review_note': '',
                            },
                            {
                                'cluster': '1', 'cell_type': 'B cells',
                                'confidence': 'medium', 'rationale': 'MS4A1 and CD79A.',
                                'review_note': 'Confirm with CD74.',
                            },
                        ],
                        'global_note': 'Research-only review required.',
                    })}
                }]
            }

    def fake_post(url, **kwargs):
        captured['url'] = url
        captured['body'] = kwargs['json']
        return FakeResponse()

    monkeypatch.setattr('requests.post', fake_post)
    result = llm_annotation.run_llm_cluster_annotation(
        [
            {
                'cluster': '0', 'n_cells': 12,
                'top_markers': ['CD3D', 'TRAC'],
                'marker_panel_candidate': 'T cells',
                'marker_panel_support': 'CD3D, TRAC',
                # The adapter only reads explicit allow-listed fields.
                'input_path': '/sensitive/project/input.h5ad',
                'donor_id': 'person-001',
            },
            {
                'cluster': '1', 'n_cells': 10,
                'top_markers': ['MS4A1', 'CD79A'],
                'marker_panel_candidate': 'B cells',
            },
        ],
        species='human', tissue_context='human PBMC', max_clusters=5,
    )

    request_text = captured['body']['messages'][1]['content']
    assert captured['url'] == 'https://llm.example/v1/chat/completions'
    assert '/sensitive/project/input.h5ad' not in request_text
    assert 'person-001' not in request_text
    assert result['annotations'][0]['cell_type'] == 'T cells'
    assert result['provider'] == 'openai'
    assert result['model'] == 'test-model'
    assert len(result['request_hash']) == 64


def test_llm_annotation_rejects_server_paths_in_tissue_context():
    from modules.llm_annotation import LLMAnnotationError, build_llm_annotation_payload

    with pytest.raises(LLMAnnotationError, match='服务器路径'):
        build_llm_annotation_payload(
            [{'cluster': '0', 'n_cells': 2, 'top_markers': ['CD3D']}],
            tissue_context='/home/private/project',
        )


def test_llm_annotation_chunks_many_clusters(monkeypatch):
    """超过单批上限时自动分批注释，所有 cluster 都获得标签。"""
    from modules import llm_annotation

    monkeypatch.setattr(llm_annotation, 'get_effective_ai_config', lambda: {
        'api_url': 'https://llm.example/v1',
        'api_key': 'test-key',
        'model': 'test-model',
        'provider': 'openai',
    })
    calls = []

    class FakeResponse:
        status_code = 200

        def json(self):
            return {
                'choices': [{
                    'message': {'content': json.dumps({
                        'annotations': [
                            {
                                'cluster': item['cluster'],
                                'cell_type': f'Type {item["cluster"]}',
                                'confidence': 'medium',
                                'rationale': 'marker evidence',
                                'review_note': '',
                            }
                            for item in self.batch
                        ],
                        'global_note': 'review all',
                    })},
                }],
            }

    def fake_post(url, **kwargs):
        body = kwargs['json']
        content = body['messages'][1]['content']
        payload = json.loads(content[content.index('{'):])
        calls.append([item['cluster'] for item in payload['clusters']])
        response = FakeResponse()
        response.batch = payload['clusters']
        return response

    monkeypatch.setattr('requests.post', fake_post)
    result = llm_annotation.run_llm_cluster_annotation(
        [
            {'cluster': str(i), 'n_cells': 5, 'top_markers': ['CD3D']}
            for i in range(5)
        ],
        max_clusters=2,
    )

    assert calls == [['0', '1'], ['2', '3'], ['4']]
    assert result['n_batches'] == 3
    assert result['max_clusters_per_batch'] == 2
    assert [item['cluster'] for item in result['annotations']] == ['0', '1', '2', '3', '4']
    assert all(item['cell_type'] == f'Type {item["cluster"]}'
               for item in result['annotations'])
    assert isinstance(result['request_hash'], list)
    assert len(result['request_hash']) == 3
    assert all(len(hash_value) == 64 for hash_value in result['request_hash'])
    assert len(result['cluster_input']) == 5


def _fake_settings(provider='openai'):
    return {
        'api_url': 'https://llm.example/v1',
        'api_key': 'test-key',
        'model': 'test-model',
        'provider': provider,
    }


def test_response_token_budget_scales_with_batch_size():
    """30 簇的批次不能再按 4096 token 请求（必然截断）。"""
    from modules.llm_annotation import (
        MAX_RESPONSE_TOKENS,
        MIN_RESPONSE_TOKENS,
        _response_token_budget,
    )

    assert _response_token_budget(1) == MIN_RESPONSE_TOKENS
    assert _response_token_budget(30) >= 30 * 240
    assert _response_token_budget(40) <= MAX_RESPONSE_TOKENS
    assert _response_token_budget(30) > 4096


def test_request_uses_dynamic_max_tokens(monkeypatch):
    from modules import llm_annotation

    monkeypatch.setattr(llm_annotation, 'get_effective_ai_config',
                        lambda: _fake_settings())
    captured = {}

    class FakeResponse:
        status_code = 200

        def json(self):
            return {'choices': [{'message': {'content': json.dumps({
                'annotations': [
                    {'cluster': item['cluster'], 'cell_type': 'T cells',
                     'confidence': 'high', 'rationale': 'CD3D', 'review_note': ''}
                    for item in self.batch
                ],
                'global_note': '',
            })}}]}

    def fake_post(url, **kwargs):
        body = kwargs['json']
        content = body['messages'][1]['content']
        payload = json.loads(content[content.index('{'):])
        captured['max_tokens'] = body['max_tokens']
        captured['temperature'] = body['temperature']
        captured['n_clusters'] = len(payload['clusters'])
        response = FakeResponse()
        response.batch = payload['clusters']
        return response

    monkeypatch.setattr('requests.post', fake_post)
    llm_annotation.run_llm_cluster_annotation(
        [{'cluster': str(i), 'n_cells': 5, 'top_markers': ['CD3D']} for i in range(4)],
        max_clusters=4,
    )

    assert captured['n_clusters'] == 4
    assert captured['max_tokens'] >= 4 * 240
    assert captured['temperature'] == 0


def test_anthropic_request_pins_temperature(monkeypatch):
    from modules import llm_annotation

    monkeypatch.setattr(llm_annotation, 'get_effective_ai_config',
                        lambda: _fake_settings(provider='anthropic'))
    captured = {}

    class FakeResponse:
        status_code = 200

        def json(self):
            return {'content': [{'type': 'text', 'text': json.dumps({
                'annotations': [
                    {'cluster': item['cluster'], 'cell_type': 'T cells',
                     'confidence': 'high', 'rationale': 'CD3D', 'review_note': ''}
                    for item in self.batch
                ],
                'global_note': '',
            })}]}

    def fake_post(url, **kwargs):
        body = kwargs['json']
        payload = json.loads(body['messages'][0]['content'].split('\n', 1)[1])
        captured['temperature'] = body.get('temperature')
        captured['max_tokens'] = body['max_tokens']
        response = FakeResponse()
        response.batch = payload['clusters']
        return response

    monkeypatch.setattr('requests.post', fake_post)
    llm_annotation.run_llm_cluster_annotation(
        [{'cluster': '0', 'n_cells': 5, 'top_markers': ['CD3D']}],
        max_clusters=4,
    )

    assert captured['temperature'] == 0


def test_truncated_batch_is_retried_with_smaller_batches(monkeypatch):
    """首轮被截断（少返回 cluster）时自动减半重试，而不是整任务失败。"""
    from modules import llm_annotation

    monkeypatch.setattr(llm_annotation, 'get_effective_ai_config',
                        lambda: _fake_settings())
    requests_seen = []

    class FakeResponse:
        status_code = 200

        def json(self):
            # 模拟截断：只返回该批次的一部分 cluster
            returned = self.batch[:1] if len(self.batch) > 1 else self.batch
            return {'choices': [{'message': {'content': json.dumps({
                'annotations': [
                    {'cluster': item['cluster'], 'cell_type': f"Type {item['cluster']}",
                     'confidence': 'low', 'rationale': 'partial', 'review_note': ''}
                    for item in returned
                ],
                'global_note': '',
            })}}]}

    def fake_post(url, **kwargs):
        content = kwargs['json']['messages'][1]['content']
        payload = json.loads(content[content.index('{'):])
        requests_seen.append([item['cluster'] for item in payload['clusters']])
        response = FakeResponse()
        response.batch = payload['clusters']
        return response

    monkeypatch.setattr('requests.post', fake_post)
    result = llm_annotation.run_llm_cluster_annotation(
        [{'cluster': str(i), 'n_cells': 5, 'top_markers': ['CD3D']} for i in range(4)],
        max_clusters=4,
    )

    assert [item['cluster'] for item in result['annotations']] == ['0', '1', '2', '3']
    assert result['n_batches'] == 1
    assert result['n_requests_sent'] > 1
    assert result['batch_splits']
    # 每次拆分都记录在案，且最终发送的是被拆小后的批次
    assert requests_seen[0] == ['0', '1', '2', '3']
    assert {cluster for batch in requests_seen[1:] for cluster in batch} == \
        {'0', '1', '2', '3'}
    assert all(len(batch) <= 2 for batch in requests_seen[1:])
