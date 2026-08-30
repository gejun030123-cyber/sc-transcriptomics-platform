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
