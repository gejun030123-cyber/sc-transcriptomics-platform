"""Regression coverage for the source-bound single-cell GO worker path."""


def test_sc_cell_go_worker_skips_irrelevant_h5ad_preflight(monkeypatch, tmp_path):
    import worker
    from modules import MODULE_REGISTRY

    completed = {}

    class FakeTask:
        id = 'task_go'

        def mark_running(self):
            return True

        def update_progress(self, *_args):
            return None

        def mark_completed(self, output_adata, result_json):
            completed['output_adata'] = output_adata
            completed['result_json'] = result_json

        def mark_failed(self, error_traceback, result_json=None):
            completed['failed'] = (error_traceback, result_json)

    class FakeConn:
        def execute(self, *_args):
            return self

        def fetchone(self):
            return {'id': 'project_go'}

        def close(self):
            return None

    class SourceBoundGO:
        def __init__(self, project_dir, params, progress_callback):
            self.project_dir = project_dir
            self.params = params
            self.progress_callback = progress_callback

        def load_adata(self, _input_path):
            raise AssertionError('sc_cell_go must not load the chained h5ad')

        def validate_input(self, _adata):
            raise AssertionError('sc_cell_go must not validate the chained h5ad')

        def run(self, _input_path):
            assert self.params['_analysis_id'] == 'task_go'
            return {'output_adata': _input_path, 'result_files': [], 'summary': {}}

    monkeypatch.setattr(worker.AnalysisTask, 'get_by_id', lambda _task_id: FakeTask())
    monkeypatch.setattr(worker, 'get_conn', lambda: FakeConn())
    monkeypatch.setattr(worker, 'register_task_outputs', lambda *_args, **_kwargs: [])
    monkeypatch.setattr(worker, '_refresh_project_status', lambda _project_id: None)
    monkeypatch.setitem(MODULE_REGISTRY, 'sc_cell_go', SourceBoundGO)
    worker._run_task(
        'task_go', 'project_go', 'sc_cell_go', {}, str(tmp_path), str(tmp_path / 'chained.h5ad'),
    )
    assert 'failed' not in completed
    assert completed['output_adata'].endswith('chained.h5ad')
