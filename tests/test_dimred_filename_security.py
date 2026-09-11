def test_dimred_filename_token_cannot_traverse_plots_directory():
    from modules.dimred import _safe_filename_token

    assert _safe_filename_token('x/../../escape') == 'x_.._.._escape'
    assert _safe_filename_token('') == 'column'
    assert '/' not in _safe_filename_token('batch/one')
    assert '\\' not in _safe_filename_token('batch\\one')
