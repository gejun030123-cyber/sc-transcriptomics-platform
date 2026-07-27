"""Regression checks for keeping disposable files off the root filesystem."""
import os
import tempfile


def test_python_tempfile_uses_configured_runtime_directory():
    from config import Config

    runtime_tmp = os.path.realpath(Config.runtime_tmp_dir())
    assert os.path.realpath(tempfile.gettempdir()) == runtime_tmp
    assert os.path.realpath(os.environ['TMPDIR']) == runtime_tmp
    for cache_dir in (os.environ['NUMBA_CACHE_DIR'], os.environ['MPLCONFIGDIR']):
        assert os.path.commonpath([runtime_tmp, os.path.realpath(cache_dir)]) == runtime_tmp

    with tempfile.TemporaryDirectory(prefix='runtime_tmp_test_') as created:
        assert os.path.commonpath([runtime_tmp, os.path.realpath(created)]) == runtime_tmp
