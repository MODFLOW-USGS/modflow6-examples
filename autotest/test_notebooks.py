from os import environ

import pytest
from modflow_devtools.markers import requires_exe
from modflow_devtools.misc import run_cmd, set_env

from conftest import NOTEBOOKS_PATH

EXCLUDE = ["ex-gwt-mt3dms-p10"]


@requires_exe("jupytext")
def test_notebooks(example):
    run, script = example
    if script.stem in EXCLUDE:
        pytest.skip(reason=f"excluded: {script.name}")
    with set_env(RUN=str(run)):
        args = [
            "jupytext",
            "--from",
            "py",
            "--to",
            "ipynb",
            "--execute",
            script,
            "-o",
            NOTEBOOKS_PATH / script.with_suffix(".ipynb").name,
        ]
        stdout, stderr, retcode = run_cmd(*args, verbose=True, env=environ)
        assert not retcode, stdout + stderr
