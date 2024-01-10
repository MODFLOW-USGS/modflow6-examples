import sys
from os import environ

from modflow_devtools.misc import run_cmd, set_env


def test_scripts(example):
    run, script = example
    with set_env(RUN=str(run)):
        args = [sys.executable, script]
        stdout, stderr, retcode = run_cmd(*args, verbose=True, env=environ)
        assert not retcode, stdout + stderr
