import sys
from os import environ

from modflow_devtools.misc import run_cmd, set_env


def test_scripts(example_script, run, plot, save):
    with set_env(RUN=str(run), PLOT=str(plot), SAVE=str(save)):
        args = [sys.executable, example_script]
        stdout, stderr, retcode = run_cmd(*args, verbose=True, env=environ)
        assert not retcode, stdout + stderr
