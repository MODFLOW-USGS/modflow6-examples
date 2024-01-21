"""Test run example scripts."""

import sys
from os import environ

from modflow_devtools.misc import run_cmd, set_env


def test_scripts(example_script, write, run, plot, plot_show, plot_save, gif):
    with set_env(
        WRITE=str(write),
        RUN=str(run),
        PLOT=str(plot),
        PLOT_SHOW=str(plot_show),
        PLOT_SAVE=str(plot_save),
        GIF=str(gif),
    ):
        args = [sys.executable, example_script]
        stdout, stderr, retcode = run_cmd(*args, verbose=True, env=environ)
        assert not retcode, stdout + stderr
