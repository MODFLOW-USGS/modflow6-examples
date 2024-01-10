from pathlib import Path

from modflow_devtools.misc import run_cmd


PROJ_ROOT = Path(__file__).parent.parent


def test_examples(example):
    # make examples directory if it does not exist
    (PROJ_ROOT / "examples").mkdir(exist_ok=True)

    # run script
    run, script = example
    args = ["RUN=True"] if run else []
    args += ["python", script]
    stdout, stderr, returncode = run_cmd(*args, verbose=True)
    assert not returncode, stdout + stderr
