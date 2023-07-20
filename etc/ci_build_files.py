import os
import shutil
import sys
from subprocess import PIPE, Popen

import pytest


def get_command_args():
    run_arg, zip_arg = False, False
    for idx, arg in enumerate(sys.argv):
        if "--run=" in arg:
            run_arg = arg.split("=")[1]
        elif arg == "--run":
            run_arg = sys.argv[idx + 1]
        if "--zip=" in arg:
            zip_arg = arg.split("=")[1]
        elif arg == "--zip":
            zip_arg = sys.argv[idx + 1]

    return run_arg, zip_arg


@pytest.fixture(scope="session")
def run(pytestconfig):
    return pytestconfig.getoption("run")


def clean_files():
    pth = os.path.join("..", "examples")
    if os.path.isdir(pth):
        shutil.rmtree(pth)


def zip_files():
    print("zipping files in the 'examples' directory")
    sdir = os.getcwd()
    os.chdir("..")
    shutil.make_archive("modflow6-examples", "zip", "examples")
    os.chdir(sdir)


def test_build_run_all(build):
    # make examples directory if it does not exist
    os.makedirs(os.path.join("..", "examples"), exist_ok=True)

    # parse run boolean and script name and build command argument list
    run, script = build
    args = ["python", script]
    if not run:
        args += ["--no_run", "--no_plot"]

    # run the script
    print(" ".join(args))
    proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.join("..", "scripts"))
    stdout, stderr = proc.communicate()
    if stdout:
        print(stdout.decode("utf-8"))
    if stderr:
        print(stderr.decode("utf-8"))
    assert proc.returncode == 0, print(stderr.decode("utf-8"))


if __name__ == "__main__":
    run_arg, zip_arg = get_command_args()
    scripts_dict = {
        file_name: (run_arg, file_name)
        for file_name in sorted(os.listdir(os.path.join("..", "scripts")))
        if file_name.endswith(".py") and file_name.startswith("ex-")
    }

    # clean the files
    clean_files()

    # build and run, if necessary, run the scripts
    for build in scripts_dict.values():
        test_build_run_all(build)

    # zip up the input files, if required
    if zip_arg:
        zip_files()
