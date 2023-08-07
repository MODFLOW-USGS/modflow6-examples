# Build files for sphinx.
import os
import shutil
import sys
from subprocess import PIPE, Popen

import flopy

# -- determine if running on CI or rtd
is_CI = "CI" in os.environ or os.environ.get("READTHEDOCS") == "True"

# -- update flopy classes ----------------------------------------------------
flopy.mf6.utils.generate_classes(owner="aprovost-usgs", branch="PRT", backup=False)

# -- update notebooks and tables ---------------------------------------------
pth = os.path.join("..", "scripts")
args = ("python", "process-scripts.py")
print(" ".join(args))
proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
stdout, stderr = proc.communicate()
if stdout:
    print(stdout.decode("utf-8"))
if stderr:
    print("Errors:\n{}".format(stderr.decode("utf-8")))

# -- run the scripts ---------------------------------------------------------
if not is_CI:
    pth = os.path.join("..", "scripts")
    py_files = [
        file_name
        for file_name in sorted(os.listdir(pth))
        if file_name.endswith(".py") and file_name.startswith("ex-")
    ]
    for file_name in py_files:
        args = ("python", file_name)
        print(" ".join(args))
        proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
        stdout, stderr = proc.communicate()
        if stdout:
            print(stdout.decode("utf-8"))
        if stderr:
            print("Errors:\n{}".format(stderr.decode("utf-8")))

# -- create and edit examples.rst and copy the figures -----------------------
if not is_CI:
    pth = os.path.join("..", "etc")
    args = ("python", "ci_create_examples.py")
    print(" ".join(args))
    proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
    stdout, stderr = proc.communicate()
    if stdout:
        print(stdout.decode("utf-8"))
    if stderr:
        print("Errors:\n{}".format(stderr.decode("utf-8")))

# -- get list of notebooks ---------------------------------------------------
pth = os.path.join("..", "notebooks")
nb_files = [
    file_name
    for file_name in sorted(os.listdir(pth))
    if file_name.endswith(".ipynb") and file_name.startswith("ex-")
]

# -- run notebooks with jupytext ---------------------------------------------
src_pth = os.path.join("..", "notebooks")
dst_pth = os.path.join("..", ".nbrun")
if os.path.isdir(dst_pth):
    shutil.rmtree(dst_pth)
os.makedirs(dst_pth)
for file_name in nb_files:
    src = os.path.join(src_pth, file_name)
    dst = os.path.join(dst_pth, file_name)
    arg = (
        "jupytext",
        "--to ipynb",
        "--from ipynb",
        "--execute",
        "-o",
        dst,
        src,
    )
    print("running command...'{}'".format(" ".join(arg)))
    os.system(" ".join(arg))

# -- remove ./_notebooks if it exists ----------------------------------------
copy_pth = os.path.join("_notebooks")
print(f"clean up {copy_pth}")
if os.path.isdir(copy_pth):
    shutil.rmtree(copy_pth)

# -- copy executed notebooks to ./_notebooks ---------------------------------
print(f"copy files in {dst_pth} -> {copy_pth}")
shutil.copytree(dst_pth, copy_pth)

# -- clean up (remove) dst_pth directory -------------------------------------
print(f"clean up {dst_pth}")
shutil.rmtree(dst_pth)
