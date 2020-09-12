# Build files for sphinx.
import os
import sys
import shutil
from subprocess import Popen, PIPE
import flopy
import pymake

# -- determine if running on CI or rtd
is_CI = 'CI' in os.environ or os.environ.get('READTHEDOCS') == 'True'

# # -- download executables ----------------------------------------------------
# pth = os.path.join("..", "bin")
# if not os.path.isdir(pth):
#     os.makedirs(pth)
# targets = pymake.usgs_program_data.get_keys(current=True)
# targets.remove("mf6")
# targets.remove("libmf6")
# targets.remove("zbud6")
# pymake.getmfexes(pth, verbose=True, exes=targets)
#
# # -- update mf6 executables with latest nightly build ------------------------
# if is_CI:
#     osname = sys.platform.lower()
#     if osname == "win32":
#         key = "win64.zip"
#     elif osname == "darwin":
#         key = "mac.zip"
#     elif osname == "linux":
#         key = "linux.zip"
#     url = pymake.get_repo_assets("MODFLOW-USGS/modflow6-nightly-build")[key]
#     pymake.download_and_unzip(url, pth, verbose=True)

# -- update flopy classes ----------------------------------------------------
flopy.mf6.utils.generate_classes(branch="develop", backup=False)

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
    py_files = [file_name for file_name in sorted(os.listdir(pth)) if
                file_name.endswith(".py") and file_name.startswith("ex-")]
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
nb_files = [file_name for file_name in sorted(os.listdir(pth)) if
            file_name.endswith(".ipynb") and file_name.startswith("ex-")]

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
print("clean up {}".format(copy_pth))
if os.path.isdir(copy_pth):
    shutil.rmtree(copy_pth)

# -- copy executed notebooks to ./_notebooks ---------------------------------
print("copy files in {} -> {}".format(dst_pth, copy_pth))
shutil.copytree(dst_pth, copy_pth)

# -- clean up (remove) dst_pth directory -------------------------------------
print("clean up {}".format(dst_pth))
shutil.rmtree(dst_pth)