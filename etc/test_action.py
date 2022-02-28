import os
import sys
import shutil
import pymake
from subprocess import Popen, PIPE

# change to root directory - local run only
starting_dir = os.getcwd()
os.chdir("..")

# Get executables
args = ("python", "etc/ci_get_exes.py")
proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=".")
stdout, stderr = proc.communicate()
if stdout:
    print(stdout.decode("utf-8"))
if stderr:
    print("Errors:\n{}".format(stderr.decode("utf-8")))

# clean up examples directory - just for local runs
expth = os.path.join("examples")
examples = [os.path.join(expth, d) for d in os.listdir(expth) if
            os.path.isdir(os.path.join(expth, d))]
for e in examples:
    shutil.rmtree(e)

# run scripts without model runs
pth = os.path.join("scripts")
scripts = [file_name for file_name in os.listdir(pth) if
           file_name.endswith(".py") and file_name.startswith("ex-")]
for s in scripts:
    args = ("python", s, "--no_run", "--no_plot")
    proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
    stdout, stderr = proc.communicate()
    if stdout:
        print(stdout.decode("utf-8"))
    if stderr:
        print("Errors:\n{}".format(stderr.decode("utf-8")))

# zip up model input files
expth = os.path.join("examples")
zpth = "modflow6-examples"
shutil.make_archive("modflow6-examples", 'zip', expth)

# run scripts plus processing script
pth = os.path.join("scripts")
scripts = [file_name for file_name in os.listdir(pth) if
           file_name.endswith(".py") and file_name.startswith("ex-")]
scripts.append("process-scripts.py")
for s in scripts:
    args = ("python", s)
    proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
    stdout, stderr = proc.communicate()
    if stdout:
        print(stdout.decode("utf-8"))
    if stderr:
        print("Errors:\n{}".format(stderr.decode("utf-8")))

# build the LaTeX document
ws = "doc"
bibnam = "mf6examples"
texnam = bibnam + ".tex"
args = (
        ("latexmk", "-c", texnam),
        ("pdflatex", texnam),
        ("bibtex", bibnam),
        ("pdflatex", texnam),
        ("pdflatex", texnam),
       )
os.chdir(ws)
for arg in args:
    print("running command...'{}'".format(" ".join(arg)))
    os.system(" ".join(arg))
os.chdir("..")

# run the notebooks
src_pth = os.path.join("notebooks")
dst_pth = os.path.join(".notebooks")
if os.path.isdir(dst_pth):
    shutil.rmtree(dst_pth)
os.makedirs(dst_pth, exist_ok=True)
nb_files = [file_name for file_name in os.listdir(src_pth) if
            file_name.endswith(".ipynb") and file_name.startswith("ex-")]
for file_name in nb_files:
    src = os.path.join(src_pth, file_name)
    dst = os.path.join(dst_pth, file_name)
    arg = (
        "jupytext",
        "--from ipynb",
        "--to ipynb",
        "--execute",
        "--output",
        dst,
        src,
    )
    print(" ".join(arg))
    os.system(" ".join(arg))

# -- remove ./_notebooks if it exists ----------------------------------------
copy_pth = os.path.join(".doc", "_notebooks")
print("clean up {}".format(copy_pth))
if os.path.isdir(copy_pth):
    shutil.rmtree(copy_pth)

# -- copy executed notebooks to ./_notebooks ---------------------------------
print("copy files in {} -> {}".format(dst_pth, copy_pth))
shutil.copytree(dst_pth, copy_pth)

# -- clean up (remove) dst_pth directory -------------------------------------
print("clean up {}".format(dst_pth))
shutil.rmtree(dst_pth)

# remove zip file - local run only
os.remove("modflow6-examples.zip")

# change back to root directory - local run only
os.chdir(starting_dir)
