import os
import sys
import shutil
import pymake
from subprocess import Popen, PIPE

# change to root directory - local run only
starting_dir = os.getcwd()
os.chdir("..")

# Get executables
pth = os.path.join("bin")
if not os.path.isdir(pth):
    os.makedirs(pth)
pymake.getmfexes(pth, verbose=True)

# Replace MODFLOW 6 executables with the latest versions
osname = sys.platform.lower()
if osname == "win32":
    key = "win64.zip"
elif osname == "darwin":
    key = "mac.zip"
elif osname == "linux":
    key = "linux.zip"
url = pymake.get_repo_assets("MODFLOW-USGS/modflow6-nightly-build")[key]
pymake.download_and_unzip(url, "bin", verbose=True)
# local copy (put your own path here)
pth = "/Users/jdhughes/Documents/Development/modflow6_git/modflow6_fork/bin/"
files = [file_name for file_name in os.listdir(pth) if not file_name.endswith(".exe")]
for f in files:
    src = os.path.join(pth, f)
    dst = os.path.join("bin", f)
    shutil.copyfile(src, dst)

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
        ("pdflatex", texnam),
        ("bibtex", bibnam),
        ("pdflatex", texnam),
        ("pdflatex", texnam),
       )
for arg in args:
    print("running command...'{}'".format(" ".join(arg)))
    with Popen(arg, stdout=PIPE, stderr=PIPE, cwd=ws) as proc:
        stdout, stderr = proc.communicate(timeout=5)
        if stdout:
            print(stdout.decode("utf-8"))
        if stderr:
            print("Errors:\n{}".format(stderr.decode("utf-8")))

# remove zip file - local run only
os.remove("modflow6-examples.zip")

# change back to root directory - local run only
os.chdir(starting_dir)
