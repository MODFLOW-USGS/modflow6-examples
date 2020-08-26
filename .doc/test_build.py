# Build files for sphinx.
import os
import sys
import shutil
from subprocess import Popen, PIPE
import flopy
import pymake

# -- download executables ----------------------------------------------------
pth = os.path.join("..", "bin")
if not os.path.isdir(pth):
    os.makedirs(pth)
pymake.getmfexes(pth, verbose=True)

# -- update mf6 executables with latest nightly build ------------------------
osname = sys.platform.lower()
if osname == "win32":
    key = "win64.zip"
elif osname == "darwin":
    key = "mac.zip"
elif osname == "linux":
    key = "linux.zip"
url = pymake.get_repo_assets("MODFLOW-USGS/modflow6-nightly-build")[key]
pymake.download_and_unzip(url, pth, verbose=True)

# -- update flopy classes ----------------------------------------------------
flopy.mf6.utils.generate_classes(branch="develop", backup=False)

# -- update notebooks --------------------------------------------------------
pth = os.path.join("..", "scripts", "process-scripts.py")
args = ("python", pth)
print(" ".join(args))
proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(pth))
stdout, stderr = proc.communicate()
if stdout:
    print(stdout.decode("utf-8"))
if stderr:
    print("Errors:\n{}".format(stderr.decode("utf-8")))

# -- get list of notebooks ---------------------------------------------------
pth = os.path.join("..", "notebooks")
nb_files = [os.path.join(pth, file_name) for file_name in sorted(os.listdir(pth)) if
            file_name.endswith(".ipynb")]

# -- run notebooks with nbconvert --------------------------------------------
output_pth = os.path.join("_notebooks")
if os.path.isdir(output_pth):
    shutil.rmtree(output_pth)
os.makedirs(output_pth)
for fpth in nb_files:
    args = ("jupyter",
            "nbconvert",
            "--ExecutePreprocessor.timeout=600",
            "--to",
            "notebook",
            "--execute",
            fpth,
            "--output-dir",
            output_pth,
            "--output",
            os.path.basename(fpth)
            )
    print(" ".join(args))
    os.system(" ".join(args))
