import os
import sys
import matplotlib.pyplot as plt
from IPython import get_ipython

# Setup working directories
work_directories = (
    os.path.join("..", "examples"),
    os.path.join("..", "figures"),
    os.path.join("..", "tables"),
)
for work_dir in work_directories:
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)

# run settings
buildModel = True
writeModel = True
runModel = True
plotModel = True
plotSave = True


# Test if being run as a script
def is_notebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False


# common figure settings
figure_ext = ".png"
plt.rcParams['image.cmap'] = "jet_r"

# parse command line arguments
if is_notebook():
    if "CI" in os.environ:
        plotSave = True
    else:
        plotSave = False
else:
    for idx, arg in enumerate(sys.argv):
        if arg in ("-nr", "--no_run"):
            runModel = False
        elif arg in ("-nw", "--no_write"):
            writeModel = False
        elif arg in ("-np", "--no_plot"):
            plotModel = False
        elif arg in ("-fe", "--figure_extension"):
            if idx + 1 < len(sys.argv):
                extension = sys.argv[idx + 1]
                if not extension.startswith("."):
                    extension = "." + extension
                figure_exts = tuple(plt.gcf().canvas.get_supported_filetypes().keys())
                if extension.lower() in figure_exts:
                    figure_ext = extension

# base example workspace
base_ws = os.path.join("..", "examples")
for idx, arg in enumerate(sys.argv):
    if arg in ("--destination"):
        if idx + 1 < len(sys.argv):
            base_ws = sys.argv[idx + 1]
            base_ws = os.path.abspath(base_ws)
assert os.path.isdir(base_ws)

# data files required for examples
data_ws = os.path.join("..", "data")

# set executable extension
eext = ""
if sys.platform.lower() == "win32":
    eext = ".exe"

# paths to executables
mf6_exe = os.path.abspath(os.path.join("..", "bin", "mf6" + eext))
mf2005_exe = os.path.abspath(os.path.join("..", "bin", "mf2005" + eext))
mf2005dbl_exe = os.path.abspath(os.path.join("..", "bin", "mf2005dbl" + eext))
mfnwt_exe = os.path.abspath(os.path.join("..", "bin", "mfnwt" + eext))
mt3dms_exe = os.path.abspath(os.path.join("..", "bin", "mt3dms" + eext))
mt3dusgs_exe = os.path.abspath(os.path.join("..", "bin", "mt3dusgs" + eext))   
