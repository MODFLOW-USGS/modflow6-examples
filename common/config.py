import os
import sys
from IPython import get_ipython

# Setup working directories
work_directories = (os.path.join("..", "examples"),
                    os.path.join("..", "figures"),
                    os.path.join("..", "tables"),)
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
        if shell == 'ZMQInteractiveShell':
            return True  # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False


# parse command line arguments
if is_notebook():
    plotSave = False
else:
    for arg in sys.argv:
        if arg in ("-nr", "--no_run"):
            runModel = False
        elif arg in ("-nw", "--no_write"):
            writeModel = False
        elif arg in ("-np", "--no_plot"):
            plotModel = False

# common figure settings
figure_ext = ".png"

# paths to executables
mf6_exe = os.path.abspath(os.path.join("..", "bin", "mf6"))
