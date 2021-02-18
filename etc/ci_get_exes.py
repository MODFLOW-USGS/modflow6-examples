import os
import sys
import numpy as np
import matplotlib as mpl
import flopy
import pymake

# Print python package versions
flopypth = flopy.__path__[0]
pymakepth = pymake.__path__[0]
print("python version:     {}".format(sys.version))
print("numpy version:      {}".format(np.__version__))
print("matplotlib version: {}".format(mpl.__version__))
print("flopy version:      {}".format(flopy.__version__))
print("pymake version:     {}".format(pymake.__version__))
print("")
print("flopy is installed in:  {}".format(flopypth))
print("pymake is installed in: {}".format(pymakepth))

# Update flopy from GitHub repo
flopy.mf6.utils.generate_classes(branch="develop", backup=False)

# Get executables
bin_pth = os.path.join("..", "bin")
if not os.path.isdir(bin_pth):
    os.makedirs(bin_pth)
pymake.getmfexes(pth=bin_pth, verbose=True)

# Replace MODFLOW 6 executables with the latest versions
pymake.getmfnightly(pth=bin_pth, exes=["mf6", "libmf6"], verbose=True)
