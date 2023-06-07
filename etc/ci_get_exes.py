import os
import sys
import numpy as np
import matplotlib as mpl
import pathlib as pl
import flopy

# Print python package versions
flopypth = flopy.__path__[0]
print("python version:     {}".format(sys.version))
print("numpy version:      {}".format(np.__version__))
print("matplotlib version: {}".format(mpl.__version__))
print("flopy version:      {}".format(flopy.__version__))
print("")
print("flopy is installed in:  {}".format(flopypth))

# Update flopy from GitHub repo
flopy.mf6.utils.generate_classes(branch="develop", backup=False)

# Set the directory for the binary files and make the
# directory if it does not exist
bin_pth = pl.Path("../bin")
bin_pth.mkdir(parents=True, exist_ok=True)

# Get executables
# get the latest executable release
os.system(f"get-modflow {bin_pth}")

# Replace MODFLOW 6 executables with the latest versions
os.system(f"get-modflow --repo modflow6-nightly-build {bin_pth}")
