import os
import sys
import numpy as np
import matplotlib as mpl
import pathlib as pl

# Print python package versions
print("python version:     {}".format(sys.version))
print("numpy version:      {}".format(np.__version__))
print("matplotlib version: {}".format(mpl.__version__))
print("")

# Set the directory for the binary files and make the
# directory if it does not exist
bin_pth = pl.Path("../bin")
bin_pth.mkdir(parents=True, exist_ok=True)

# Get executables
# get the latest executable release
os.system(f"get-modflow {bin_pth}")

# Replace MODFLOW 6 executables with the latest versions
os.system(f"get-modflow --repo modflow6-nightly-build {bin_pth}")

# install triangle in python bin directory
os.system(f"get-modflow --subset 'triangle,' :python")
