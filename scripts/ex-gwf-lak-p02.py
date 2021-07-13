# ## Lake package (LAK) Package problem 2
#
# This is the lake package example problem (test 2) from the
# Lake Package documentation (Merritt and Konikow, 2000).
#

# ### LAK Package problem 2 Setup
#
# Imports

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import flopy
import shapefile as shp

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set figure properties specific to the

figure_size = (6.3, 5.6)
masked_values = (0, 1e30, -1e30)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-lak-p02"

# Model units

length_units = "feet"
time_units = "days"

# Table LAK Package problem 2 Model Parameters

nper = 1  # Number of periods
nlay = 5  # Number of layers
nrow = 27  # Number of rows
ncol = 17  # Number of columns
top = 200.0  # Top of the model ($ft$)
botm_str = "102., 97., 87., 77., 67."  # Bottom elevations ($ft$)
strt = 115.0  # Starting head ($ft$)
k11 = 30.0  # Horizontal hydraulic conductivity ($ft/d$)
k33 = 30.0  # Vertical hydraulic conductivity ($ft/d$)
ss = 3e-4  # Specific storage ($1/d$)
sy = 0.2  # Specific yield (unitless)
H1 = 160.0  # Constant head on left side of model ($ft$)
H2 = 140.0  # Constant head on right side of model ($ft$)
recharge = 0.0116  # Aereal recharge rate ($ft/d$)
etvrate = 0.0141  # Maximum evapotranspiration rate ($ft/d$)
etvdepth = 15.0  # Evapotranspiration extinction depth ($ft$)
lak_strt = 130.0  # Starting lake stage ($ft$)
lak_etrate = 0.0103  # Lake evaporation rate ($ft/d$)

# parse parameter strings into tuples

botm = [float(value) for value in botm_str.split(",")]

# Static temporal data used by TDIS file

tdis_ds = ((1500.0, 200, 1.005),)

# define delr and delc
delr = np.array(
    [
        250.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        500.0,
        500.0,
        500.0,
        500.0,
        500.00,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        250.0,
    ]
)
delc = np.array(
    [
        250.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        500.0,
        500.0,
        500.0,
        500.0,
        500.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        500.0,
        500.0,
        500.0,
        500.0,
        500.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        250.0,
    ]
)

# Define dimensions
extents = (0.0, delr.sum(), 0.0, delc.sum())
shape2d = (nrow, ncol)
shape3d = (nlay, nrow, ncol)

# Load the idomain arrays

data_pth = os.path.join("..", "data", sim_name)
fpth = os.path.join(data_pth, "idomain-01.txt")
idomain0 = np.loadtxt(fpth, dtype=int)
fpth = os.path.join(data_pth, "idomain-02.txt")
idomain1 = np.loadtxt(fpth, dtype=int)
idomain = [idomain0, idomain1, 1, 1, 1]

# create linearly varying evapotranspiration surface

xlen = delr.sum() - 0.5 * (delr[0] + delr[-1])
x = 0.0
s1d = H1 * np.ones(ncol, dtype=float)
for idx in range(1, ncol):
    x += 0.5 * (delr[idx - 1] + delr[idx])
    frac = x / xlen
    s1d[idx] = H1 + (H2 - H1) * frac
surf = np.tile(s1d, (nrow, 1))
surf[idomain0 == 0] = botm[0] - 2
surf[idomain1 == 0] = botm[1] - 2

# ### Create LAK Package problem 2 Model Boundary Conditions
#
# Constant head boundary conditions
#

chd_spd = []
for k in range(nlay):
    chd_spd += [[k, i, 0, H1] for i in range(nrow)]
    chd_spd += [[k, i, ncol - 1, H2] for i in range(nrow)]

# LAK Package

lak_time_conv = 86400.0
lak_len_conv = 3.28081

lak0_conn = [
    [0, 0, 0, 6, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 1, 0, 7, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 2, 0, 8, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 3, 0, 9, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 4, 0, 10, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 5, 0, 5, 6, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 6, 1, 6, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 7, 1, 7, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 8, 1, 7, 6, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 9, 1, 8, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 10, 1, 8, 6, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 11, 1, 9, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 12, 1, 9, 6, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 13, 1, 10, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 14, 0, 11, 6, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 15, 0, 5, 7, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 16, 1, 6, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 17, 1, 6, 7, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 18, 2, 7, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 19, 2, 8, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 20, 2, 9, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 21, 1, 10, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 22, 1, 10, 7, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 23, 0, 11, 7, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 24, 0, 5, 8, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 25, 1, 6, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 26, 1, 6, 8, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 27, 2, 7, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 28, 2, 8, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 29, 2, 9, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 30, 1, 10, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 31, 1, 10, 8, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 32, 0, 11, 8, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 33, 0, 5, 9, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 34, 1, 6, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 35, 1, 6, 9, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 36, 2, 7, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 37, 2, 8, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 38, 2, 9, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 39, 1, 10, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 40, 1, 10, 9, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 41, 0, 11, 9, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 42, 0, 5, 10, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 43, 1, 6, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 44, 1, 7, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 45, 1, 7, 10, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 46, 1, 8, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 47, 1, 8, 10, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 48, 1, 9, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 49, 1, 9, 10, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [0, 50, 1, 10, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [0, 51, 0, 11, 10, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 52, 0, 6, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 53, 0, 7, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 54, 0, 8, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 55, 0, 9, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [0, 56, 0, 10, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
]

lak1_conn = [
    [1, 0, 0, 16, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 1, 0, 17, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 2, 0, 18, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 3, 0, 19, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 4, 0, 20, 5, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 5, 0, 15, 6, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 6, 1, 16, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 7, 1, 17, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 8, 1, 17, 6, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 9, 1, 18, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 10, 1, 18, 6, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 11, 1, 19, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 12, 1, 19, 6, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 13, 1, 20, 6, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 14, 0, 21, 6, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 15, 0, 15, 7, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 16, 1, 16, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 17, 1, 16, 7, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 18, 2, 17, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 19, 2, 18, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 20, 2, 19, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 21, 1, 20, 7, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 22, 1, 20, 7, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 23, 0, 21, 7, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 24, 0, 15, 8, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 25, 1, 16, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 26, 1, 16, 8, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 27, 2, 17, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 28, 2, 18, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 29, 2, 19, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 30, 1, 20, 8, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 31, 1, 20, 8, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 32, 0, 21, 8, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 33, 0, 15, 9, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 34, 1, 16, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 35, 1, 16, 9, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 36, 2, 17, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 37, 2, 18, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 38, 2, 19, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 39, 1, 20, 9, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 40, 1, 20, 9, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 41, 0, 21, 9, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 42, 0, 15, 10, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 43, 1, 16, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 44, 1, 17, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 45, 1, 17, 10, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 46, 1, 18, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 47, 1, 18, 10, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 48, 1, 19, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 49, 1, 19, 10, "HORIZONTAL", 0.1, 0, 0, 250, 500],
    [1, 50, 1, 20, 10, "VERTICAL", 0.1, 0, 0, 0, 0],
    [1, 51, 0, 21, 10, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 52, 0, 16, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 53, 0, 17, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 54, 0, 18, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 55, 0, 19, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
    [1, 56, 0, 20, 11, "HORIZONTAL", 0.1, 0, 0, 500, 500],
]

lak_packagedata = [
    [0, lak_strt, len(lak0_conn)],
    [1, lak_strt, len(lak1_conn)],
]

lak_outlets = [
    [0, 0, -1, "manning", 114.85, 5.0, 0.05, 8.206324419006205e-4],
    [1, 1, -1, "manning", 109.4286, 5.0, 0.05, 9.458197164349258e-4],
]

lak_spd = [
    [0, "rainfall", recharge],
    [0, "evaporation", lak_etrate],
    [1, "rainfall", recharge],
    [1, "evaporation", lak_etrate],
]

# SFR package

sfr_conv = 128390.4

sfr_pakdata = [
    [
        0,
        0,
        1,
        4,
        1000,
        5,
        0.001103448,
        123.94827,
        0.5,
        0.5,
        0.050000001,
        1,
        1,
        0,
    ],
    [
        1,
        0,
        2,
        4,
        1000,
        5,
        0.001103448,
        122.84483,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        2,
        0,
        3,
        4,
        1000,
        5,
        0.001103448,
        121.74138,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        3,
        0,
        3,
        5,
        1000,
        5,
        0.001103448,
        120.63793,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        4,
        0,
        3,
        6,
        500,
        5,
        0.001103448,
        119.81035,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        5,
        0,
        3,
        7,
        750,
        5,
        0.001103448,
        119.12069,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        6,
        0,
        4,
        7,
        1000,
        5,
        0.001103448,
        118.15517,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        7,
        0,
        5,
        7,
        1000,
        5,
        0.001103448,
        117.05173,
        0.5,
        0.5,
        0.050000001,
        1,
        1,
        0,
    ],
    [
        8,
        0,
        11,
        8,
        1000,
        5,
        0.000820632,
        114.43968,
        0.5,
        0.5,
        0.050000001,
        1,
        1,
        0,
    ],
    [
        9,
        0,
        12,
        8,
        1000,
        5,
        0.000820632,
        113.61905,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        10,
        0,
        13,
        9,
        559,
        5,
        0.000820632,
        112.97937,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        11,
        0,
        13,
        9,
        559,
        5,
        0.000820632,
        112.52063,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        12,
        0,
        14,
        9,
        1000,
        5,
        0.000820632,
        111.88095,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        13,
        0,
        15,
        9,
        1000,
        5,
        0.000820632,
        111.06032,
        0.5,
        0.5,
        0.050000001,
        1,
        1,
        0,
    ],
    [
        14,
        0,
        21,
        9,
        1000,
        5,
        0.00094582,
        108.95569,
        0.5,
        0.5,
        0.050000001,
        1,
        1,
        0,
    ],
    [
        15,
        0,
        22,
        9,
        750,
        5,
        0.00094582,
        108.1281,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        16,
        0,
        22,
        10,
        500,
        5,
        0.00094582,
        107.53696,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        17,
        0,
        22,
        11,
        1000,
        5,
        0.00094582,
        106.82759,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        18,
        0,
        22,
        12,
        1000,
        5,
        0.00094582,
        105.88177,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        19,
        0,
        22,
        13,
        1000,
        5,
        0.00094582,
        104.93595,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        20,
        0,
        22,
        14,
        1000,
        5,
        0.00094582,
        103.99014,
        0.5,
        0.5,
        0.050000001,
        2,
        1,
        0,
    ],
    [
        21,
        0,
        22,
        15,
        1000,
        5,
        0.00094582,
        103.04431,
        0.5,
        0.5,
        0.050000001,
        1,
        1,
        0,
    ],
]

sfr_conn = [
    [0, -1],
    [1, 0, -2],
    [2, 1, -3],
    [3, 2, -4],
    [4, 3, -5],
    [5, 4, -6],
    [6, 5, -7],
    [7, 6],
    [8, -9],
    [9, 8, -10],
    [10, 9, -11],
    [11, 10, -12],
    [12, 11, -13],
    [13, 12],
    [14, -15],
    [15, 14, -16],
    [16, 15, -17],
    [17, 16, -18],
    [18, 17, -19],
    [19, 18, -20],
    [20, 19, -21],
    [21, 20],
]

sfr_spd = [[0, "inflow", 691200.0]]

# MVR package

mvr_paks = [
    ["SFR-1"],
    ["LAK-1"],
]

mvr_spd = [
    ["SFR-1", 7, "LAK-1", 0, "FACTOR", 1.0],
    ["LAK-1", 0, "SFR-1", 8, "FACTOR", 1.0],
    ["SFR-1", 13, "LAK-1", 1, "FACTOR", 1.0],
    ["LAK-1", 1, "SFR-1", 14, "FACTOR", 0.5],
]

# Solver parameters

nouter = 500
ninner = 100
hclose = 1e-9
rclose = 1e-6


# ### Functions to build, write, run, and plot the MODFLOW 6 LAK Package problem 2 model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model():
    if config.buildModel:
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
        )
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_ds, time_units=time_units
        )
        flopy.mf6.ModflowIms(
            sim,
            print_option="summary",
            linear_acceleration="bicgstab",
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord="{} strict".format(rclose),
        )
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=sim_name, newtonoptions="newton", save_flows=True
        )
        flopy.mf6.ModflowGwfdis(
            gwf,
            length_units=length_units,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            idomain=idomain,
            top=top,
            botm=botm,
        )
        obs_file = "{}.gwf.obs".format(sim_name)
        csv_file = obs_file + ".csv"
        obslist = [
            ["A", "head", (0, 3, 3)],
            ["B", "head", (0, 13, 8)],
            ["C", "head", (0, 23, 13)],
        ]
        obsdict = {csv_file: obslist}
        flopy.mf6.ModflowUtlobs(
            gwf, filename=obs_file, print_input=False, continuous=obsdict
        )

        flopy.mf6.ModflowGwfnpf(
            gwf,
            icelltype=1,
            k=k11,
            k33=k33,
            save_specific_discharge=True,
        )
        flopy.mf6.ModflowGwfsto(
            gwf,
            iconvert=1,
            sy=sy,
            ss=ss,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
        flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
        flopy.mf6.ModflowGwfevta(
            gwf, surface=surf, rate=etvrate, depth=etvdepth
        )
        lak = flopy.mf6.ModflowGwflak(
            gwf,
            pname="LAK-1",
            time_conversion=lak_time_conv,
            length_conversion=lak_len_conv,
            mover=True,
            print_stage=True,
            nlakes=2,
            noutlets=len(lak_outlets),
            packagedata=lak_packagedata,
            connectiondata=lak0_conn + lak1_conn,
            outlets=lak_outlets,
            perioddata=lak_spd,
        )
        obs_file = "{}.lak.obs".format(sim_name)
        csv_file = obs_file + ".csv"
        obs_dict = {
            csv_file: [
                ("lake1", "stage", (0,)),
                ("lake2", "stage", (1,)),
            ]
        }
        lak.obs.initialize(
            filename=obs_file, digits=10, print_input=True, continuous=obs_dict
        )
        flopy.mf6.ModflowGwfsfr(
            gwf,
            pname="SFR-1",
            unit_conversion=sfr_conv,
            mover=True,
            print_stage=True,
            print_flows=True,
            nreaches=len(sfr_pakdata),
            packagedata=sfr_pakdata,
            connectiondata=sfr_conn,
            perioddata=sfr_spd,
        )
        flopy.mf6.ModflowGwfmvr(
            gwf,
            maxmvr=4,
            maxpackages=2,
            packages=mvr_paks,
            perioddata=mvr_spd,
        )

        head_filerecord = "{}.hds".format(sim_name)
        budget_filerecord = "{}.cbc".format(sim_name)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        )
        return sim
    return None


# Function to write MODFLOW 6 LAK Package problem 2 model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the LAK Package problem 2 model.
# True is returned if the model runs successfully
#


@config.timeit
def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)
    return success


# Function to plot grid


def plot_grid(gwf, silent=True):
    sim_ws = os.path.join(ws, sim_name)

    # create lake array
    ilake = gwf.dis.idomain.array
    ilake[ilake == 0] = 100
    ilake[ilake == 1] = 0
    ilake[ilake == 100] = 1
    for k in range(nlay):
        for i in range(16, nrow):
            for j in range(ncol):
                if ilake[k, i, j] == 1:
                    ilake[k, i, j] = 2

    # get edges and centers of cells
    xedges, yedges = gwf.modelgrid.xyedges[0], gwf.modelgrid.xyedges[1]
    xcenters, ycenters = gwf.modelgrid.xycenters[0], gwf.modelgrid.xycenters[1]

    # create sfr network
    poly0 = [
        [xcenters[4], yedges[1]],
        [xcenters[4], ycenters[4]],
        [xcenters[7], ycenters[4]],
        [xcenters[7], yedges[6]],
    ]
    poly1 = [
        [xcenters[8], yedges[11]],
        [xcenters[8], yedges[13]],
        [xcenters[9], yedges[14]],
        [xcenters[9], yedges[16]],
    ]
    poly2 = [
        [xcenters[9], yedges[21]],
        [xcenters[9], ycenters[22]],
        [xedges[16], ycenters[22]],
    ]
    parts = [poly0, poly1, poly2]
    shape_pth = os.path.join(ws, sim_name, "sfr.shp")
    w = shp.Writer(target=shape_pth, shapeType=shp.POLYLINE)
    w.field("no", "C")
    w.line([poly0])
    w.record(["1"])
    w.line([poly1])
    w.record(["2"])
    w.line([poly2])
    w.record(["3"])
    w.close()
    sfr = shp.Reader(shape_pth)

    # load the observations
    fpth = os.path.join(ws, sim_name, "{}.lak.obs.csv".format(sim_name))
    lak_results = np.genfromtxt(fpth, delimiter=",", names=True)

    # create MODFLOW 6 head object
    file_name = gwf.oc.head_filerecord.get_data()[0][0]
    fpth = os.path.join(sim_ws, file_name)
    hobj = flopy.utils.HeadFile(fpth)

    # create MODFLOW 6 cell-by-cell budget object
    file_name = gwf.oc.budget_filerecord.get_data()[0][0]
    fpth = os.path.join(sim_ws, file_name)
    cobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    kstpkper = hobj.get_kstpkper()

    head = hobj.get_data(kstpkper=kstpkper[0])
    spdis = cobj.get_data(text="DATA-SPDIS", kstpkper=kstpkper[0])

    # add lake stage to heads
    head[ilake == 1] = lak_results["LAKE1"][-1]
    head[ilake == 2] = lak_results["LAKE2"][-1]

    # observation locations
    p1 = (xcenters[3], ycenters[3])
    p2 = (xcenters[8], ycenters[13])
    p3 = (xcenters[13], ycenters[23])

    # lake text locations
    pl1 = (xcenters[8], ycenters[8])
    pl2 = (xcenters[8], ycenters[18])

    fs = USGSFigure(figure_type="map", verbose=False)
    fig = plt.figure(
        figsize=(4, 6.9),
        tight_layout=True,
    )
    plt.axis("off")

    nrows, ncols = 10, 1
    axes = [fig.add_subplot(nrows, ncols, (1, 8))]

    for idx, ax in enumerate(axes):
        ax.set_xlim(extents[:2])
        ax.set_ylim(extents[2:])
        ax.set_aspect("equal")

    # legend axis
    axes.append(fig.add_subplot(nrows, ncols, (9, 10)))

    # set limits for legend area
    ax = axes[-1]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # get rid of ticks and spines for legend area
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.patch.set_alpha(0.0)

    ax = axes[0]
    mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
    mm.plot_bc("CHD", color="cyan")
    for shape in sfr.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        ax.plot(x, y, color="#3BB3D0", lw=1.5, zorder=1)
    mm.plot_inactive(color_noflow="#5DBB63")
    mm.plot_grid(lw=0.5, color="black")
    cv = mm.contour_array(
        head,
        levels=np.arange(120, 160, 5),
        linewidths=0.75,
        linestyles="-",
        colors="blue",
        masked_values=masked_values,
    )
    plt.clabel(cv, fmt="%1.0f")
    mm.plot_specific_discharge(spdis, normalize=True, color="0.75")
    ax.plot(p1[0], p1[1], marker="o", mfc="red", mec="black", ms=4)
    ax.plot(p2[0], p2[1], marker="o", mfc="red", mec="black", ms=4)
    ax.plot(p3[0], p3[1], marker="o", mfc="red", mec="black", ms=4)
    ax.set_xlabel("x-coordinate, in feet")
    ax.set_ylabel("y-coordinate, in feet")
    fs.add_text(
        ax,
        "A",
        x=p1[0] + 150,
        y=p1[1] + 150,
        transform=False,
        bold=False,
        color="red",
        ha="left",
        va="bottom",
    )
    fs.add_text(
        ax,
        "B",
        x=p2[0] + 150,
        y=p2[1] + 150,
        transform=False,
        bold=False,
        color="red",
        ha="left",
        va="bottom",
    )
    fs.add_text(
        ax,
        "C",
        x=p3[0] + 150,
        y=p3[1] + 150,
        transform=False,
        bold=False,
        color="red",
        ha="left",
        va="bottom",
    )
    fs.add_text(
        ax,
        "Lake 1",
        x=pl1[0],
        y=pl1[1],
        transform=False,
        italic=False,
        color="white",
        ha="center",
        va="center",
    )
    fs.add_text(
        ax,
        "Lake 2",
        x=pl2[0],
        y=pl2[1],
        transform=False,
        italic=False,
        color="white",
        ha="center",
        va="center",
    )
    fs.remove_edge_ticks(ax)

    # legend
    ax = axes[-1]
    ax.plot(
        -10000,
        -10000,
        lw=0,
        marker="s",
        ms=10,
        mfc="#5DBB63",
        mec="black",
        markeredgewidth=0.5,
        label="Lake boundary",
    )
    ax.plot(
        -10000,
        -10000,
        lw=0,
        marker="s",
        ms=10,
        mfc="cyan",
        mec="black",
        markeredgewidth=0.5,
        label="Constant-head boundary",
    )
    ax.plot(
        -10000,
        -10000,
        lw=1.5,
        color="#3BB3D0",
        label="Stream network",
    )
    ax.plot(
        -10000,
        -10000,
        lw=0,
        marker="o",
        ms=4,
        mfc="red",
        mec="black",
        markeredgewidth=0.5,
        label="Observation well",
    )
    ax.plot(
        -10000,
        -10000,
        lw=0.75,
        ls="-",
        color="blue",
        label=r"Head contour, $ft$",
    )
    ax.plot(
        -10000,
        -10000,
        lw=0,
        marker=u"$\u2192$",
        ms=10,
        mfc="0.75",
        mec="0.75",
        label="Normalized specific discharge",
    )
    fs.graph_legend(ax, loc="lower center", ncol=2)

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-grid{}".format(sim_name, config.figure_ext),
        )
        fig.savefig(fpth)

    return


# Function to plot the lake results


def plot_lak_results(gwf, silent=True):
    fs = USGSFigure(figure_type="graph", verbose=False)

    # load the observations
    fpth = os.path.join(ws, sim_name, "{}.lak.obs.csv".format(sim_name))
    lak_results = np.genfromtxt(fpth, delimiter=",", names=True)
    fpth = os.path.join(ws, sim_name, "{}.gwf.obs.csv".format(sim_name))
    gwf_results = np.genfromtxt(fpth, delimiter=",", names=True)

    dtype = [
        ("time", float),
        ("LAKE1", float),
        ("LAKE2", float),
        ("A", float),
        ("B", float),
        ("C", float),
    ]

    results = np.zeros((lak_results.shape[0] + 1), dtype=dtype)
    results["time"][1:] = lak_results["time"]
    results["LAKE1"][0] = lak_strt
    results["LAKE1"][1:] = lak_results["LAKE1"]
    results["LAKE2"][0] = lak_strt
    results["LAKE2"][1:] = lak_results["LAKE2"]
    results["A"][0] = strt
    results["A"][1:] = gwf_results["A"]
    results["B"][0] = strt
    results["B"][1:] = gwf_results["B"]
    results["C"][0] = strt
    results["C"][1:] = gwf_results["C"]

    # create the figure
    fig, axes = plt.subplots(
        ncols=1,
        nrows=2,
        sharex=True,
        figsize=(6.3, 4.3),
        constrained_layout=True,
    )

    ax = axes[0]
    ax.set_xlim(0, 1500)
    ax.set_ylim(110, 130)
    ax.plot(
        results["time"],
        results["LAKE1"],
        lw=0.75,
        ls="--",
        color="black",
        label="Lake 1 stage",
    )
    ax.plot(
        results["time"],
        results["LAKE2"],
        lw=0.75,
        ls="-.",
        color="black",
        label="Lake 2 stage",
    )
    ax.set_xticks([0, 250, 500, 750, 1000, 1250, 1500])
    ax.set_yticks([110, 115, 120, 125, 130])
    ax.set_ylabel("Lake stage, in feet")
    fs.graph_legend(ax, loc="upper right")
    fs.heading(ax, idx=0)

    ax = axes[1]
    ax.set_xlim(0, 1500)
    ax.set_ylim(110, 160)
    ax.plot(
        results["time"],
        results["A"],
        lw=0.75,
        ls="-",
        color="0.5",
        label="Point A",
    )
    ax.plot(
        results["time"],
        results["B"],
        lw=0.75,
        ls="-",
        color="black",
        label="Point B",
    )
    ax.plot(
        results["time"],
        results["C"],
        lw=0.75,
        ls="-.",
        color="black",
        label="Point C",
    )
    ax.set_xticks([0, 250, 500, 750, 1000, 1250, 1500])
    ax.set_xlabel("Simulation time, in days")
    ax.set_yticks([110, 120, 130, 140, 150, 160])
    ax.set_ylabel("Head, in feet")
    fs.graph_legend(ax, loc="upper left")
    fs.heading(ax, idx=1)

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-01{}".format(sim_name, config.figure_ext),
        )
        fig.savefig(fpth)

    return


# Function to plot the LAK Package problem 2 model results.


def plot_results(sim, silent=True):
    if config.plotModel:
        gwf = sim.get_model(sim_name)

        plot_grid(gwf, silent=silent)

        plot_lak_results(gwf, silent=silent)


# Function that wraps all of the steps for the LAK Package problem 2 model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(silent=True):
    sim = build_model()

    write_model(sim, silent=silent)

    success = run_model(sim, silent=silent)
    assert success, "could not run...{}".format(sim_name)

    if success:
        plot_results(sim, silent=silent)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(silent=False)


# nosetest end

if __name__ == "__main__":
    # ### LAK Package problem 2 Simulation
    #

    simulation()
