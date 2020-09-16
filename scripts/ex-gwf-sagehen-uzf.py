# ## Sagehen example with UZF, SFR, and MVR advanced packages activated
#
# This script reproduces example 1 in the UZF1 Techniques and Methods (Niswonger et al., 2006).

# ### MODFLOW 6 Sagehen Problem Setup

# Append to system path to include the common subdirectory

import os
import sys

sys.path.append(os.path.join("..", "common"))

# Imports

import flopy
import numpy as np
import pandas as pd
import config
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
from figspecs import USGSFigure
from datetime import datetime

sys.path.append(os.path.join("..", "data", "sagehen-uzf"))
import build_sagehen_uzf_helper_funcs as sageBld

mf6exe = os.path.abspath(config.mf6_exe)
assert os.path.isfile(mf6exe)
print(mf6exe)
exe_name_mf = config.mf2005_exe
exe_name_mt = config.mt3dms_exe

# Set figure properties specific to this problem

figure_size = (6, 6)
figure_size_ts = (6, 4)

# Base simulation and model name and workspace

ws = config.base_ws
example_name = "ex-gwf-sagehen"

# Model units

length_units = "meters"
time_units = "days"

# Table Sagehen Model Parameters

nlay = 1             # Number of layers in parent model
nrow = 73            # Number of rows in parent model
ncol = 81            # Number of columns in parent model
delr = 90.           # Parent model column width ($m$)
delc = 90.           # Parent model row width ($m$)
k11 = "0.005 - 0.3"  # Horizontal hydraulic conductivity ($m/d$)
k33 = "0.01 - 0.3"   # Vertical hydraulic conductivity ($m/d$)
sydum = "0.1 - 0.2"  # Specific Yield
dum1 = "0.15 - 0.25" # Saturated water content
dum2 = 0.10          # Extinction water content
dum3 = "0.005 - 0.3" # Saturated hydraulic conductivity of unsaturated zone
eps = 4.0            # Brooks-Corey Epsilon
sfrdum1 = 213        # Number of SFR stream reaches
sfrdum2 = 0.3        # Hydraulic conductivity of streambed ($m/d$)
sfrdum3 = 3          # Width of stream reaches ($m$)
sfrdum4 = 1          # Streambed thickness ($m$)


# Additional model input preparation

# Time related variables
num_ts = 399
perlen = [1] * num_ts
nper = len(perlen)
nstp = [1] * num_ts
tsmult = [1.0] * num_ts

# Further parent model grid discretization

# from mf-nwt .dis file
top = np.loadtxt(os.path.join("..","data","sagehen-uzf","dis_support","top1.txt"))
bot1 = np.loadtxt(os.path.join("..","data","sagehen-uzf","dis_support","bot1.txt"))
# from mf-nwt .bas file
idomain1 = np.loadtxt(os.path.join("..","data","sagehen-uzf","bas_support","ibnd1.txt"))
strt = np.loadtxt(os.path.join("..","data","sagehen-uzf","bas_support","strt1.txt"))
# peel out locations of negative values for setting constant head data
tmp1 = np.where(idomain1 < 0)
listOfChdCoords = list(zip(np.zeros_like(tmp1[0]), tmp1[0], tmp1[1]))
# get the corresponding constant head values
if(len(listOfChdCoords) > 0):
    chd_lay1 = list(np.take(strt , np.ravel_multi_index(tmp1, strt.shape)))
chdspd = []
for i in np.arange(len(listOfChdCoords)):
    chdspd.append([listOfChdCoords[i], chd_lay1[i]])
# finally, get rid of the negative values in idomain since mf6 treats negatives like zeros
idomain = np.abs(idomain1)

# from mf-nwt .upw file
k11 = np.loadtxt(os.path.join("..","data","sagehen-uzf","upw_support","kh1.txt"))
sy = np.loadtxt(os.path.join("..","data","sagehen-uzf","upw_support","sy1.txt"))
k33 = np.loadtxt(os.path.join("..","data","sagehen-uzf","upw_support","kv1.txt"))
icelltype = 1  # Water table resides in layer 1
iconvert = np.ones_like(strt)

# Solver settings

nouter, ninner = 300, 500
hclose, rclose, relax = 3e-2, 3e-2, 0.97

# #### Prepping input for SFR package 
# Package_data information

# Define the connections
conns = sageBld.gen_mf6_sfr_connections()

# These are zero based
sfrcells = [
    (0, 44, 8),
    (0, 43, 8),
    (0, 42, 8),
    (0, 42, 9),
    (0, 41, 9),
    (0, 41, 10),
    (0, 41, 11),
    (0, 40, 11),
    (0, 40, 12),
    (0, 39, 12),
    (0, 39, 13),
    (0, 38, 13),
    (0, 38, 14),
    (0, 38, 15),
    (0, 38, 16),
    (0, 37, 16),
    (0, 37, 17),
    (0, 37, 18),
    (0, 37, 19),
    (0, 37, 20),
    (0, 37, 21),
    (0, 37, 22),
    (0, 36, 23),
    (0, 36, 24),
    (0, 35, 24),
    (0, 35, 25),
    (0, 34, 25),
    (0, 34, 26),
    (0, 35, 27),
    (0, 35, 28),
    (0, 35, 29),
    (0, 35, 30),
    (0, 35, 31),
    (0, 34, 31),
    (0, 34, 32),
    (0, 33, 32),
    (0, 33, 33),
    (0, 32, 33),
    (0, 32, 34),
    (0, 31, 34),
    (0, 31, 35),
    (0, 31, 36),
    (0, 52, 36),
    (0, 51, 36),
    (0, 51, 37),
    (0, 50, 37),
    (0, 49, 37),
    (0, 49, 38),
    (0, 48, 38),
    (0, 47, 38),
    (0, 46, 38),
    (0, 46, 39),
    (0, 45, 39),
    (0, 44, 39),
    (0, 44, 38),
    (0, 43, 38),
    (0, 42, 37),
    (0, 41, 37),
    (0, 40, 37),
    (0, 39, 37),
    (0, 39, 38),
    (0, 38, 38),
    (0, 37, 38),
    (0, 36, 39),
    (0, 35, 39),
    (0, 35, 40),
    (0, 34, 40),
    (0, 33, 40),
    (0, 32, 40),
    (0, 31, 40),
    (0, 30, 41),
    (0, 30, 32),
    (0, 30, 33),
    (0, 30, 34),
    (0, 30, 35),
    (0, 30, 36),
    (0, 47, 47),
    (0, 46, 47),
    (0, 45, 47),
    (0, 45, 46),
    (0, 44, 46),
    (0, 43, 46),
    (0, 42, 46),
    (0, 41, 46),
    (0, 40, 46),
    (0, 40, 47),
    (0, 39, 47),
    (0, 38, 47),
    (0, 37, 46),
    (0, 36, 46),
    (0, 35, 47),
    (0, 34, 47),
    (0, 34, 48),
    (0, 33, 48),
    (0, 33, 49),
    (0, 32, 49),
    (0, 54, 71),
    (0, 53, 71),
    (0, 52, 71),
    (0, 51, 71),
    (0, 50, 71),
    (0, 49, 72),
    (0, 48, 72),
    (0, 47, 72),
    (0, 47, 73),
    (0, 46, 73),
    (0, 45, 74),
    (0, 44, 74),
    (0, 44, 75),
    (0, 43, 75),
    (0, 44, 61),
    (0, 43, 61),
    (0, 42, 61),
    (0, 42, 62),
    (0, 41, 62),
    (0, 40, 62),
    (0, 39, 62),
    (0, 23, 54),
    (0, 24, 54),
    (0, 24, 55),
    (0, 25, 55),
    (0, 26, 55),
    (0, 27, 56),
    (0, 28, 56),
    (0, 29, 56),
    (0, 30, 56),
    (0, 31, 56),
    (0, 32, 56),
    (0, 32, 57),
    (0, 33, 57),
    (0, 33, 58),
    (0, 34, 58),
    (0, 35, 58),
    (0, 36, 59),
    (0, 22, 70),
    (0, 23, 70),
    (0, 24, 70),
    (0, 25, 70),
    (0, 26, 71),
    (0, 26, 72),
    (0, 27, 72),
    (0, 28, 72),
    (0, 29, 72),
    (0, 30, 72),
    (0, 31, 72),
    (0, 32, 72),
    (0, 33, 72),
    (0, 33, 73),
    (0, 34, 73),
    (0, 35, 73),
    (0, 35, 72),
    (0, 36, 72),
    (0, 37, 71),
    (0, 38, 71),
    (0, 39, 71),
    (0, 40, 71),
    (0, 41, 71),
    (0, 41, 72),
    (0, 30, 37),
    (0, 30, 38),
    (0, 30, 39),
    (0, 30, 40),
    (0, 30, 41),
    (0, 29, 41),
    (0, 29, 42),
    (0, 29, 43),
    (0, 28, 43),
    (0, 28, 44),
    (0, 28, 45),
    (0, 28, 46),
    (0, 29, 46),
    (0, 29, 47),
    (0, 30, 48),
    (0, 31, 49),
    (0, 31, 50),
    (0, 32, 51),
    (0, 32, 52),
    (0, 33, 52),
    (0, 33, 53),
    (0, 34, 53),
    (0, 34, 54),
    (0, 34, 55),
    (0, 35, 56),
    (0, 35, 57),
    (0, 35, 58),
    (0, 36, 58),
    (0, 36, 59),
    (0, 37, 59),
    (0, 37, 60),
    (0, 37, 61),
    (0, 37, 62),
    (0, 38, 62),
    (0, 38, 63),
    (0, 38, 64),
    (0, 39, 64),
    (0, 39, 65),
    (0, 39, 66),
    (0, 39, 67),
    (0, 40, 68),
    (0, 40, 69),
    (0, 41, 70),
    (0, 41, 71),
    (0, 41, 72),
    (0, 41, 72),
    (0, 42, 72),
    (0, 42, 73),
    (0, 42, 74),
    (0, 43, 74),
    (0, 43, 75),
    (0, 43, 76),
    (0, 43, 77),
    (0, 43, 78),
    (0, 44, 78)
]

rlen = [
     90.,
     90.,
     75.,
     75.,
     75.,
     90.,
     75.,
     75.,
     75.,
     75.,
     75.,
     75.,
     90.,
     90.,
     60.,
     30.,
    102.,
     90.,
     90.,
     90.,
    102.,
    102.,
    102.,
     72.,
     30.,
     72.,
     30.,
     90.,
    102.,
     90.,
     90.,
    102.,
     30.,
     72.,
     30.,
     72.,
     30.,
     60.,
     72.,
     30.,
    102.,
     90.,
     90.,
     72.,
     30.,
    102.,
     60.,
     30.,
     90.,
    102.,
     60.,
     30.,
     90.,
     30.,
     60.,
    102.,
    102.,
     90.,
    102.,
     30.,
     60.,
    102.,
    102.,
    102.,
     60.,
     30.,
    102.,
     90.,
     90.,
    102.,
    114.,
     90.,
     90.,
    102.,
    102.,
     90.,
     90.,
    102.,
     30.,
     60.,
     90.,
    102.,
     90.,
     90.,
     60.,
     30.,
     90.,
     90.,
     90.,
     90.,
    102.,
     30.,
     60.,
     72.,
     30.,
    114.,
     90.,
     90.,
     90.,
    102.,
    102.,
    102.,
    102.,
     30.,
     60.,
    102.,
    102.,
     30.,
     72.,
     30.,
     90.,
    102.,
     30.,
     60.,
    102.,
     90.,
    102.,
     60.,
     30.,
     60.,
    102.,
     90.,
     90.,
     90.,
     90.,
     90.,
    102.,
     72.,
     30.,
     72.,
     30.,
    102.,
     90.,
    114.,
     90.,
    102.,
     90.,
    102.,
    102.,
     30.,
    102.,
     90.,
     90.,
     90.,
     90.,
     90.,
     30.,
     60.,
     90.,
     30.,
     60.,
    102.,
     90.,
     90.,
     90.,
     90.,
     30.,
     30.,
     90.,
     90.,
     90.,
     90.,
     72.,
     30.,
     90.,
     60.,
     30.,
     90.,
     90.,
     30.,
     60.,
    102.,
    114.,
    114.,
     90.,
    102.,
     30.,
     72.,
     60.,
     30.,
    102.,
    102.,
    114.,
     90.,
     60.,
     30.,
     72.,
     30.,
    102.,
    102.,
     30.,
     60.,
    120.,
     60.,
     30.,
     90.,
     90.,
     90.,
     90.,
    114.,
    102.,
     90.,
     30.,
     30.,
     30.,
    102.,
     60.,
     30.,
     90.,
     90.,
     90.,
     60.,
     30.
]

rwid = 3.0
rgrd = [
    0.150,
    0.120,
    0.180,
    0.160,
    0.130,
    0.111,
    0.047,
    0.060,
    0.040,
    0.100,
    0.227,
    0.090,
    0.042,
    0.064,
    0.083,
    0.081,
    0.062,
    0.065,
    0.089,
    0.097,
    0.141,
    0.186,
    0.217,
    0.203,
    0.206,
    0.206,
    0.144,
    0.147,
    0.135,
    0.129,
    0.124,
    0.136,
    0.162,
    0.147,
    0.157,
    0.147,
    0.115,
    0.117,
    0.111,
    0.120,
    0.099,
    0.073,
    0.037,
    0.038,
    0.060,
    0.048,
    0.024,
    0.029,
    0.032,
    0.028,
    0.024,
    0.029,
    0.033,
    0.038,
    0.032,
    0.038,
    0.051,
    0.047,
    0.037,
    0.063,
    0.063,
    0.049,
    0.069,
    0.077,
    0.063,
    0.045,
    0.037,
    0.043,
    0.048,
    0.054,
    0.065,
    0.067,
    0.091,
    0.091,
    0.071,
    0.073,
    0.021,
    0.031,
    0.045,
    0.033,
    0.029,
    0.042,
    0.075,
    0.103,
    0.092,
    0.095,
    0.087,
    0.083,
    0.094,
    0.102,
    0.093,
    0.081,
    0.099,
    0.077,
    0.057,
    0.056,
    0.044,
    0.050,
    0.075,
    0.076,
    0.074,
    0.074,
    0.071,
    0.072,
    0.056,
    0.060,
    0.048,
    0.043,
    0.049,
    0.039,
    0.042,
    0.056,
    0.081,
    0.071,
    0.068,
    0.068,
    0.068,
    0.044,
    0.078,
    0.071,
    0.051,
    0.054,
    0.056,
    0.056,
    0.050,
    0.048,
    0.038,
    0.022,
    0.049,
    0.059,
    0.043,
    0.043,
    0.045,
    0.049,
    0.042,
    0.031,
    0.016,
    0.010,
    0.012,
    0.015,
    0.012,
    0.011,
    0.022,
    0.044,
    0.056,
    0.060,
    0.114,
    0.100,
    0.067,
    0.086,
    0.127,
    0.141,
    0.118,
    0.100,
    0.083,
    0.087,
    0.100,
    0.067,
    0.056,
    0.083,
    0.100,
    0.076,
    0.045,
    0.020,
    0.053,
    0.042,
    0.038,
    0.047,
    0.047,
    0.057,
    0.040,
    0.032,
    0.045,
    0.053,
    0.042,
    0.049,
    0.094,
    0.085,
    0.036,
    0.027,
    0.030,
    0.033,
    0.024,
    0.017,
    0.025,
    0.021,
    0.015,
    0.010,
    0.010,
    0.012,
    0.018,
    0.022,
    0.017,
    0.019,
    0.010,
    0.013,
    0.022,
    0.017,
    0.021,
    0.043,
    0.044,
    0.038,
    0.050,
    0.033,
    0.021,
    0.020,
    0.024,
    0.029,
    0.020,
    0.011,
    0.024,
    0.033,
    0.022
]

rtp = [
    2458.0,
    2449.0,
    2337.0,
    2416.0,
    2413.0,
    2397.0,
    2393.0,
    2390.0,
    2386.0,
    2384.0,
    2367.0,
    2360.0,
    2355.0,
    2351.0,
    2344.0,
    2341.0,
    2335.0,
    2331.0,
    2323.0,
    2315.0,
    2305.0,
    2287.0,
    2267.0,
    2246.0,
    2239.0,
    2225.0,
    2218.0,
    2209.0,
    2195.0,
    2183.0,
    2171.0,
    2160.0,
    2149.0,
    2141.0,
    2134.0,
    2125.0,
    2119.0,
    2114.0,
    2106.0,
    2101.0,
    2092.0,
    2085.0,
    2146.0,
    2143.0,
    2141.0,
    2136.0,
    2134.0,
    2133.0,
    2131.0,
    2128.0,
    2126.0,
    2125.0,
    2123.0,
    2121.0,
    2119.0,
    2117.0,
    2112.0,
    2107.0,
    2103.0,
    2101.0,
    2096.0,
    2093.0,
    2087.0,
    2079.0,
    2073.0,
    2071.0,
    2068.0,
    2065.0,
    2060.0,
    2056.0,
    2048.0,
    2115.0,
    2109.0,
    2098.0,
    2091.0,
    2084.0,
    2112.0,
    2110.0,
    2108.0,
    2106.0,
    2104.0,
    2101.0,
    2096.0,
    2087.0,
    2079.0,
    2076.0,
    2069.0,
    2063.0,
    2054.0,
    2046.0,
    2035.0,
    2031.0,
    2026.0,
    2020.0,
    2017.0,
    2013.0,
    2000.0,
    1996.0,
    1991.0,
    1982.0,
    1976.0,
    1967.0,
    1961.0,
    1955.0,
    1953.0,
    1948.0,
    1942.0,
    1940.0,
    1937.0,
    1935.0,
    2003.0,
    1999.0,
    1994.0,
    1990.0,
    1985.0,
    1978.0,
    1972.0,
    2032.0,
    2030.0,
    2025.0,
    2021.0,
    2016.0,
    2011.0,
    2006.0,
    2001.0,
    1997.0,
    1992.0,
    1990.0,
    1989.0,
    1985.0,
    1983.0,
    1980.0,
    1976.0,
    1971.0,
    2051.0,
    2047.0,
    2045.0,
    2044.0,
    2043.0,
    2042.0,
    2041.0,
    2040.0,
    2039.0,
    2036.0,
    2031.0,
    2026.0,
    2022.0,
    2014.0,
    2010.0,
    2005.0,
    2001.0,
    1989.0,
    1976.0,
    1967.0,
    1958.0,
    1952.0,
    1945.0,
    1943.0,
    2076.0,
    2071.0,
    2061.0,
    2053.0,
    2048.0,
    2047.0,
    2041.0,
    2037.0,
    2036.0,
    2033.0,
    2029.0,
    2026.0,
    2023.0,
    2021.0,
    2017.0,
    2011.0,
    2006.0,
    2002.0,
    1998.0,
    1991.0,
    1988.0,
    1987.0,
    1985.0,
    1982.0,
    1978.0,
    1977.0,
    1975.0,
    1974.0,
    1973.0,
    1972.0,
    1970.0,
    1969.0,
    1968.0,
    1967.0,
    1966.0,
    1965.0,
    1964.0,
    1963.0,
    1961.0,
    1959.0,
    1958.0,
    1955.0,
    1949.0,
    1946.0,
    1943.0,
    1942.0,
    1941.0,
    1940.0,
    1938.0,
    1937.0,
    1935.0,
    1934.0,
    1933.0,
    1930.0,
    1929.0
]

rbth = 1.0
rhk = 0.6
man = 0.04
ustrf = 1.0
ndv = 0
pkdat = []
for i in np.arange(len(rlen)):
    ncon = len(conns[i]) - 1
    pkdat.append(
        (
            i,
            sfrcells[i],
            rlen[i],
            rwid,
            rgrd[i],
            rtp[i],
            rbth,
            rhk,
            man,
            ncon,
            ustrf,
            ndv,
        )
    )

# #### Prepping input for DRN package

# Instantiating MODFLOW 6 drain package
# Here, the drain (DRN) package is used to simulate groundwater discharge to land surface to keep this
# water separate from rejected infiltrated simulated by the UZF package. 
# Need to cycle through all land surface cells and create a drain for handling groundwater discharge to land surface
drn_spd  = []
drn_dict = {}
drn_dict_rev = {}
cond     = 10000  # Use an arbitrarily high conductance term to avoid impeding groundwater discharge
ddrn     = -0.25   # See definition of auxdepthname in drain package documentation to learn more about this parameter
idrnno   = 0
for i in np.arange(0, top.shape[0]):
    for j in np.arange(0, top.shape[1]):
        # Don't add drains to sfr and chd cells:
        sfrCell_bool = 1 if len([itm for itm in sfrcells if itm[1]==i and itm[2]==j]) > 0 else 0
        chdCell_bool = 1 if len([itm for itm in listOfChdCoords if itm[1]==i and itm[2]==j]) > 0 else 0
        if idomain1[i, j] and not sfrCell_bool and not chdCell_bool:
            drn_spd.append([(0, i, j), top[i, j], cond, ddrn])  #  'ddrn',
            # append dictionary of drain indices
            drn_dict.update({(i,j): idrnno})
            drn_dict_rev.update({idrnno: (i, j)})
            idrnno += 1


# #### Prepping input for UZF package 
# Package_data information

iuzbnd = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","iuzbnd.txt"))
thts = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","thts.txt"))
uzk33 = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","vks.txt"))
finf_grad = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","finf_gradient.txt"))
# next, load time series of multipliers
extwc_ts = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","extwc.txt"))
finf_ts = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","finf_ts.txt"))
pet_ts = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","pet_ts.txt"))
rootdp_ts = np.loadtxt(os.path.join("..","data","sagehen-uzf","uzf_support","rooting_depth.txt"))

# Need to set iuzbnd inactive where there are constant head cells, or where the
# model grid is inactive
cells2inactivate = idomain1 - iuzbnd
iuzbnd = iuzbnd + cells2inactivate
for chd_cell in listOfChdCoords:
    iuzbnd[chd_cell[1], chd_cell[2]] = 0

ha = 0.
hroot = 0.
rootact = 0.

uzf_packagedata = []
pd0             = []
iuzno_cell_dict = {}
iuzno_dict_rev  = {}
iuzno           = 0
surfdep         = 0.1
# Set up the UZF static variables
nuzfcells = 0
for k in range(nlay):
    for i in range(0, iuzbnd.shape[0] - 1):
        for j in range(0,iuzbnd.shape[1] - 1):
            if iuzbnd[i, j] != 0:
                nuzfcells += 1
                if k == 0:
                    lflag = 1
                    iuzno_cell_dict.update({(i, j): iuzno})  # establish new dictionary entry for current cell 
                                                             # addresses & iuzno connections are both 0-based
                    iuzno_dict_rev.update({iuzno: (i, j)})   # For post-processing the mvr output, need a dict with iuzno as key
                else:
                    lflag = 0
                
                # Set the vertical connection, which is the cell below, 
                # but in this 1 layer model set to -1 which flopy adjusts to 0 (no vert conn)
                ivertcon = -1
                
                vks = uzk33[i, j]
                thtr = extwc_ts[0]
                thtsx = thts[i, j]
                thti = thtr + 0.01
                eps = 4.0
                
                # Set the boundname for the land surface cells
                bndnm = 'sage'
                
                # <iuzno> <cellid(ncelldim)> <landflag> <ivertcon> <surfdep> <vks> <thtr> <thts> <thti> <eps> [<boundname>]
                uz = [iuzno,      (k, i, j),     lflag,  ivertcon,  surfdep,  vks,  thtr,  thtsx,  thti,  eps,   bndnm]
                uzf_packagedata.append(uz)
                
                iuzno += 1

# Next prepare the stress period data for UZF
# Store the steady state uzf stresses in dictionary
uzf_perioddata = {}
for t in range(num_ts):
    iuzno = 0
    spdx = []
    for i in range(0, iuzbnd.shape[0] - 1):
        for j in range(0, iuzbnd.shape[1] - 1):
            if iuzbnd[i, j] != 0:
                finf = finf_grad[i, j] * finf_ts[t]
                pet = pet_ts[t]
                extdp = rootdp_ts[t]
                extwc = extwc_ts[t]
                #mf6io.pdf: <iuzno> <finf> <pet> <extdp> <extwc> <ha> <hroot> <rootact>
                spdx.append([iuzno,  finf,  pet,  extdp,  extwc,  ha,  hroot,  rootact])
                iuzno += 1
    uzf_perioddata.update({t: spdx})

# Set up runoff connections, which relies on a helper function inside a companion script

# This function uses the top elevation array and SFR locations to calculate the irunbnd array from the UZF1 package.
# Of course in MF6, MVR now handles these runoff connections  
pth = os.path.join("..","data","sagehen-uzf")
irunbnd = sageBld.determine_runoff_conns_4mvr(pth, nrow, ncol)  # at this point, irunbnd is 1-based, compensate below

iuzno = 0
k     = 0             # Hard-wire the layer no.
first0ok = True
static_mvrperioddata = []
for i in range(0, iuzbnd.shape[0]):
    for j in range(0,iuzbnd.shape[1]):
        if irunbnd[i, j] > 0:           # This is a uzf -> sfr connection
            iuzno = iuzno_cell_dict.get((i, j))
            if iuzno or first0ok:
                static_mvrperioddata.append(('UZF-1', iuzno, 'SFR-1', irunbnd[i, j] - 1, 'FACTOR', 1.))
        
        drn_idx = drn_dict.get((i, j))
        if drn_idx:
            static_mvrperioddata.append(('DRN-1', drn_idx, 'SFR-1', irunbnd[i, j] - 1, 'FACTOR', 1.))
            first0ok = False

mvrspd = {0: static_mvrperioddata}           
mvrpack = [['UZF-1'], ['SFR-1'], ['DRN-1']]
maxpackages = len(mvrpack)
maxmvr = 10000    # Something arbitrarily high

# ### Function to build models
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model

def build_model(sim_name, silent=False):
    if config.buildModel:
        
        # Instantiate the MODFLOW 6 simulation
        sim_ws = os.path.join(ws, example_name)
        sim = flopy.mf6.MFSimulation(
            sim_name=example_name,
            version="mf6",
            sim_ws=sim_ws,
            exe_name=mf6exe
        )
        
        # Instantiating MODFLOW 6 time discretization
        tdis_rc = []
        for i in range(len(perlen)):
            tdis_rc.append((perlen[i], nstp[i], tsmult[i]))
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_rc, time_units=time_units
        )
        
        # Instantiating MODFLOW 6 groundwater flow model
        gwfname = example_name
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=gwfname,
            save_flows=True,
            newtonoptions=True,
            model_nam_file="{}.nam".format(gwfname),
        )
        
        # Instantiating MODFLOW 6 solver for flow model
        imsgwf = flopy.mf6.ModflowIms(
            sim,
            print_option="summary",
            complexity="complex",
            outer_dvclose=hclose,
            outer_maximum=nouter,
            under_relaxation="dbd",
            linear_acceleration="BICGSTAB",
            under_relaxation_theta=0.7,
            under_relaxation_kappa=0.08,
            under_relaxation_gamma=0.05,
            under_relaxation_momentum=0.0,
            backtracking_number=20,
            backtracking_tolerance=2.0,
            backtracking_reduction_factor=0.2,
            backtracking_residual_limit=5.0e-4,
            inner_dvclose=hclose,
            rcloserecord=[0.0001, "relative_rclose"],
            inner_maximum=ninner,
            relaxation_factor=relax,
            number_orthogonalizations=2,
            preconditioner_levels=8,
            preconditioner_drop_tolerance=0.001,
            filename="{}.ims".format(gwfname)
        )
        sim.register_ims_package(imsgwf, [gwf.name])
        
        # Instantiating MODFLOW 6 discretization package
        flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=bot1,
            idomain=idomain,
            filename="{}.dis".format(gwfname)
        )
        
        # Instantiating MODFLOW 6 initial conditions package for flow model
        flopy.mf6.ModflowGwfic(
            gwf, 
            strt=strt, 
            filename="{}.ic".format(gwfname)
        )
        
        # Instantiating MODFLOW 6 node-property flow package
        flopy.mf6.ModflowGwfnpf(
            gwf,
            save_flows=False,
            alternative_cell_averaging="AMT-HMK",
            icelltype=icelltype,
            k=k11,
            k33=k33,
            save_specific_discharge=False,
            filename="{}.npf".format(gwfname)
        )
        
        # Instantiate MODFLOW 6 storage package 
        flopy.mf6.ModflowGwfsto(
            gwf, 
            ss=2e-6, 
            sy=sy,
            iconvert=iconvert,
            steady_state={0:True},
            transient={1:True},
            filename='{}.sto'.format(gwfname)
        )
        
        # Instantiating MODFLOW 6 output control package for flow model
        flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord="{}.bud".format(gwfname),
            head_filerecord="{}.hds".format(gwfname),
            headprintrecord=[
                ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
            ],
            saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
            printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        )
        
        # Instantiating MODFLOW 6 constant head package
        chdspdx = {0: chdspd}
        flopy.mf6.ModflowGwfchd(
            gwf,
            maxbound=len(chdspd),
            stress_period_data=chdspdx,
            save_flows=False,
            pname="CHD-1",
            filename="{}.chd".format(gwfname),
        )
        
        maxbound = len(drn_spd)   # The total number 
        spd = {0: drn_spd}
        flopy.mf6.ModflowGwfdrn(
            gwf, 
            pname='DRN-1',
            auxiliary=['ddrn'],
            auxdepthname='ddrn',
            print_input=False, 
            print_flows=False,
            maxbound=maxbound,
            mover=True,
            stress_period_data=spd,   # wel_spd established in the MVR setup
            boundnames=False, 
            save_flows=True,
            filename='{}.drn'.format(gwfname)
        )
        
        # Instantiating MODFLOW 6 streamflow routing package
        flopy.mf6.ModflowGwfsfr(
            gwf,
            print_stage=False,
            print_flows=False,
            budget_filerecord=gwfname + ".sfr.bud",
            save_flows=True,
            mover=True,
            pname="SFR-1",
            unit_conversion=86400.0,
            boundnames=True,
            nreaches=len(conns),
            packagedata=pkdat,
            connectiondata=conns,
            perioddata=None,
            filename="{}.sfr".format(gwfname),
        )
        
        # Instantiating MODFLOW 6 unsaturated zone flow package
        flopy.mf6.ModflowGwfuzf(
            gwf, 
            nuzfcells=nuzfcells, 
            boundnames=True,
            mover=True,
            ntrailwaves=15, 
            nwavesets=150, 
            print_flows=False,
            save_flows=True,
            simulate_et=True, 
            linear_gwet=True,
            packagedata=uzf_packagedata, 
            perioddata=uzf_perioddata,
            budget_filerecord='{}.uzf.bud'.format(gwfname),
            pname='UZF-1',
            filename='{}.uzf'.format(gwfname)
        )
        
        flopy.mf6.ModflowGwfmvr(
            gwf, 
            pname="MVR-1",
            maxmvr=maxmvr, 
            print_flows=True,
            maxpackages=maxpackages, 
            packages=mvrpack, 
            perioddata=mvrspd,
            budget_filerecord=gwfname + '.mvr.bud',
            filename='{}.mvr'.format(gwfname)
        )
        
        return sim
    return None

# Function to write model files

def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)

# Function to run the model. True is returned if the model runs successfully

def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success = False
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)
    return success

# Function to plot the model results

def plot_results(mf6, idx):
    if config.plotModel:
        print("Plotting model results...")
        sim_name = mf6.name
        fs = USGSFigure(figure_type="graph", verbose=False)
        
        # Generate a plot of FINF distribution
        finf_plt = finf_grad.copy()
        finf_plt[idomain1==0] = np.nan
        
        fig = plt.figure(figsize=figure_size, dpi=300, tight_layout=True)  
        ax = fig.add_subplot(1, 1, 1)          
        plt.imshow(finf_plt, cmap='jet') 
        title = "Precipitation distribution"         
        cbar = plt.colorbar(shrink=0.5)        
        cbar.ax.set_title('Infiltration\nrate\nfactor', pad=20)  
        plt.xlabel("Column Number") 
        plt.ylabel("Row Number") 
        
        letter = chr(ord("@") + idx + 1)
        fs.heading(letter=letter, heading=title)
        
        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}{}".format(sim_name + "-finfFact", config.figure_ext)
            )
            fig.savefig(fpth)
        
        # Start by retrieving some output
        mf6_out_pth = mf6.simulation_data.mfpath.get_sim_path()
        mod_bud_file = list(mf6.model_names)[0] + ".bud"
        sfr_bud_file = list(mf6.model_names)[0] + ".sfr.bud"
        uzf_bud_file = list(mf6.model_names)[0] + ".uzf.bud"
        hed_file = list(mf6.model_names)[0] + ".hds"
        mod_out = os.path.join(mf6_out_pth, mod_bud_file)
        sfr_out = os.path.join(mf6_out_pth, sfr_bud_file)
        uzf_out = os.path.join(mf6_out_pth, uzf_bud_file)
        hds_out = os.path.join(mf6_out_pth, hed_file)
        modobj = bf.CellBudgetFile(mod_out, precision="double")
        sfrobj = bf.CellBudgetFile(sfr_out, precision="double")
        uzfobj = bf.CellBudgetFile(uzf_out, precision="double")
        hdsobj = bf.HeadFile(hds_out)
        
        ckstpkper = modobj.get_kstpkper()
        
        hds = []
        depths = []
        hd_tmp = hdsobj.get_data(kstpkper = ckstpkper[0])
        hd_tmp = np.where(hd_tmp==1e30, 0, hd_tmp)
        hds.append(hd_tmp)
        depths.append(top - hd_tmp[0,:,:])
        depths = np.array(depths)  
        depths = depths[0,:,:]
        
        # Generate a plot of the gw table depths for the steady state stress period
        depths[idomain1==0] = np.nan
        
        fig = plt.figure(figsize=figure_size, dpi=300, tight_layout=True)  
        ax = fig.add_subplot(1, 1, 1)          
        plt.imshow(depths, vmin=0, vmax=25, cmap='jet') 
        cbar = plt.colorbar(shrink=0.5, ticks=[0,5,10,15,20,25])
        cbar.ax.set_yticklabels(['0','5','10','15','20','> 25'])
        cbar.ax.set_title('Depth to\nwater\ntable,\nin m', pad=20)  
        plt.xlabel("Column Number") 
        plt.ylabel("Row Number") 
        title = "Depth To Groundwater"
        
        letter = chr(ord("@") + idx + 2)
        fs.heading(letter=letter, heading=title)
        
        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}{}".format(sim_name + "-gwDepth", config.figure_ext)
            )
            fig.savefig(fpth)
              
        drn_disQ = []
        uzrech   = []
        finf_tot = []
        sfr_gwsw = []
        sfr_flow = []
        rejinf   = []  
        uzet     = []
        gwet     = []
        outflow  = []
        
        for kstpkper in ckstpkper:
            # 1. Compile groundwater discharge to land surface
            drn_tmp = modobj.get_data(kstpkper=kstpkper, text = "      DRN-TO-MVR")
            drn_arr = np.zeros_like(top)
            for itm in drn_tmp[0]:
                i, j = drn_dict_rev[itm[1] - 1]
                drn_arr[i, j] = itm[2]
            drn_disQ.append(drn_arr)
                        
            # 2. Compile groundwater discharge to stream cells
            sfr_tmp = sfrobj.get_data(kstpkper=kstpkper, text = "             GWF")
            sfr_arr = np.zeros_like(top)
            for x, itm in enumerate(sfr_tmp[0]):
                i = sfrcells[x][1]
                j = sfrcells[x][2]
                sfr_arr[i, j] = itm[2]
            sfr_gwsw.append(sfr_arr)
            
            # 3. Compile Infiltrated amount
            uzf_tmp = uzfobj.get_data(kstpkper=kstpkper, text = "    INFILTRATION")
            finf_arr = np.zeros_like(top)
            for itm in uzf_tmp[0]:
                i, j = iuzno_dict_rev[itm[0] - 1]
                finf_arr[i, j] = itm[2]
            finf_tot.append(finf_arr)
            
            # 4. Compile recharge from UZF
            uzrch_tmp = uzfobj.get_data(kstpkper=kstpkper, text = "             GWF")
            uzrch_arr = np.zeros_like(top)
            for itm in uzrch_tmp[0]:
                i, j = iuzno_dict_rev[itm[0] - 1]
                uzrch_arr[i, j] = itm[2]
            uzrech.append(uzrch_arr)
            
            # 5. Compile rejected infiltration
            rejinf_tmp = uzfobj.get_data(kstpkper=kstpkper, text = "  REJ-INF-TO-MVR")
            rejinf_arr = np.zeros_like(top)
            for itm in rejinf_tmp[0]:
                i, j = iuzno_dict_rev[itm[0] - 1]
                rejinf_arr[i, j] = itm[2]
            rejinf.append(rejinf_arr)
            
            # 6. Compile unsat ET
            uzet_tmp = uzfobj.get_data(kstpkper=kstpkper, text = "            UZET")
            uzet_arr = np.zeros_like(top)
            for itm in uzet_tmp[0]:
                i, j = iuzno_dict_rev[itm[0] - 1]
                uzet_arr[i, j] = itm[2]
            uzet.append(uzet_arr) 
            
            # 7. Compile groundwater ET
            gwet_tmp = modobj.get_data(kstpkper=kstpkper, text = "        UZF-GWET")
            gwet_arr = np.zeros_like(top)
            for itm in gwet_tmp[0]:
                i, j = iuzno_dict_rev[itm[1] - 1]
                gwet_arr[i, j] = itm[2]
            gwet.append(gwet_arr)
            
            # 8. Get flows at outlet
            outletQ = sfrobj.get_data(kstpkper=kstpkper, text = "    FLOW-JA-FACE")
            outflow.append(outletQ[0][-1][2])
            
            print('Finished processing stress period ' + str(kstpkper[1] + 1))
        
        drn_disQ = np.array(drn_disQ)
        sfr_gwsw = np.array(sfr_gwsw)
        finf_tot = np.array(finf_tot)
        uzrech = np.array(uzrech)
        rejinf = np.array(rejinf)
        uzet = np.array(uzet)
        gwet = np.array(gwet)  
        outflow = np.array(outflow)
        
        drn_disQ_ts = drn_disQ.sum(axis=-1).sum(axis=-1)
        sfr_gwsw_ts = sfr_gwsw.sum(axis=-1).sum(axis=-1)
        finf_tot_ts = finf_tot.sum(axis=-1).sum(axis=-1)
        uzrech_ts = uzrech.sum(axis=-1).sum(axis=-1)
        rejinf_ts = rejinf.sum(axis=-1).sum(axis=-1)
        uzet_ts = uzet.sum(axis=-1).sum(axis=-1)
        gwet_ts = gwet.sum(axis=-1).sum(axis=-1)
        
        rng = pd.date_range(start='11/1/1999', end='10/31/2000', freq='D')
        data = {'Recharge': abs(uzrech_ts[0:366]),
        	      'Infiltration': abs(finf_tot_ts[0:366]),
        	      'GroundwaterET': abs(gwet_ts[0:366]),
        	      'UnsaturatedZoneET': abs(uzet_ts[0:366])}
        vals1 = pd.DataFrame(data, index=rng)
        
        fig = plt.figure(figsize=figure_size_ts, dpi=300, tight_layout=True)
        ax = fig.add_subplot(1, 1, 1)    
        vals1.Infiltration.plot(style='k-', linewidth=0.3)
        vals1.Recharge.plot(style='b-', linewidth=0.7)
        vals1.GroundwaterET.plot(style='-', color='darkgreen', linewidth=0.7, label='Groundwater ET')
        vals1.UnsaturatedZoneET.plot(style='-', color='orange', linewidth=1, label='Unsaturated Zone ET')
        ax.set_ylim(0, 700000)
        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=True) # labels along the bottom edge are off    
        plt.xlabel('Month', fontsize=8)
        plt.ylabel('Volumetric Rate, $m^3$ per day', fontsize=8)
        plt.legend(prop={'size': 6})
        
        title = "Unsaturated Zone Flow Budget"
        letter = chr(ord("@") + idx + 3)
        fs.heading(letter=letter, heading=title)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}{}".format(sim_name + "-uzFlow", config.figure_ext)
            )
            fig.savefig(fpth)
        
        data_sfr = {'GroundwaterDischarge': abs(drn_disQ_ts[0:366]),
        	          'RejectedInfiltration': abs(rejinf_ts[0:366]),
        	          'gwswExchange': abs(sfr_gwsw_ts[0:366]),
        	          'outflow': outflow[0:366]}
        vals2 = pd.DataFrame(data_sfr, index=rng)
        
        fig = plt.figure(figsize=figure_size_ts, dpi=300, tight_layout=True)
        ax = fig.add_subplot(1, 1, 1)    
        vals2.outflow.plot(style='-', linewidth=0.5, color='royalblue', label='Streamflow at Outlet')
        vals2.GroundwaterDischarge.plot(style='k-', linewidth=0.7, label='Groundwater Discharge')
        vals2.RejectedInfiltration.plot(style='-', linewidth=0.7, color='brown', label='Rejected Infiltration')
        vals2.gwswExchange.plot(style='-', color='darkgreen', linewidth=0.7, label='Groundwater Discharge to Streams')
        ax.set_ylim(0, 350000)
        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=True) # labels along the bottom edge are off    
        plt.xlabel('Month', fontsize=8)
        plt.ylabel('Volumetric Rate, $m^3$ per day', fontsize=8)
        plt.legend(prop={'size': 6})
        
        title = "Surface Water Flow"
        letter = chr(ord("@") + idx + 4)
        fs.heading(letter=letter, heading=title)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}{}".format(sim_name + "-swFlow", config.figure_ext)
            )
            fig.savefig(fpth)


# Function that wraps all of the steps for each scenario
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def scenario(idx, silent=True):
    sim = build_model(example_name)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)

    if success:
        plot_results(sim, idx)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Mehl and Hill (2013) results
    #
    # Two-dimensional transport in a uniform flow field

    scenario(0)
