# ## Reilly Multi-Aquifer Well Problem,
#
# This is the Multi-Aquifer Well from Reilly and others (1989).
#

# ### Reilly MAW Problem Setup
#
# Imports

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set figure properties specific to the

figure_size = (6.3, 4.3)
masked_values = (0, 1e30, -1e30)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-maw-p03"

# Model units

length_units = "feet"
time_units = "days"

# Scenario parameters

parameters = {
    "ex-gwf-maw-p03a": {
        "simulation": "regional",
    },
    "ex-gwf-maw-p03b": {
        "simulation": "multi-aquifer well",
    },
    "ex-gwf-maw-p03c": {
        "simulation": "high conductivity",
    },
}


# function to calculate the well connection conductance

def calc_cond(area, l1, l2, k1, k2):
    c1 = area * k1 / l1
    c2 = area * k2 / l2
    return c1 * c2 / (c1 + c2)


# Table Reilly MAW Problem Model Parameters

nper = 1  # Number of periods
nlay_r = 21  # Number of layers (regional)
nrow_r = 1  # Number of rows (regional)
ncol_r = 200  # Number of columns (regional)
nlay = 41  # Number of layers (local)
nrow = 16  # Number of rows (local)
ncol = 27  # Number of columns (local)
delr_r = 50.0  # Regional column width ($ft$)
delc_r = 1.0  # Regional row width ($ft$)
top = 10.0  # Top of the model ($ft$)
aq_bottom = -205.  # Model bottom elevation ($ft$)
strt = 10.  # Starting head ($ft$)
k11 = 250.  # Horizontal hydraulic conductivity ($ft/d$)
k33 = 50.  # Vertical hydraulic conductivity ($ft/d$)
recharge = 0.004566  # Areal recharge ($ft/d$)
H1 = 0.  # Regional downgradient constant head ($ft$)
maw_loc = (15, 13)  # Row, column location of well
maw_lay = (1, 12)  # Layers with well screen
maw_radius = 0.1333  # Well radius ($ft$)
maw_bot = -65.  # Bottom of the well ($ft$)
maw_highK = 1e+9  # Hydraulic conductivity for well ($ft/d$)

# set delr and delc for the local model

delr = [
    1.000e+1,
    1.000e+1,
    9.000e+0,
    6.000e+0,
    4.000e+0,
    5.000e+0,
    2.000e+0,
    1.330e+0,
    1.250e+0,
    1.000e+0,
    1.000e+0,
    7.500e-1,
    5.000e-1,
    3.330e-1,
    5.000e-1,
    7.500e-1,
    1.000e+0,
    1.000e+0,
    1.250e+0,
    1.333e+0,
    2.000e+0,
    3.000e+0,
    4.000e+0,
    6.000e+0,
    9.000e+0,
    1.000e+1,
    1.000e+1,
]

delc = [
    10,
    9.38,
    9,
    6,
    4,
    3,
    2,
    1.33,
    1.25,
    1,
    1,
    0.75,
    0.5,
    0.375,
    0.25,
    0.1665,
]

# Static temporal data used by TDIS file

tdis_ds = ((1.0, 1, 1.),)

# Define dimensions

extents = (0.0, np.array(delr).sum(), 0.0, np.array(delc).sum())
shape2d = (nrow, ncol)
shape3d = (nlay, nrow, ncol)

# ### Create Reilly MAW Problem Model Boundary Conditions

# MAW Package

nconn = 2 + 3 * (maw_lay[1] - maw_lay[0] + 1)
maw_packagedata = [[0, maw_radius, maw_bot, strt, "SPECIFIED", nconn]]

# Build the MAW connection data

maw_conn = []
i, j = maw_loc
iconn = 0
for k in range(maw_lay[0], maw_lay[1] + 1, 1):
    # connection to layer below
    if k == maw_lay[0]:
        area = delc[i] * delr[j]
        l1 = 2.5
        l2 = 2.5
        cond = calc_cond(area, l1, l2, k33, maw_highK)
        maw_conn.append([0, iconn, k - 1, i, j, top, maw_bot, cond, -999.0])
        iconn += 1

    # connection to layer below
    if k == maw_lay[1]:
        area = delc[i] * delr[j]
        l1 = 2.5
        l2 = 2.5
        cond = calc_cond(area, l1, l2, maw_highK, k33)
        maw_conn.append([0, iconn, k + 1, i, j, top, maw_bot, cond, -999.0])
        iconn += 1

    # connection to left
    area = delc[i] * 5.
    l1 = 0.5 * delr[j]
    l2 = 0.5 * delr[j - 1]
    cond = calc_cond(area, l1, l2, maw_highK, k11)
    maw_conn.append([0, iconn, k, i, j - 1, top, maw_bot, cond, -999.0])
    iconn += 1

    # connection to north
    area = delr[j] * 5.
    l1 = 0.5 * delc[i]
    l2 = 0.5 * delc[i - 1]
    cond = calc_cond(area, l1, l2, maw_highK, k11)
    maw_conn.append([0, iconn, k, i - 1, j, top, maw_bot, cond, -999.0])
    iconn += 1

    # connection to right
    area = delc[i] * 5.
    l1 = 0.5 * delr[j]
    l2 = 0.5 * delr[j + 1]
    cond = calc_cond(area, l1, l2, maw_highK, k11)
    maw_conn.append([0, iconn, k, i, j + 1, top, maw_bot, cond, -999.0])
    iconn += 1

# Solver parameters

nouter = 500
ninner = 100
hclose = 1e-9
rclose = 1e-4


# ### Functions to build, write, run, and plot the MODFLOW 6 Reilly MAW Problem model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(name, simulation="regional"):
    if config.buildModel:
        if simulation == "regional":
            sim = build_regional(name)
        else:
            sim = build_local(name, simulation)

    return sim


def build_regional(name):
    sim_ws = os.path.join(ws, name)
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
    )
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_maximum=nouter,
        outer_dvclose=hclose,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=[rclose, "strict"],
    )
    botm = np.arange(-5, aq_bottom - 10., -10.)
    icelltype = [1] + [0 for k in range(1, nlay_r)]
    gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay_r,
        nrow=nrow_r,
        ncol=ncol_r,
        delr=delr_r,
        delc=delc_r,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=icelltype,
        k=k11,
        k33=k33,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[0, 0, ncol_r - 1, 0.]])
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    head_filerecord = "{}.hds".format(sim_name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        saverecord=[("HEAD", "LAST")],
        printrecord=[("BUDGET", "LAST")],
    )

    return sim


def build_local(name, simulation):
    # get regional heads for constant head boundaries
    pth = list(parameters.keys())[0]
    fpth = os.path.join(ws, pth, "{}.hds".format(sim_name))
    hobj = flopy.utils.HeadFile(fpth)
    h = hobj.get_data()

    # calculate factor for constant heads
    f1 = 0.5 * (delr_r + delr[0]) / delc_r
    f2 = 0.5 * (delr_r + delr[-1]) / delc_r

    # build chd for model
    regional_range = [k for k in range(nlay_r)]
    local_range = [[0]] + [[k, k + 1] for k in range(1, nlay, 2)]
    chd_spd = []
    for kr, kl in zip(regional_range, local_range):
        h1 = h[kr, 0, 0]
        h2 = h[kr, 0, 1]
        hi1 = h1 + f1 * (h2 - h1)
        h1 = h[kr, 0, 2]
        h2 = h[kr, 0, 3]
        hi2 = h1 + f2 * (h2 - h1)
        for k in kl:
            for il in range(nrow):
                chd_spd.append([k, il, 0, hi1])
                chd_spd.append([k, il, ncol - 1, hi2])

    sim_ws = os.path.join(ws, name)
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
    )
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_maximum=nouter,
        outer_dvclose=hclose,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=[rclose, "strict"],
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True)

    botm = np.arange(-5, aq_bottom - 5., -5.)
    icelltype = [1] + [0 for k in range(1, nlay)]
    i, j = maw_loc
    if simulation == "multi-aquifer well":
        k11_sim = k11
        k33_sim = k33
        idomain = np.ones(shape3d, dtype=np.float)
        for k in range(maw_lay[0], maw_lay[1] + 1, 1):
            idomain[k, i, j] = 0
    else:
        k11_sim = np.ones(shape3d, dtype=np.float) * k11
        k33_sim = np.ones(shape3d, dtype=np.float) * k33
        idomain = 1
        for k in range(maw_lay[0], maw_lay[1] + 1, 1):
            k11_sim[k, i, j] = maw_highK
            k33_sim[k, i, j] = maw_highK

    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=icelltype,
        k=k11_sim,
        k33=k33_sim,
        save_specific_discharge=True,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    if simulation == "multi-aquifer well":
        maw = flopy.mf6.ModflowGwfmaw(
            gwf,
            no_well_storage=True,
            nmawwells=1,
            packagedata=maw_packagedata,
            connectiondata=maw_conn,
        )
        obs_file = "{}.maw.obs".format(sim_name)
        csv_file = obs_file + ".csv"
        obs_dict = {
            csv_file: [
                ("head", "head", (0,)),
                ("Q1", "maw", (0,), (0,)),
                ("Q2", "maw", (0,), (1,)),
            ]
        }
        maw.obs.initialize(
            filename=obs_file, digits=10, print_input=True, continuous=obs_dict
        )

    flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("BUDGET", "LAST")],
    )
    return sim


# Function to write MODFLOW 6 Reilly MAW Problem model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the Reilly MAW Problem model.
# True is returned if the model runs successfully
#


def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)
    return success


# Function to plot the lake results


def plot_maw_results(silent=True):
    fs = USGSFigure(figure_type="graph", verbose=False)

    # load the observations
    name = list(parameters.keys())[0]
    fpth = os.path.join(ws, name, "{}.maw.obs.csv".format(sim_name))
    maw0 = np.genfromtxt(fpth, delimiter=",", names=True)
    name = list(parameters.keys())[1]
    fpth = os.path.join(ws, name, "{}.maw.obs.csv".format(sim_name))
    maw1 = np.genfromtxt(fpth, delimiter=",", names=True)

    time = maw0["time"] * 86400.0

    tmin = time[0]
    tmax = time[-1]

    # create the figure
    fig, axes = plt.subplots(
        ncols=1,
        nrows=2,
        sharex=True,
        figsize=figure_size,
        constrained_layout=True,
    )

    ax = axes[0]
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(-1000, 1000)
    ax.semilogx(
        time,
        maw0["Q1"],
        lw=0.75,
        ls="-",
        color="blue",
        label="Upper aquifer",
    )
    ax.semilogx(
        time,
        maw0["Q2"],
        lw=0.75,
        ls="-",
        color="red",
        label="Lower aquifer",
    )
    ax.axhline(0, lw=0.5, color="0.5")
    ax.set_ylabel(" ")
    fs.heading(ax, heading="Non-pumping case", idx=0)
    fs.graph_legend(ax, loc="upper right", ncol=2)

    ax = axes[1]
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(-500, 2500)
    ax.semilogx(
        time,
        maw1["Q1"],
        lw=0.75,
        ls="-",
        color="blue",
        label="Upper aquifer",
    )
    ax.semilogx(
        time,
        maw1["Q2"],
        lw=0.75,
        ls="-",
        color="red",
        label="Lower aquifer",
    )
    ax.axhline(0, lw=0.5, color="0.5")
    ax.set_xlabel(" ")
    ax.set_ylabel(" ")
    for axis in (ax.xaxis,):
        axis.set_major_formatter(mpl.ticker.ScalarFormatter())
    fs.heading(ax, heading="Pumping case", idx=1)

    # add y-axis label that spans both subplots
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # get rid of ticks and spines for legend area
    # ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.patch.set_alpha(0.0)

    ax.set_xlabel("Simulation time, in seconds")
    ax.set_ylabel("Discharge rate, in cubic meters per day")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-01{}".format(sim_name, config.figure_ext),
        )
        fig.savefig(fpth)

    return


# Plot the grid


def plot_grid(silent=True):
    if silent:
        verbosity = 0
    else:
        verbosity = 1
    name = list(parameters.keys())[0]
    sim_ws = os.path.join(ws, name)
    sim = flopy.mf6.MFSimulation.load(
        sim_name=sim_name, sim_ws=sim_ws, verbosity_level=verbosity
    )
    gwf = sim.get_model(sim_name)

    # get regional heads for constant head boundaries
    fpth = os.path.join(ws, name, "{}.hds".format(sim_name))
    hobj = flopy.utils.HeadFile(fpth)
    h = hobj.get_data()

    fs = USGSFigure(figure_type="map", verbose=False)
    fig = plt.figure(
        figsize=(6.3, 3.5),
    )
    plt.axis("off")

    nrows, ncols = 10, 1
    axes = [fig.add_subplot(nrows, ncols, (1, 6))]

    # for idx, ax in enumerate(axes):
    #     ax.set_xlim(extents[:2])
    #     # ax.set_ylim(extents[2:])
    #     # ax.set_aspect("equal")

    # legend axis
    axes.append(fig.add_subplot(nrows, ncols, (7, 10)))

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
    mm = flopy.plot.PlotCrossSection(gwf, ax=ax, line={"row": 0})
    ca = mm.plot_array(h, head=h)
    mm.plot_bc("CHD", color="cyan", head=h)
    mm.plot_grid(lw=0.5, color="0.5")
    cv = mm.contour_array(
        h,
        levels=np.arange(0, 6, 0.5),
        linewidths=0.5,
        linestyles="-",
        colors="black",
        masked_values=masked_values,
    )
    plt.clabel(cv, fmt="%1.1f")
    ax.plot([50, 150, 150, 50, 50], [10, 10, aq_bottom, aq_bottom, 10],
            lw=1.5, color="#39FF14")
    ax.set_xlabel("x-coordinate, in feet")
    ax.set_ylabel("Elevation, in feet")

    # legend
    ax = axes[-1]
    ax.plot(
        -10000,
        -10000,
        lw=0,
        marker="s",
        ms=10,
        mfc="cyan",
        mec="0.5",
        markeredgewidth=0.5,
        label="Constant head",
    )
    ax.plot(
        -10000,
        -10000,
        lw=0,
        marker="s",
        ms=10,
        mfc="none",
        mec="#39FF14",
        markeredgewidth=1.5,
        label="Local model domain",
    )
    ax.plot(
        -10000,
        -10000,
        lw=0.5,
        color="black",
        label="Head contour, $ft$",
    )
    cbar = plt.colorbar(
        ca, shrink=0.5, orientation="horizontal", ax=ax
    )
    cbar.ax.tick_params(size=0)
    cbar.ax.set_xlabel(r"Head, $ft$", fontsize=9)
    fs.graph_legend(ax, loc="lower center", ncol=3)

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-grid{}".format(sim_name, config.figure_ext),
        )
        fig.savefig(fpth)


# Function to plot the Reilly MAW Problem model results.


def plot_results(silent=True):
    if config.plotModel:
        plot_grid(silent=silent)
        # plot_maw_results(silent=silent)
        pass


# Function that wraps all of the steps for the Reilly MAW Problem model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(idx=0, silent=True):
    key = list(parameters.keys())[idx]
    params = parameters[key].copy()

    sim = build_model(key, **params)

    write_model(sim, silent=silent)

    success = run_model(sim, silent=silent)
    assert success, "could not run...{}".format(sim_name)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(idx=0, silent=False)


def test_02():
    simulation(idx=1, silent=False)


def test_03():
    simulation(idx=2, silent=False)


def test_plot():
    plot_results()


# nosetest end

if __name__ == "__main__":
    # ### Reilly MAW Problem Simulation
    #
    # Regional model

    simulation(0)

    # Local model with MAW well

    simulation(1)

    # Local model with high K well

    simulation(2)

    # Plot the results

    plot_results()
