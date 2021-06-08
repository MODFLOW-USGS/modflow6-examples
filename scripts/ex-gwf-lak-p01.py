# ## Lake package (LAK) Package problem 1
#
# This is the lake package example problem (test 1) from the
# Lake Package documentation (Merritt and Konikow, 2000).
#

# ### LAK Package Problem 1 Setup
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

figure_size = (6.3, 5.6)
masked_values = (0, 1e30, -1e30)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-lak-p01"

# Model units

length_units = "feet"
time_units = "days"

# Table LAK Package Problem 1 Model Parameters

nper = 1  # Number of periods
nlay = 5  # Number of layers
nrow = 17  # Number of rows
ncol = 17  # Number of columns
top = 500.0  # Top of the model ($ft$)
botm_str = "107., 97., 87., 77., 67."  # Bottom elevations ($ft$)
strt = 115.0  # Starting head ($ft$)
k11 = 30.0  # Horizontal hydraulic conductivity ($ft/d$)
k33_str = (
    "1179., 30., 30., 30., 30."  # Vertical hydraulic conductivity ($ft/d$)
)
ss = 3e-4  # Specific storage ($1/d$)
sy = 0.2  # Specific yield (unitless)
H1 = 160.0  # Constant head on left side of model ($ft$)
H2 = 140.0  # Constant head on right side of model ($ft$)
recharge = 0.0116  # Aereal recharge rate ($ft/d$)
etvrate = 0.0141  # Maximum evapotranspiration rate ($ft/d$)
etvdepth = 15.0  # Evapotranspiration extinction depth ($ft$)
lak_strt = 110.0  # Starting lake stage ($ft$)
lak_etrate = 0.0103  # Lake evaporation rate ($ft/d$)

# parse parameter strings into tuples

botm = [float(value) for value in botm_str.split(",")]
k33 = [float(value) for value in k33_str.split(",")]

# Static temporal data used by TDIS file

tdis_ds = ((5000.0, 100, 1.02),)

# define delr and delc
delr = np.array(
    [
        250.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        500.00,
        500.00,
        500.00,
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
        500.00,
        500.00,
        500.00,
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

# ### Create LAK Package Problem 1 Model Boundary Conditions
#
# Constant head boundary conditions
#

chd_spd = []
for k in range(nlay):
    chd_spd += [[k, i, 0, H1] for i in range(nrow)]
    chd_spd += [[k, i, ncol - 1, H2] for i in range(nrow)]

# LAK Package

lak_conn = [
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

lak_packagedata = [[0, lak_strt, len(lak_conn)]]

lak_spd = [
    [0, "rainfall", recharge],
    [0, "evaporation", lak_etrate],
]

# Solver parameters

nouter = 500
ninner = 100
hclose = 1e-9
rclose = 1e-6


# ### Functions to build, write, run, and plot the MODFLOW 6 LAK Package Problem 1 model
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
            sim, modelname=sim_name, newtonoptions="NEWTON", save_flows=True
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
            ["B", "head", (0, 13, 13)],
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
            print_stage=True,
            nlakes=1,
            noutlets=0,
            packagedata=lak_packagedata,
            connectiondata=lak_conn,
            perioddata=lak_spd,
        )
        obs_file = "{}.lak.obs".format(sim_name)
        csv_file = obs_file + ".csv"
        obs_dict = {
            csv_file: [
                ("stage", "stage", (0,)),
            ]
        }
        lak.obs.initialize(
            filename=obs_file, digits=10, print_input=True, continuous=obs_dict
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


# Function to write MODFLOW 6 LAK Package Problem 1 model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the LAK Package Problem 1 model.
# True is returned if the model runs successfully
#


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
    head[head == 1e30] = lak_results["STAGE"][-1]

    # observation locations
    xcenters, ycenters = gwf.modelgrid.xycenters[0], gwf.modelgrid.xycenters[1]
    p1 = (xcenters[3], ycenters[3])
    p2 = (xcenters[13], ycenters[13])

    fs = USGSFigure(figure_type="map", verbose=False)
    fig = plt.figure(
        figsize=(4, 6.9),
        tight_layout=True,
    )
    plt.axis("off")

    nrows, ncols = 10, 1
    axes = [fig.add_subplot(nrows, ncols, (1, 5))]
    axes.append(fig.add_subplot(nrows, ncols, (6, 8), sharex=axes[0]))

    for idx, ax in enumerate(axes):
        ax.set_xlim(extents[:2])
        if idx == 0:
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
    mm.plot_inactive(color_noflow="#5DBB63")
    mm.plot_grid(lw=0.5, color="black")
    cv = mm.contour_array(
        head,
        levels=np.arange(140, 160, 2),
        linewidths=0.75,
        linestyles="-",
        colors="blue",
        masked_values=masked_values,
    )
    plt.clabel(cv, fmt="%1.0f")
    mm.plot_specific_discharge(spdis, normalize=True, color="0.75")
    ax.plot(p1[0], p1[1], marker="o", mfc="red", mec="black", ms=4)
    ax.plot(p2[0], p2[1], marker="o", mfc="red", mec="black", ms=4)
    ax.set_xlabel("x-coordinate, in feet")
    ax.set_ylabel("y-coordinate, in feet")
    fs.heading(ax, heading="Map view", idx=0)
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
    fs.remove_edge_ticks(ax)

    ax = axes[1]
    xs = flopy.plot.PlotCrossSection(gwf, ax=ax, line={"row": 8})
    xs.plot_array(np.ones(shape3d), head=head, cmap="jet")
    xs.plot_bc("CHD", color="cyan", head=head)
    xs.plot_ibound(color_noflow="#5DBB63", head=head)
    xs.plot_grid(lw=0.5, color="black")
    ax.set_xlabel("x-coordinate, in feet")
    ax.set_ylim(67, 160)
    ax.set_ylabel("Elevation, in feet")
    fs.heading(ax, heading="Cross-section view", idx=1)
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
        lw=0,
        marker="s",
        ms=10,
        mfc="blue",
        mec="black",
        markeredgewidth=0.5,
        label="Water table",
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
        ("STAGE", float),
        ("A", float),
        ("B", float),
    ]

    results = np.zeros((lak_results.shape[0] + 1), dtype=dtype)
    results["time"][1:] = lak_results["time"]
    results["STAGE"][0] = 110.0
    results["STAGE"][1:] = lak_results["STAGE"]
    results["A"][0] = 115.0
    results["A"][1:] = gwf_results["A"]
    results["B"][0] = 115.0
    results["B"][1:] = gwf_results["B"]

    # create the figure
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        sharex=True,
        figsize=(6.3, 3.15),
        constrained_layout=True,
    )

    ax.set_xlim(0, 3000)
    ax.set_ylim(110, 160)
    ax.plot(
        results["time"],
        results["STAGE"],
        lw=0.75,
        ls="--",
        color="black",
        label="Lake stage",
    )
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
    ax.set_xlabel("Simulation time, in days")
    ax.set_ylabel("Head or stage, in feet")
    fs.graph_legend(ax, loc="lower right")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-01{}".format(sim_name, config.figure_ext),
        )
        fig.savefig(fpth)

    return


# Function to plot the LAK Package Problem 1 model results.


def plot_results(sim, silent=True):
    if config.plotModel:
        gwf = sim.get_model(sim_name)

        plot_grid(gwf, silent=silent)

        plot_lak_results(gwf, silent=silent)


# Function that wraps all of the steps for the LAK Package Problem 1 model
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
    # ### LAK Package Problem 1 Simulation
    #

    simulation()
