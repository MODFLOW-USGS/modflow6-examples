# ## FHB example
#
# This example shows how the time series capability in MODFLOW 6 can be
# combined with the constant head and well packages to replicate the
# functionality of the Flow and Head Boundary (FHB) Package in previous
# versions of MODFLOW.
#

# ### FHB Problem Setup
#
# Imports

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import flopy

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set default figure properties

figure_size = (4, 4)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-fhb"

# Model units

length_units = "meters"
time_units = "days"

# Table FHB Model Parameters

nper = 3  # Number of periods
nlay = 1  # Number of layers
ncol = 10  # Number of columns
nrow = 3  # Number of rows
delr = 1000.0  # Column width ($m$)
delc = 1000.0  # Row width ($m$)
top = 50.0  # Top of the model ($m$)
botm_str = "-200.0"  # Layer bottom elevations ($m$)
strt = 0.0  # Starting head ($m$)
icelltype_str = "0"  # Cell conversion type
k11_str = "20.0"  # Horizontal hydraulic conductivity ($m/d$)
ss = 0.01  # Specific storage ($/m$)

# Static temporal data used by TDIS file
# Simulation has 1 steady stress period (1 day)
# and 3 transient stress periods (10 days each).
# Each transient stress period has 120 2-hour time steps.

perlen = [400.0, 200.0, 400.0]
nstp = [10, 4, 6]
tsmult = [1.0, 1.0, 1.0]
tdis_ds = list(zip(perlen, nstp, tsmult))

# parse parameter strings into tuples

botm = [float(value) for value in botm_str.split(",")]
k11 = [float(value) for value in k11_str.split(",")]
icelltype = [int(value) for value in icelltype_str.split(",")]

# Solver parameters

nouter = 50
ninner = 100
hclose = 1e-9
rclose = 1e-6

# ### Functions to build, write, run, and plot the MODFLOW 6 FHB model
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
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord="{} strict".format(rclose),
        )
        gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True)
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
        )
        flopy.mf6.ModflowGwfnpf(
            gwf,
            icelltype=icelltype,
            k=k11,
            save_specific_discharge=True,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)
        flopy.mf6.ModflowGwfsto(
            gwf,
            storagecoefficient=True,
            iconvert=0,
            ss=1.0e-6,
            sy=None,
            transient={0: True},
        )

        chd_spd = []
        chd_spd += [[0, i, 9, "CHDHEAD"] for i in range(3)]
        chd_spd = {0: chd_spd}
        tsdata = [(0.0, 0.0), (307.0, 1.0), (791.0, 5.0), (1000.0, 2.0)]
        tsdict = {
            "timeseries": tsdata,
            "time_series_namerecord": "CHDHEAD",
            "interpolation_methodrecord": "LINEAREND",
        }
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chd_spd,
            timeseries=tsdict,
            pname="CHD",
        )

        wel_spd = []
        wel_spd += [[0, 1, 0, "FLOWRATE"]]
        wel_spd = {0: wel_spd}
        tsdata = [
            (0.0, 2000.0),
            (307.0, 6000.0),
            (791.0, 5000.0),
            (1000.0, 9000.0),
        ]
        tsdict = {
            "timeseries": tsdata,
            "time_series_namerecord": "FLOWRATE",
            "interpolation_methodrecord": "LINEAREND",
        }
        flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data=wel_spd,
            timeseries=tsdict,
            pname="WEL",
        )

        head_filerecord = "{}.hds".format(sim_name)
        budget_filerecord = "{}.cbc".format(sim_name)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        obsdict = {}
        obslist = [
            ["h1_2_1", "head", (0, 1, 0)],
            ["h1_2_10", "head", (0, 1, 9)],
        ]
        obsdict["{}.obs.head.csv".format(sim_name)] = obslist
        obslist = [["icf1", "flow-ja-face", (0, 1, 1), (0, 1, 0)]]
        obsdict["{}.obs.flow.csv".format(sim_name)] = obslist
        obs = flopy.mf6.ModflowUtlobs(
            gwf, print_input=False, continuous=obsdict
        )

        return sim
    return None


# Function to write MODFLOW 6 FHB model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the FHB model.
# True is returned if the model runs successfully
#


@config.timeit
def run_model(sim, silent=False):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent, report=True)
        if not success:
            print(buff)
    return success


# Function to plot the FHB model results.
#
def plot_grid(sim):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model(sim_name)

    fig = plt.figure(figsize=(4, 3.0))
    fig.tight_layout()

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
    pmv.plot_grid()
    pmv.plot_bc(name="CHD")
    pmv.plot_bc(name="WEL")
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-grid{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


def plot_ts(sim):
    fs = USGSFigure(figure_type="graph", verbose=False)
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model(sim_name)
    obsnames = gwf.obs.output.obs_names
    obs_list = [
        gwf.obs.output.obs(f=obsnames[0]),
        gwf.obs.output.obs(f=obsnames[1]),
    ]
    ylabel = ["head (m)", "flow ($m^3/d$)"]
    obs_fig = ("obs-head", "obs-flow", "ghb-obs")
    for iplot, obstype in enumerate(obs_list):
        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_subplot()
        tsdata = obstype.data
        for name in tsdata.dtype.names[1:]:
            ax.plot(tsdata["totim"], tsdata[name], label=name, marker="o")
        ax.set_xlabel("time (d)")
        ax.set_ylabel(ylabel[iplot])
        fs.graph_legend(ax)
        if config.plotSave:
            fpth = os.path.join(
                "..",
                "figures",
                "{}-{}{}".format(
                    sim_name, obs_fig[iplot], config.figure_ext
                ),
            )
            fig.savefig(fpth)
    return


def plot_results(sim, silent=True):
    if config.plotModel:
        plot_grid(sim)
        plot_ts(sim)
    return


# Function that wraps all of the steps for the FHB model
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
    if success:
        plot_results(sim, silent=silent)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(silent=False)


# nosetest end

if __name__ == "__main__":
    # ### FHB Simulation
    #
    # Model grid and simulation results

    simulation()
