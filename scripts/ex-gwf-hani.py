# ## Hani example
#
# Simple steady state model using a regular MODFLOW grid to simulate the
# response of an anisotropic confined aquifer to a pumping well. A
# constant-head boundary condition surrounds the active domain.  K22 is set
# to 0.01.  Drawdown is more pronounced in the K11 direction.
#

# ### Hani Problem Setup
#
# Imports

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import flopy
import flopy.utils.cvfdutil

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set default figure properties

figure_size = (3.5, 3.5)

# Base simulation and model name and workspace

ws = config.base_ws

# Model units

length_units = "meters"
time_units = "days"

# Scenario parameters

parameters = {
    "ex-gwf-hanir": {"angle1": 0, "xt3d": False},
    "ex-gwf-hanix": {"angle1": 25, "xt3d": True},
    "ex-gwf-hanic": {"angle1": 90, "xt3d": False},
}

# Table Hani Model Parameters

nper = 1  # Number of periods
nlay = 1  # Number of layers
nrow = 51  # Number of rows
ncol = 51  # Number of columns
delr = 10.0  # Spacing along rows ($m$)
delc = 10.0  # Spacing along columns ($m$)
top = 0.0  # Top of the model ($m$)
botm = -10.0  # Layer bottom elevations ($m$)
strt = 0.0  # Starting head ($m$)
icelltype = 0  # Cell conversion type
k11 = 1.0  # Horizontal hydraulic conductivity in the 11 direction ($m/d$)
k22 = 0.01  # Horizontal hydraulic conductivity in the 22 direction ($m/d$)
pumping_rate = -1.0  # Pumping rate ($m^3/d$)

# Static temporal data used by TDIS file
# Simulation has 1 steady stress period (1 day)
# and 3 transient stress periods (10 days each).
# Each transient stress period has 120 2-hour time steps.

perlen = [1.0]
nstp = [1]
tsmult = [1.0, 1.0, 1.0]
tdis_ds = list(zip(perlen, nstp, tsmult))

nouter = 50
ninner = 100
hclose = 1e-9
rclose = 1e-6

# ### Functions to build, write, run, and plot the MODFLOW 6 Hani model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(sim_name, angle1, xt3d):
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
            linear_acceleration="bicgstab",
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord=[rclose, "strict"],
        )
        gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True)
        flopy.mf6.ModflowGwfdis(
            gwf,
            length_units=length_units,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            top=top,
            botm=botm,
        )
        flopy.mf6.ModflowGwfnpf(
            gwf,
            icelltype=icelltype,
            k=k11,
            k22=k22,
            angle1=angle1,
            save_specific_discharge=True,
            xt3doptions=xt3d,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)

        ibd = -1 * np.ones((nrow, ncol), dtype=np.int)
        ibd[1:-1, 1:-1] = 1
        chdrow, chdcol = np.where(ibd == -1)
        chd_spd = [[0, i, j, 0.0] for i, j in zip(chdrow, chdcol)]
        flopy.mf6.ModflowGwfchd(
            gwf, stress_period_data=chd_spd, pname="CHD",
        )
        flopy.mf6.ModflowGwfwel(
            gwf, stress_period_data=[0, 25, 25, pumping_rate], pname="WEL",
        )
        head_filerecord = "{}.hds".format(sim_name)
        budget_filerecord = "{}.cbc".format(sim_name)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        return sim
    return None


# Function to write MODFLOW 6 Hani model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the FHB model.
# True is returned if the model runs successfully
#


def run_model(sim, silent=False):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent, report=True)
        if not success:
            print(buff)
    return success


# Function to plot the Hani model results.
#
def plot_grid(idx, sim):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = list(parameters.keys())[idx]
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model(sim_name)

    fig = plt.figure(figsize=figure_size)
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


def plot_head(idx, sim):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = list(parameters.keys())[idx]
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model(sim_name)

    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    fname = os.path.join(sim_ws, "{}.hds".format(sim_name))
    hdobj = flopy.utils.HeadFile(fname)
    head = hdobj.get_data()

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
    cb = pmv.plot_array(0 - head, cmap="jet", alpha=0.25)
    cs = pmv.contour_array(0 - head, levels=np.arange(0.1, 1, 0.1))
    cbar = plt.colorbar(cb, shrink=0.25)
    cbar.ax.set_xlabel(r"Drawdown, ($m$)")
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-head{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


def plot_results(idx, sim, silent=True):
    if config.plotModel:
        if idx == 0:
            plot_grid(idx, sim)
        plot_head(idx, sim)
    return


# Function that wraps all of the steps for the FHB model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(idx, silent=True):
    key = list(parameters.keys())[idx]
    params = parameters[key].copy()
    sim = build_model(key, **params)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(idx, sim, silent=silent)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(0, silent=False)


def test_02():
    simulation(1, silent=False)


def test_03():
    simulation(2, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Hani Simulation
    #
    # Simulated heads in the Hani model with anisotropy in x direction.

    simulation(0)

    # Simulated heads in the Hani model with anisotropy in y direction.

    simulation(1)

    # Simulated heads in the Hani model with anisotropy rotated 15 degrees.

    simulation(2)
