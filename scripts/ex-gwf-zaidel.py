# ## Zaidel (2013) example
#
# This problem is described in Zaidel (2013) and represents a discontinuous
# water table configuration over a stairway impervious base.
#

# ### Zaidel (2013) Problem Setup
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

# Set figure properties specific to the

figure_size = (6.3, 2.5)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-zaidel"

# Model units

length_units = "meters"
time_units = "days"
# Scenario parameters

parameters = {
    "ex-gwf-zaidel-p01a": {
        "H2": 1.0,
    },
    "ex-gwf-zaidel-p02a": {
        "H2": 10.0,
    },
}

# Table

nper = 1  # Number of periods
nlay = 1  # Number of layers
nrow = 1  # Number of rows
ncol = 200  # Number of columns
delr = 5.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
top = 25.0  # Top of the model ($m$)
strt = 23.0  # Starting head ($m$)
icelltype = 1  # Cell conversion type
k11 = 0.0001  # Horizontal hydraulic conductivity ($m/day$)
H1 = 23.0  # Constant head in column 1 ($m$)

# Static temporal data used by TDIS file

tdis_ds = ((1.0, 1, 1.0),)

# Build stairway bottom

botm = np.zeros((nlay, nrow, ncol), dtype=float)
base = 20.0
for j in range(ncol):
    botm[0, :, j] = base
    if j + 1 in (40, 80, 120, 160):
        base -= 5

# Solver parameters

nouter = 500
ninner = 50
hclose = 1e-9
rclose = 1e-6

# ### Functions to build, write, run, and plot the Zaidel model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(H2=1.0):
    if config.buildModel:
        # Constant head cells are specified on the left and right edge of the model
        chd_spd = [
            [0, 0, 0, H1],
            [0, 0, ncol - 1, H2],
        ]

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
            rcloserecord="{} strict".format(rclose),
        )
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=sim_name, newtonoptions="NEWTON"
        )
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
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

        head_filerecord = "{}.hds".format(sim_name)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            saverecord=[("HEAD", "ALL")],
        )
        return sim
    return None


# Function to write Zaidel model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the Zaidel model.
# True is returned if the model runs successfully
#


def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)

    return success


# Function to plot the Zaidel model results.
#


def plot_results(idx, sim, silent=True):
    verbose = not silent
    if config.plotModel:
        fs = USGSFigure(figure_type="map", verbose=verbose)
        sim_ws = os.path.join(ws, sim_name)
        gwf = sim.get_model(sim_name)
        xedge = gwf.modelgrid.xvertices[0]
        zedge = np.array([botm[0, 0, 0]] + botm.flatten().tolist())

        # create MODFLOW 6 head object
        file_name = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        hobj = flopy.utils.HeadFile(fpth)

        # extract heads
        head = hobj.get_data()
        vmin, vmax = 0, 25

        # Create figure for simulation
        extents = (0, ncol * delr, -1, 25.0)
        fig, ax = plt.subplots(
            ncols=1,
            nrows=1,
            figsize=figure_size,
            dpi=300,
            constrained_layout=True,
            sharey=True,
        )

        ax.set_xlim(extents[:2])
        ax.set_ylim(extents[2:])

        fmp = flopy.plot.PlotCrossSection(
            model=gwf, ax=ax, extent=extents, line={"row": 0}
        )
        ax.fill_between(xedge, zedge, y2=-1, color="0.75", step="pre", lw=0.0)
        plot_obj = fmp.plot_array(head, head=head, vmin=vmin, vmax=vmax)
        fmp.plot_bc("CHD", color="cyan", head=head)
        ax.set_xlabel("x-coordinate, in meters")
        ax.set_ylabel("Elevation, in meters")

        # create legend
        ax.plot(
            -10000,
            -10000,
            lw=0,
            marker="s",
            ms=10,
            mfc="cyan",
            mec="cyan",
            label="Constant Head",
        )
        ax.plot(
            -10000,
            -10000,
            lw=0,
            marker="s",
            ms=10,
            mfc="0.75",
            mec="0.75",
            label="Model Base",
        )
        fs.graph_legend(ax, ncol=2, loc="upper right")

        # plot colorbar
        cax = plt.axes([0.62, 0.76, 0.325, 0.025])
        cbar = plt.colorbar(
            plot_obj, shrink=0.8, orientation="horizontal", cax=cax
        )
        cbar.ax.tick_params(size=0)
        cbar.ax.set_xlabel(r"Head, $m$", fontsize=9)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..",
                "figures",
                "{}-{:02d}{}".format(sim_name, idx + 1, config.figure_ext),
            )
            fig.savefig(fpth)


# Function that wraps all of the steps for the TWRI model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(idx, silent=True):
    key = list(parameters.keys())[idx]
    params = parameters[key].copy()

    sim = build_model(**params)

    write_model(sim, silent=silent)

    success = run_model(sim, silent=silent)

    if success:
        plot_results(idx, sim, silent=silent)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(0, silent=False)


def test_02():
    simulation(1, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Zaidel Simulation
    #
    # Simulated heads in the Zaidel model with H2 = 1.

    simulation(0)

    # Simulated heads in the Zaidel model with H2 = 10.

    simulation(1)
