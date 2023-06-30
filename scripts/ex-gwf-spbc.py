# ## SPBC example
#
# Periodic boundary condition problem is based on Laattoe and others (2014).
# A MODFLOW 6 GWF-GWF Exchange is used to connect the left column with the
# right column.
#

# ### SPBC Problem Setup
#
# Imports

import os
import sys

import flopy
import matplotlib.pyplot as plt
import numpy as np

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set default figure properties

figure_size = (6, 4)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-spbc"

# Model units

length_units = "meters"
time_units = "days"

# Table SPBC Model Parameters

nper = 1  # Number of periods
nlay = 190  # Number of layers
ncol = 100  # Number of columns
nrow = 1  # Number of rows
delr = 0.06  # Column width ($m$)
delc = 1.0  # Row width ($m$)
delv = 0.03  # Layer thickness ($m$)
top = 0.0  # Top of the model ($m$)
strt = 0.0  # Starting head ($m$)
icelltype = 0  # Cell conversion type
hydraulic_conductivity = 1.0  # Horizontal hydraulic conductivity ($m/d$)

# Static temporal data used by TDIS file
# Simulation has 1 steady stress period (1 day)
# and 3 transient stress periods (10 days each).
# Each transient stress period has 120 2-hour time steps.

perlen = [1.0]
nstp = [1]
tsmult = [1.0]
tdis_ds = list(zip(perlen, nstp, tsmult))

# assign botm

botm = [top - k * delv for k in range(1, nlay + 1)]

# Solver parameters

nouter = 50
ninner = 100
hclose = 1e-9
rclose = 1e-6

# ### Functions to build, write, run, and plot the MODFLOW 6 SPBC model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model():
    if config.buildModel:
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6"
        )
        flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
        flopy.mf6.ModflowIms(
            sim,
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord=f"{rclose} strict",
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
        ihc, cl1, cl2, hwva = 1, delr / 2.0, delr / 2.0, delc
        angldegx = 90.0
        cdist = delr
        exgdata = [
            [(k, 0, 0), (k, 0, ncol - 1), ihc, cl1, cl2, hwva, angldegx, cdist]
            for k in range(nlay)
        ]
        exg = flopy.mf6.ModflowGwfgwf(
            sim,
            exgtype="GWF6-GWF6",
            nexg=len(exgdata),
            auxiliary=["ANGLDEGX", "CDIST"],
            exgmnamea=sim_name,
            exgmnameb=sim_name,
            exchangedata=exgdata,
        )
        flopy.mf6.ModflowGwfnpf(
            gwf,
            icelltype=icelltype,
            k=hydraulic_conductivity,
            save_specific_discharge=True,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)

        hm = 1.0
        lmbda = ncol * delr
        wv = 2 * np.pi / lmbda
        x = gwf.modelgrid.xcellcenters
        chd_head = hm * np.sin(wv * x)
        chd_spd = []
        for j in range(ncol):
            chd_spd.append([0, 0, j, chd_head[0, j]])
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chd_spd,
            pname="CHD",
        )
        head_filerecord = f"{sim_name}.hds"
        budget_filerecord = f"{sim_name}.cbc"
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )
        return sim
    return None


# Function to write MODFLOW 6 SPBC model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the SPBC model.
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


# Function to plot the SPBC model results.
#
def plot_grid(sim):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model(sim_name)

    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    # create MODFLOW 6 head object
    head = gwf.output.head().get_data()

    # create MODFLOW 6 cell-by-cell budget object
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
        gwf.output.budget().get_data(text="DATA-SPDIS", totim=1.0)[0],
        gwf,
    )

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 0})
    # pxs.plot_grid()
    pxs.plot_bc(name="CHD")
    pxs.plot_array(head, cmap="jet")
    levels = np.arange(-1, 1, 0.1)
    cs = pxs.contour_array(
        head, levels=levels, colors="k", linewidths=1.0, linestyles="-"
    )
    pxs.plot_vector(qx, qy, qz, normalize=False, kstep=5, hstep=5)
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("z position (m)")
    ax.set_ylim(-3, 0)
    fs.remove_edge_ticks(ax)
    plt.clabel(cs, fmt="%3.1f", fontsize=5)

    # save figure
    if config.plotSave:
        fpth = os.path.join("..", "figures", f"{sim_name}-grid{config.figure_ext}")
        fig.savefig(fpth)
    return


def plot_results(sim, silent=True):
    if config.plotModel:
        plot_grid(sim)
    return


# Function that wraps all of the steps for the SPBC model
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
    # ### SPBC Simulation
    #
    # Model grid and simulation results

    simulation()
