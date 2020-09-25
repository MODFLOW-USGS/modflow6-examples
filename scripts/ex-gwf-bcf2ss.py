# ## Original BCF2SS MODFLOW example
#
# This problem is described in McDonald and Harbaugh (1988) and duplicated in
# Harbaugh and McDonald (1996). This problem is also is distributed with
# MODFLOW-2005 (Harbaugh, 2005) and MODFLOW 6 (Langevin and others, 2017).
#

# ### BCF2SS Problem Setup
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

figure_size = (6, 6)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-bcf2ss"

# Model units

length_units = "feet"
time_units = "days"

# Table BCF2SS Model Parameters

nper = 2  # Number of periods
nlay = 2  # Number of layers
nrow = 10  # Number of rows
ncol = 15  # Number of columns
delr = 500.0  # Column width ($ft$)
delc = 500.0  # Row width ($ft$)
top = 150.0  # Top of the model ($ft$)
botm_str = "50.0, -50."  # Layer bottom elevations ($ft$)
strt = 0.0  # Starting head ($ft$)
wetfct = 1 # WETDRY factor
iwetit = 1  # WETDRY iteration interval
ihdwet = 0  # WETDRY equation (h = BOT + WETFCT (hm - BOT))
icelltype_str = "1, 0"  # Cell conversion type
k11_str = "10.0, 5.0"  # Horizontal hydraulic conductivity ($ft/d$)
k33 = 0.1  # Vertical hydraulic conductivity ($ft/d$)
recharge = 0.004  # Recharge rate ($ft/d$)

# Static temporal data used by TDIS file

tdis_ds = (
    (1., 1., 1),
    (1., 1., 1),
)

# parse parameter strings into tuples

botm = [float(value) for value in botm_str.split(",")]
icelltype = [int(value) for value in icelltype_str.split(",")]
k11 = [float(value) for value in k11_str.split(",")]


# Load the wetdry array for layer 1

pth = os.path.join("..", "data", sim_name, "wetdry01.txt")
wetdry_layer0 = np.loadtxt(pth, )
wetdry = [wetdry_layer0, 0]


# ### Create BCF2SS Model Boundary Conditions

# Well boundary conditions

wel_spd = {
    1: [
        [1, 2, 3, -35000.0],
        [1, 7, 3, -35000.0],
    ]
}

# River boundary conditions

riv_spd = {
    0: [[1, i, 14, 0., 10000., -5] for i in range(nrow)]
}

# Solver parameters

nouter = 500
ninner = 100
hclose = 1e-6
rclose = 1e-3
relax = 0.97

# ### Functions to build, write, run, and plot the MODFLOW 6 TWRI model
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
            rcloserecord=[rclose, "strict"],
            relaxation_factor=relax,
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
            rewet_record=["wetfct", wetfct, "iwetit", iwetit, "ihdwet", ihdwet],
            wetdry=wetdry,
            icelltype=icelltype,
            k=k11,
            k33=k33,
            save_specific_discharge=True,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)
        flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd)
        flopy.mf6.ModflowGwfwel(gwf, stress_period_data=wel_spd)
        flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
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


# Function to write MODFLOW 6 TWRI model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the TWRI model.
# True is returned if the model runs successfully
#


def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)

    return success


# Function to plot the BCF2SS model results.
#


def plot_results(sim, silent=True):
    if config.plotModel:

        fs = USGSFigure(figure_type="map", verbose=False)
        sim_ws = os.path.join(ws, sim_name)
        gwf = sim.get_model(sim_name)

        # create MODFLOW 6 head object
        file_name = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        hobj = flopy.utils.HeadFile(fpth)

        # create MODFLOW 6 cell-by-cell budget object
        file_name = gwf.oc.budget_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        cobj = flopy.utils.CellBudgetFile(fpth, precision="double")

        # extract heads
        head = hobj.get_data(totim=1)

        # plot grid
        fig = plt.figure(figsize=(6.8, 3), constrained_layout=True)
        gs = mpl.gridspec.GridSpec(7, 10, figure=fig, wspace=5)
        plt.axis("off")

        ax = fig.add_subplot(gs[:, 0:7])
        ax.set_aspect('equal')
        mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
        mm.plot_grid(lw=0.5, color="0.5")
        mm.plot_bc(ftype="WEL", kper=1, plotAll=True)
        mm.plot_bc(ftype="RIV", color="green", plotAll=True)
        ax.set_ylabel("y-coordinate, in feet")
        ax.set_xlabel("x-coordinate, in feet")
        fs.heading(ax, letter="A", heading="Map view")
        fs.remove_edge_ticks(ax)

        ax = fig.add_subplot(gs[0:5, 7:])
        mm = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 7})
        mm.plot_grid(lw=0.5, color="0.5")
        mm.plot_array(np.ones((nlay, nrow, ncol)), head=head, cmap="jet")
        mm.plot_bc(ftype="WEL", kper=1)
        mm.plot_bc(ftype="RIV", color="green", head=head)

        # items for legend
        mm.ax.plot(
            -1000,
            -1000,
            "s",
            ms=5,
            color="green",
            mec="black",
            mew=0.5,
            label="River",
        )
        mm.ax.plot(
            -1000,
            -1000,
            "s",
            ms=5,
            color="red",
            mec="black",
            mew=0.5,
            label="Well",
        )
        mm.ax.plot(
            -1000,
            -1000,
            "s",
            ms=5,
            color="blue",
            mec="black",
            mew=0.5,
            label="Steady-state\nwater level",
        )
        fs.graph_legend(
            mm.ax,
            ncol=2,
            bbox_to_anchor=(0.5, -0.75),
            borderaxespad=0,
            frameon=False,
            loc="lower center",
        )

        ax.set_ylabel("Elevation, in feet")
        ax.set_xlabel("x-coordinate along model row 8, in feet")
        fs.heading(ax, letter="B", heading="Cross-section view")
        fs.remove_edge_ticks(ax)

        # figure with heads in each layer

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}-grid{}".format(sim_name, config.figure_ext)
            )
            fig.savefig(fpth)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}{}".format(sim_name, config.figure_ext)
            )
            # fig.savefig(fpth)


# Function that wraps all of the steps for the TWRI model
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
    # ### BCF2SS Simulation
    #
    # Simulated heads in model the unconfined, middle, and lower aquifers (model layers
    # 1, 3, and 5) are shown in the figure below. MODFLOW-2005 results for a quasi-3D
    # model are also shown. The location of drain (green) and well (gray) boundary
    # conditions, normalized specific discharge, and head contours (25 ft contour
    # intervals) are also shown.

    simulation()
