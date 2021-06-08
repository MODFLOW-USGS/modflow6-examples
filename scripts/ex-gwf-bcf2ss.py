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

# Load the wetdry array for layer 1

pth = os.path.join("..", "data", sim_name, "wetdry01.txt")
wetdry_layer0 = np.loadtxt(
    pth,
)


# Scenario parameters

parameters = {
    "ex-gwf-bcf2ss-p01a": {
        "rewet": True,
        "wetfct": 1.0,
        "iwetit": 1,
        "ihdwet": 0,
        "linear_acceleration": "cg",
        "newton": None,
    },
    "ex-gwf-bcf2ss-p02a": {
        "rewet": False,
        "wetfct": None,
        "iwetit": None,
        "ihdwet": None,
        "linear_acceleration": "bicgstab",
        "newton": "NEWTON",
    },
}

# Table BCF2SS Model Parameters

nper = 2  # Number of periods
nlay = 2  # Number of layers
nrow = 10  # Number of rows
ncol = 15  # Number of columns
delr = 500.0  # Column width ($ft$)
delc = 500.0  # Row width ($ft$)
top = 150.0  # Top of the model ($ft$)
botm_str = "50.0, -50."  # Layer bottom elevations ($ft$)
icelltype_str = "1, 0"  # Cell conversion type
k11_str = "10.0, 5.0"  # Horizontal hydraulic conductivity ($ft/d$)
k33 = 0.1  # Vertical hydraulic conductivity ($ft/d$)
strt = 0.0  # Starting head ($ft$)
recharge = 0.004  # Recharge rate ($ft/d$)

# Static temporal data used by TDIS file

tdis_ds = (
    (1.0, 1.0, 1),
    (1.0, 1.0, 1),
)

# parse parameter strings into tuples

botm = [float(value) for value in botm_str.split(",")]
icelltype = [int(value) for value in icelltype_str.split(",")]
k11 = [float(value) for value in k11_str.split(",")]


# ### Create BCF2SS Model Boundary Conditions

# Well boundary conditions

wel_spd = {
    1: [
        [1, 2, 3, -35000.0],
        [1, 7, 3, -35000.0],
    ]
}

# River boundary conditions

riv_spd = {0: [[1, i, 14, 0.0, 10000.0, -5] for i in range(nrow)]}

# Solver parameters

nouter = 500
ninner = 100
hclose = 1e-6
rclose = 1e-3
relax = 0.97

# ### Functions to build, write, run, and plot the MODFLOW 6 TWRI model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(
    name, rewet, wetfct, iwetit, ihdwet, linear_acceleration, newton
):
    if config.buildModel:
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
        )
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_ds, time_units=time_units
        )
        flopy.mf6.ModflowIms(
            sim,
            linear_acceleration=linear_acceleration,
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord="{} strict".format(rclose),
            relaxation_factor=relax,
        )
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=sim_name, save_flows=True, newtonoptions=newton
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
        if rewet:
            rewet_record = [
                "wetfct",
                wetfct,
                "iwetit",
                iwetit,
                "ihdwet",
                ihdwet,
            ]
            wetdry = [wetdry_layer0, 0]
        else:
            rewet_record = None
            wetdry = None

        flopy.mf6.ModflowGwfnpf(
            gwf,
            rewet_record=rewet_record,
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


# Function to plot the BCF2SS model results with heads in each layer.
#


def plot_simulated_results(num, gwf, ho, co, silent=True):
    verbose = not silent
    fs = USGSFigure(figure_type="map", verbose=verbose)

    botm_arr = gwf.dis.botm.array

    fig = plt.figure(figsize=(6.8, 6), constrained_layout=False)
    gs = mpl.gridspec.GridSpec(ncols=10, nrows=7, figure=fig, wspace=5)
    plt.axis("off")

    ax1 = fig.add_subplot(gs[:3, :5])
    ax2 = fig.add_subplot(gs[:3, 5:], sharey=ax1)
    ax3 = fig.add_subplot(gs[3:6, :5], sharex=ax1)
    ax4 = fig.add_subplot(gs[3:6, 5:], sharex=ax1, sharey=ax1)
    ax5 = fig.add_subplot(gs[6, :])
    axes = [ax1, ax2, ax3, ax4, ax5]

    labels = ("A", "B", "C", "D")
    aquifer = ("Upper aquifer", "Lower aquifer")
    cond = ("natural conditions", "pumping conditions")
    vmin, vmax = -10, 140
    masked_values = [1e30, -1e30]
    levels = [
        np.arange(0, 130, 10),
        (10, 20, 30, 40, 50, 55, 60),
    ]
    plot_number = 0
    for idx, totim in enumerate(
        (
            1,
            2,
        )
    ):
        head = ho.get_data(totim=totim)
        head[head < botm_arr] = -1e30
        spdis = co.get_data(text="DATA-SPDIS", kstpkper=(0, totim - 1))

        for k in range(nlay):
            ax = axes[plot_number]
            ax.set_aspect("equal")
            mm = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=k)
            mm.plot_grid(lw=0.5, color="0.5")
            cm = mm.plot_array(
                head, masked_values=masked_values, vmin=vmin, vmax=vmax
            )
            mm.plot_bc(ftype="WEL", kper=totim - 1)
            mm.plot_bc(ftype="RIV", color="green", kper=0)
            mm.plot_specific_discharge(spdis, normalize=True, color="0.75")
            cn = mm.contour_array(
                head,
                masked_values=masked_values,
                levels=levels[idx],
                colors="black",
                linewidths=0.5,
            )
            plt.clabel(cn, fmt="%3.0f")
            heading = "{} under {}".format(aquifer[k], cond[totim - 1])
            fs.heading(ax, letter=labels[plot_number], heading=heading)
            fs.remove_edge_ticks(ax)

            plot_number += 1

    # set axis labels
    ax1.set_ylabel("y-coordinate, in feet")
    ax3.set_ylabel("y-coordinate, in feet")
    ax3.set_xlabel("x-coordinate, in feet")
    ax4.set_xlabel("x-coordinate, in feet")

    # legend
    ax = axes[-1]
    ax.set_ylim(1, 0)
    ax.set_xlim(-5, 5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.patch.set_alpha(0.0)

    # items for legend
    ax.plot(
        -1000,
        -1000,
        "s",
        ms=5,
        color="green",
        mec="black",
        mew=0.5,
        label="River",
    )
    ax.plot(
        -1000,
        -1000,
        "s",
        ms=5,
        color="red",
        mec="black",
        mew=0.5,
        label="Well",
    )
    ax.plot(
        -1000,
        -1000,
        "s",
        ms=5,
        color="none",
        mec="black",
        mew=0.5,
        label="Dry cell",
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
    ax.plot(
        -1000,
        -1000,
        lw=0.5,
        color="black",
        label="Head, in feet",
    )
    fs.graph_legend(
        ax,
        ncol=5,
        frameon=False,
        loc="upper center",
    )

    cbar = plt.colorbar(cm, ax=ax, shrink=0.5, orientation="horizontal")
    cbar.ax.set_xlabel("Head, in feet")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-{:02d}{}".format(sim_name, num, config.figure_ext),
        )
        fig.savefig(fpth)


# Function to plot simulated results for a simulation


def plot_results(silent=True):
    if config.plotModel:
        verbose = not silent
        if silent:
            verbosity_level = 0
        else:
            verbosity_level = 1

        fs = USGSFigure(figure_type="map", verbose=verbose)
        name = list(parameters.keys())[0]
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation.load(
            sim_name=sim_name, sim_ws=sim_ws, verbosity_level=verbosity_level
        )
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
        fig = plt.figure(figsize=(6.8, 3.5), constrained_layout=True)
        gs = mpl.gridspec.GridSpec(nrows=8, ncols=10, figure=fig, wspace=5)
        plt.axis("off")

        ax = fig.add_subplot(gs[:7, 0:7])
        ax.set_aspect("equal")
        mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
        mm.plot_bc(ftype="WEL", kper=1, plotAll=True)
        mm.plot_bc(ftype="RIV", color="green", plotAll=True)
        mm.plot_grid(lw=0.5, color="0.5")
        ax.set_ylabel("y-coordinate, in feet")
        ax.set_xlabel("x-coordinate, in feet")
        fs.heading(ax, letter="A", heading="Map view")
        fs.remove_edge_ticks(ax)

        ax = fig.add_subplot(gs[:5, 7:])
        mm = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 7})
        mm.plot_array(np.ones((nlay, nrow, ncol)), head=head, cmap="jet")
        mm.plot_bc(ftype="WEL", kper=1)
        mm.plot_bc(ftype="RIV", color="green", head=head)
        mm.plot_grid(lw=0.5, color="0.5")
        ax.set_ylabel("Elevation, in feet")
        ax.set_xlabel("x-coordinate along model row 8, in feet")
        fs.heading(ax, letter="B", heading="Cross-section view")
        fs.remove_edge_ticks(ax)

        # items for legend
        ax = fig.add_subplot(gs[7, :])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_color("none")
        ax.spines["bottom"].set_color("none")
        ax.spines["left"].set_color("none")
        ax.spines["right"].set_color("none")
        ax.patch.set_alpha(0.0)
        ax.plot(
            -1000,
            -1000,
            "s",
            ms=5,
            color="green",
            mec="black",
            mew=0.5,
            label="River",
        )
        ax.plot(
            -1000,
            -1000,
            "s",
            ms=5,
            color="red",
            mec="black",
            mew=0.5,
            label="Well",
        )
        ax.plot(
            -1000,
            -1000,
            "s",
            ms=5,
            color="blue",
            mec="black",
            mew=0.5,
            label="Steady-state water table",
        )
        fs.graph_legend(
            ax,
            ncol=3,
            frameon=False,
            loc="upper center",
        )

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..",
                "figures",
                "{}-grid{}".format(sim_name, config.figure_ext),
            )
            fig.savefig(fpth)

        # figure with wetdry array
        fig = plt.figure(figsize=(4.76, 3), constrained_layout=True)
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect("equal")
        mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
        wd = mm.plot_array(wetdry_layer0)
        mm.plot_grid(lw=0.5, color="0.5")
        cbar = plt.colorbar(wd, shrink=0.5)
        cbar.ax.set_ylabel("WETDRY parameter")
        ax.set_ylabel("y-coordinate, in feet")
        ax.set_xlabel("x-coordinate, in feet")
        fs.remove_edge_ticks(ax)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}-01{}".format(sim_name, config.figure_ext)
            )
            fig.savefig(fpth)

        # plot simulated rewetting results
        plot_simulated_results(2, gwf, hobj, cobj)

        # plot simulated newton results
        name = list(parameters.keys())[1]
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation.load(
            sim_name=sim_name, sim_ws=sim_ws, verbosity_level=verbosity_level
        )
        gwf = sim.get_model(sim_name)

        # create MODFLOW 6 head object
        file_name = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        hobj = flopy.utils.HeadFile(fpth)

        # create MODFLOW 6 cell-by-cell budget object
        file_name = gwf.oc.budget_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        cobj = flopy.utils.CellBudgetFile(fpth, precision="double")

        # plot the newton results
        plot_simulated_results(3, gwf, hobj, cobj)


# Function that wraps all of the steps for the TWRI model
#
# 1. build_model,
# 2. write_model, and
# 3. run_model
# 4. plot_results.
#


def simulation(idx, silent=True):
    key = list(parameters.keys())[idx]
    params = parameters[key].copy()

    sim = build_model(key, **params)

    write_model(sim, silent=silent)

    success = run_model(sim, silent=silent)

    assert success, "could not run...{}".format(key)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(0, silent=False)


def test_02():
    simulation(1, silent=False)


def test_plot():
    plot_results(silent=False)


# nosetest end

if __name__ == "__main__":
    # ### BCF2SS Simulation
    #
    #  Node Property Flow Package with rewetting option

    simulation(0)

    # Newton-Raphson formulation

    simulation(1)

    # Simulated water levels and normalized specific discharge vectors in the
    # upper and lower aquifers under natural and pumping conditions using (1) the
    # rewetting option in the Node Property Flow (NPF) Package with the
    # Standard Conductance Formulation and (2) the Newton-Raphson formulation.
    # A. Upper aquifer results under natural conditions. B. Lower aquifer results
    # under natural conditions C. Upper aquifer results under pumping conditions.
    # D. Lower aquifer results under pumping conditions

    plot_results()
