# ## MODFLOW-NWT Problem 3 example
#
#

# ### MODFLOW-NWT Problem 3 Setup
#
# Imports

import os
import sys

import flopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set figure properties

figure_size = (6.3, 5.6)
masked_values = (1e30, -1e30)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-nwt-p03"

# Model units

length_units = "meters"
time_units = "days"

# Scenario parameters

parameters = {
    "ex-gwf-nwt-p03a": {
        "recharge": "high",
    },
    "ex-gwf-nwt-p03b": {
        "recharge": "low",
    },
}

# Table

nper = 1  # Number of periods
nlay = 1  # Number of layers
nrow = 80  # Number of rows
ncol = 80  # Number of columns
delr = 100.0  # Cell size in the x-direction ($m$)
delc = 100.0  # Cell size in y-direction ($m$)
top = 200.0  # Top of the model ($m$)
k11 = 1.0  # Horizontal hydraulic conductivity ($m/day$)
H1 = 24.0  # Constant head water level ($m$)

# plotting ranges and contour levels

vmin, vmax = 20, 60
smin, smax = 0, 25
bmin, bmax = 0, 90
vlevels = np.arange(vmin, vmax + 5, 5)
slevels = np.arange(smin, smax + 5, 5)
blevels = np.arange(bmin + 10, bmax, 10)
vcolor = "black"
scolor = "black"
bcolor = "black"


# Static temporal data used by TDIS file

tdis_ds = ((365.0, 1, 1.0),)

# Calculate extents, and shape3d

extents = (0, delr * ncol, 0, delc * nrow)
shape3d = (nlay, nrow, ncol)
ticklabels = np.arange(0, 10000, 2000)

# Load the bottom

fpth = os.path.join("..", "data", sim_name, "bottom.txt")
botm = np.loadtxt(fpth).reshape(shape3d)

# Set the starting heads

strt = botm + 20.0

# Load the high recharge rate

fpth = os.path.join("..", "data", sim_name, "recharge_high.txt")
rch_high = np.loadtxt(fpth)

# Generate the low recharge rate from the high recharge rate

rch_low = rch_high.copy() * 1e-3

# Constant head boundary conditions
chd_spd = [
    [0, i, ncol - 1, H1]
    for i in (
        45,
        46,
        47,
    )
]

# Solver parameters

nouter = 500
ninner = 500
hclose = 1e-9
rclose = 1e-6

# ### Functions to build, write, run, and plot the model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(
    name,
    recharge="high",
):
    if config.buildModel:
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6"
        )
        flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
        flopy.mf6.ModflowIms(
            sim,
            print_option="all",
            complexity="simple",
            linear_acceleration="bicgstab",
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord=rclose,
        )
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=sim_name,
            newtonoptions="newton under_relaxation",
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
            icelltype=1,
            k=k11,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

        if recharge == "high":
            rch = rch_high
        elif recharge == "low":
            rch = rch_low
        flopy.mf6.ModflowGwfrcha(gwf, recharge=rch)

        head_filerecord = f"{sim_name}.hds"
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            saverecord=[("HEAD", "ALL")],
        )
        return sim
    return None


# Function to write MODFLOW-NWT Problem 3 model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the model. True is returned if the model runs successfully.
#


@config.timeit
def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)

    return success


# Function to create a figure


def create_figure(nsubs=1, size=(4, 4)):
    fig = plt.figure(figsize=size, constrained_layout=False)
    gs = mpl.gridspec.GridSpec(ncols=10, nrows=7, figure=fig, wspace=5)
    plt.axis("off")

    axes = []
    if nsubs == 1:
        axes.append(fig.add_subplot(gs[:5, :]))
    elif nsubs == 2:
        axes.append(fig.add_subplot(gs[:6, :5]))
        axes.append(fig.add_subplot(gs[:6, 5:], sharey=axes[0]))

    for ax in axes:
        ax.set_xlim(extents[:2])
        ax.set_ylim(extents[2:])
        ax.set_aspect("equal")
        ax.set_xticks(ticklabels)
        ax.set_yticks(ticklabels)

    # legend axis
    axes.append(fig.add_subplot(gs[5:, :]))

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

    return fig, axes


# Function to plot the grid


def plot_grid(gwf, silent=True):
    verbose = not silent
    fs = USGSFigure(figure_type="map", verbose=verbose)

    bot = gwf.dis.botm.array

    fig, axes = create_figure(size=(3.15, 4))
    ax = axes[0]
    mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
    bot_coll = mm.plot_array(bot, vmin=bmin, vmax=bmax)
    mm.plot_bc("CHD", color="cyan")
    cv = mm.contour_array(
        bot,
        levels=blevels,
        linewidths=0.5,
        linestyles="-",
        colors=bcolor,
    )
    plt.clabel(cv, fmt="%1.0f")
    ax.set_xlabel("x-coordinate, in meters")
    ax.set_ylabel("y-coordinate, in meters")
    fs.remove_edge_ticks(ax)

    # legend
    ax = axes[1]
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
        lw=0.5,
        ls="-",
        color=bcolor,
        label="Bottom elevation contour, m",
    )
    fs.graph_legend(ax, loc="center", ncol=2)

    cax = plt.axes([0.275, 0.125, 0.45, 0.025])
    cbar = plt.colorbar(
        bot_coll,
        shrink=0.8,
        orientation="horizontal",
        cax=cax,
    )
    cbar.ax.tick_params(size=0)
    cbar.ax.set_xlabel(r"Bottom Elevation, $m$")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            f"{sim_name}-grid{config.figure_ext}",
        )
        fig.savefig(fpth)

    return


def plot_recharge(gwf, silent=True):
    verbose = not silent
    fs = USGSFigure(figure_type="map", verbose=verbose)

    fig, axes = create_figure(nsubs=2, size=figure_size)
    ax = axes[0]
    mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
    rch_coll = mm.plot_array(rch_high)
    mm.plot_bc("CHD", color="cyan")
    cv = mm.contour_array(
        rch_high,
        levels=[1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6],
        linewidths=0.5,
        linestyles="-",
        colors="black",
    )
    plt.clabel(cv, fmt="%1.0e")
    cbar = plt.colorbar(
        rch_coll,
        shrink=0.8,
        orientation="horizontal",
        ax=ax,
        format="%.0e",
    )
    cbar.ax.tick_params(size=0)
    cbar.ax.set_xlabel(r"Recharge rate, $m/day$")
    ax.set_xlabel("x-coordinate, in meters")
    ax.set_ylabel("y-coordinate, in meters")
    fs.heading(ax, letter="A")
    fs.remove_edge_ticks(ax)

    ax = axes[1]
    mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
    rch_coll = mm.plot_array(rch_low)
    mm.plot_bc("CHD", color="cyan")
    cv = mm.contour_array(
        rch_low,
        levels=[1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9, 7e-9],
        linewidths=0.5,
        linestyles="-",
        colors="black",
    )
    plt.clabel(cv, fmt="%1.0e")
    cbar = plt.colorbar(
        rch_coll,
        shrink=0.8,
        orientation="horizontal",
        ax=ax,
        format="%.0e",
    )
    cbar.ax.tick_params(size=0)
    cbar.ax.set_xlabel(r"Recharge rate, $m/day$")
    ax.set_xlabel("x-coordinate, in meters")
    fs.heading(ax, letter="B")
    fs.remove_edge_ticks(ax)

    # legend
    ax = axes[-1]
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
        lw=0.5,
        ls="-",
        color=bcolor,
        label=r"Recharge rate contour, $m/day$",
    )
    fs.graph_legend(ax, loc="center", ncol=2)

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            f"{sim_name}-01{config.figure_ext}",
        )
        fig.savefig(fpth)

    return


# Function to plot the model results.


def plot_results(idx, sim, silent=True):
    verbose = not silent
    if config.plotModel:
        fs = USGSFigure(figure_type="map", verbose=verbose)
        name = list(parameters.keys())[idx]
        sim_ws = os.path.join(ws, name)
        gwf = sim.get_model(sim_name)

        bot = gwf.dis.botm.array

        if idx == 0:
            plot_grid(gwf, silent=silent)
            plot_recharge(gwf, silent=silent)

        # create MODFLOW 6 head object
        hobj = gwf.output.head()

        # get times
        times = hobj.get_times()

        # extract heads and specific discharge
        head = hobj.get_data(totim=times[0])
        imask = head <= bot + 0.001
        head[imask] = -1e30
        sat_thick = head - botm
        sat_thick[imask] = -1e30

        # Create figure for simulation
        fig, axes = create_figure(nsubs=2, size=figure_size)

        ax = axes[0]
        mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
        h_coll = mm.plot_array(
            head, vmin=vmin, vmax=vmax, masked_values=masked_values, zorder=10
        )
        cv = mm.contour_array(
            head,
            masked_values=masked_values,
            levels=vlevels,
            linewidths=0.5,
            linestyles="-",
            colors=vcolor,
            zorder=10,
        )
        plt.clabel(cv, fmt="%1.0f", zorder=10)
        mm.plot_bc("CHD", color="cyan", zorder=11)
        cbar = plt.colorbar(
            h_coll,
            shrink=0.8,
            orientation="horizontal",
            ax=ax,
            format="%.0f",
        )
        cbar.ax.tick_params(size=0)
        cbar.ax.set_xlabel(r"Water level, $m$")
        ax.set_xlabel("x-coordinate, in meters")
        ax.set_ylabel("y-coordinate, in meters")
        fs.heading(ax, letter="A")
        fs.remove_edge_ticks(ax)

        ax = axes[1]
        mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
        s_coll = mm.plot_array(
            sat_thick,
            vmin=smin,
            vmax=smax,
            masked_values=masked_values,
            zorder=10,
        )
        cv = mm.contour_array(
            sat_thick,
            masked_values=masked_values,
            levels=slevels,
            linewidths=0.5,
            linestyles=":",
            colors=scolor,
            zorder=10,
        )
        plt.clabel(cv, fmt="%1.0f", zorder=10)
        mm.plot_bc("CHD", color="cyan", zorder=11)
        cbar = plt.colorbar(
            s_coll,
            shrink=0.8,
            orientation="horizontal",
            ax=ax,
            format="%.0f",
        )
        cbar.ax.tick_params(size=0)
        cbar.ax.set_xlabel(r"Saturated thickness, $m$")
        ax.set_xlabel("x-coordinate, in meters")
        # ax.set_ylabel("y-coordinate, in meters")
        fs.heading(ax, letter="B")
        fs.remove_edge_ticks(ax)

        # create legend
        ax = axes[-1]
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
            lw=0.5,
            ls="-",
            color=vcolor,
            label="Head contour, m",
        )
        ax.plot(
            -10000,
            -10000,
            lw=0.5,
            ls=":",
            color=scolor,
            label="Saturated thickness contour, m",
        )
        fs.graph_legend(ax, loc="center", ncol=3)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..",
                "figures",
                f"{sim_name}-{idx + 2:02d}{config.figure_ext}",
            )
            fig.savefig(fpth)


# Function that wraps all of the steps for MODFLOW-NWT Problem 3 model
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


# nosetest end

if __name__ == "__main__":
    # ### MODFLOW-NWT Problem 3 Simulation
    #
    # Simulated heads in the MODFLOW-NWT Problem 3 model with high recharge.

    simulation(0)

    # Simulated heads in the MODFLOW-NWT Problem 3 model with low recharge.

    simulation(1)
