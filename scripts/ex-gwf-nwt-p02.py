# ## Flow diversion example
#
#

# ### Flow diversion Problem Setup
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

# Set figure properties

figure_size = (4, 5.33)
masked_values = (1e30, -1e30)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-nwt-p02"

# Model units

length_units = "feet"
time_units = "days"

# Scenario parameters

parameters = {
    "ex-gwf-nwt-p02a": {
        "newton": True,
    },
    "ex-gwf-bump-p02b": {
        "rewet": True,
        "wetfct": 1.0,
        "iwetit": 1,
        "ihdwet": 0,
        "wetdry": 2.0,
    },
}

# Table

nper = 4  # Number of periods
nlay = 14  # Number of layers
nrow = 40  # Number of rows
ncol = 40  # Number of columns
delr = 125.0  # Column width ($ft$)
delc = 125.0  # Row width ($ft$)
top = 80.0  # Top of the model ($ft$)
k11 = 5.0  # Horizontal hydraulic conductivity ($ft/day$)
k33 = 0.25  # Horizontal hydraulic conductivity ($ft/day$)
ss = 0.0002  # Specific storage ($1/day$)
sy = 0.2  # Specific yield (unitless)
H1 = 25  # Constant head along left and lower edges and starting head ($ft$)
rech = 0.05 # Recharge rate ($ft/day$)


# Static temporal data used by TDIS file

tdis_ds = (
    (190.0, 10, 1.0),
    (518.0, 2, 1.0),
    (1921.0, 17, 1.0),
    (1.0, 1, 1.0),
)

# Calculate extents, and shape3d

extents = (0, delr * ncol, 0, delc * nrow)
shape3d = (nlay, nrow, ncol)

# Create the bottom

botm = np.arange(65., -5., 5.)

# Create icelltype (which is the same as iconvert)

icelltype

fpth = os.path.join("..", "data", sim_name, "bottom.txt")
botm = np.loadtxt(fpth).reshape(shape3d)

# create a cylinder

cylinder = botm.copy()
cylinder[cylinder < 7.5] = 0.0
cylinder[cylinder >= 7.5] = 20.0

# Constant head boundary conditions
chd_spd = [[0, i, 0, H1] for i in range(nrow)]
chd_spd += [[0, i, ncol - 1, H2] for i in range(nrow)]

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
    newton=False,
    rewet=False,
    cylindrical=False,
    wetfct=None,
    iwetit=None,
    ihdwet=None,
    wetdry=None,
):
    if config.buildModel:
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
        )
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_ds, time_units=time_units
        )
        if newton:
            linear_acceleration = "bicgstab"
            newtonoptions = "under_relaxation"
        else:
            linear_acceleration = "cg"
            newtonoptions = None

        flopy.mf6.ModflowIms(
            sim,
            print_option="ALL",
            linear_acceleration=linear_acceleration,
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord=rclose,
        )
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=sim_name,
            newtonoptions=newtonoptions,
            save_flows=True,
        )
        if cylindrical:
            bot = cylinder
        else:
            bot = botm
        flopy.mf6.ModflowGwfdis(
            gwf,
            length_units=length_units,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=bot,
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
        else:
            rewet_record = None
        flopy.mf6.ModflowGwfnpf(
            gwf,
            rewet_record=rewet_record,
            icelltype=1,
            k=k11,
            wetdry=wetdry,
            save_specific_discharge=True,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=H1)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

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


# Function to write flow diversion model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the model. True is returned if the model runs successfully.
#


def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)

    return success


# Function to create a figure


def create_figure():
    fig = plt.figure(figsize=figure_size, constrained_layout=False)
    gs = mpl.gridspec.GridSpec(ncols=10, nrows=7, figure=fig, wspace=5)
    plt.axis("off")

    # create axes
    ax1 = fig.add_subplot(gs[:5, :])
    ax2 = fig.add_subplot(gs[5:, :])

    # set limits for map figure
    ax1.set_xlim(extents[:2])
    ax1.set_ylim(extents[2:])
    ax1.set_aspect("equal")

    # set limits for legend area
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)

    # get rid of ticks and spines for legend area
    ax2.axis("off")
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines["top"].set_color("none")
    ax2.spines["bottom"].set_color("none")
    ax2.spines["left"].set_color("none")
    ax2.spines["right"].set_color("none")
    ax2.patch.set_alpha(0.0)

    axes = [ax1, ax2]

    return fig, axes


# Function to plot the grid


def plot_grid(gwf, silent=True):
    verbose = not silent
    fs = USGSFigure(figure_type="map", verbose=verbose)

    bot = gwf.dis.botm.array

    fig, axes = create_figure()
    ax = axes[0]
    mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
    bot_coll = mm.plot_array(bot, cmap="viridis", vmin=bmin, vmax=bmax)
    mm.plot_bc("CHD", color="cyan")
    cv = mm.contour_array(
        bot,
        levels=blevels,
        linewidths=0.5,
        linestyles=":",
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
        ls=":",
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
            "{}-grid{}".format(sim_name, config.figure_ext),
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

        # create MODFLOW 6 head object
        file_name = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        hobj = flopy.utils.HeadFile(fpth)

        # create MODFLOW 6 cell-by-cell budget object
        file_name = gwf.oc.budget_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        cobj = flopy.utils.CellBudgetFile(fpth, precision="double")

        # extract heads and specific discharge
        head = hobj.get_data(totim=1.0)
        imask = (head > -1e30) & (head <= bot)
        head[imask] = -1e30  # botm[imask]
        spdis = cobj.get_data(text="DATA-SPDIS", totim=1.0)

        # Create figure for simulation
        fig, axes = create_figure()

        ax = axes[0]
        mm = flopy.plot.PlotMapView(gwf, ax=ax, extent=extents)
        if bot.max() < 20:
            cv = mm.contour_array(
                bot,
                levels=blevels,
                linewidths=0.5,
                linestyles=":",
                colors=bcolor,
                zorder=9,
            )
            plt.clabel(cv, fmt="%1.0f", zorder=9)
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
        mm.plot_specific_discharge(
            spdis, normalize=True, color="0.75", zorder=11
        )
        ax.set_xlabel("x-coordinate, in meters")
        ax.set_ylabel("y-coordinate, in meters")
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
            lw=0,
            marker=u"$\u2192$",
            ms=10,
            mfc="0.75",
            mec="0.75",
            label="Normalized specific discharge",
        )
        if bot.max() < 20:
            ax.plot(
                -10000,
                -10000,
                lw=0.5,
                ls=":",
                color=bcolor,
                label="Bottom elevation contour, m",
            )
        ax.plot(
            -10000,
            -10000,
            lw=0.5,
            ls="-",
            color=vcolor,
            label="Head contour, m",
        )
        fs.graph_legend(ax, loc="center", ncol=2)

        cax = plt.axes([0.275, 0.125, 0.45, 0.025])
        cbar = plt.colorbar(
            h_coll, shrink=0.8, orientation="horizontal", cax=cax
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
    # ### Zaidel Simulation
    #
    # Simulated heads in the flow diversion model with Newton-Raphson.

    simulation(0)

    # Simulated heads in the flow diversion model with rewetting.

    simulation(1)

    # Simulated heads in the flow diversion model with Newton-Raphson and
    # cylinderical topography.

    simulation(2)
