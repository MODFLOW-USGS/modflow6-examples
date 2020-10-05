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

figure_size = (6.3, 6.3)
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
    "ex-gwf-nwt-p02b": {
        "rewet": True,
        "wetfct": 0.5,
        "iwetit": 1,
        "ihdwet": 1,
        "wetdry": -0.5,
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
H1 = 25.  # Constant head along left and lower edges and starting head ($ft$)
rech = 0.05 # Recharge rate ($ft/day$)


# Static temporal data used by TDIS file

tdis_ds = (
    (190.0, 10, 1.0),
    (518.0, 2, 1.0),
    (1921.0, 17, 1.0),
    (1.0, 1, 1.0),
)

# Calculate extents, and shape3d

extents = (0, delr * ncol, 20, 65)
shape3d = (nlay, nrow, ncol)

# Create the bottom

botm = np.arange(65., -5., -5.)

# Create icelltype (which is the same as iconvert)

icelltype = 9 * [1] + 5 * [0]

# Constant head boundary conditions

chd_spd = []
for k in range(9, nlay, 1):
    chd_spd += [[k, i, ncol-1, H1] for i in range(nrow-1)]
    chd_spd += [[k, nrow - 1, j, H1] for j in range(ncol)]


# Recharge boundary conditions

rch_spd = []
for i in range(0, 2, 1):
    for j in range(0, 2, 1):
        rch_spd.append([0, i, j, rech])

# Solver parameters

nouter = 500
ninner = 100
hclose = 1e-6
rclose = 1000.

# ### Functions to build, write, run, and plot the model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(
    name,
    newton=False,
    rewet=False,
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
            newtonoptions = " "
            no_ptc = "ALL"
            complexity = "complex"
        else:
            newtonoptions = None
            no_ptc = None
            complexity = "simple"

        flopy.mf6.ModflowIms(
            sim,
            complexity=complexity,
            print_option="SUMMARY",
            no_ptcrecord=no_ptc,
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
            wetdry = 9 * [wetdry] + 5 * [0]
        else:
            rewet_record = None
        flopy.mf6.ModflowGwfnpf(
            gwf,
            rewet_record=rewet_record,
            icelltype=icelltype,
            k=k11,
            k33=k33,
            wetdry=wetdry,
        )
        flopy.mf6.ModflowGwfsto(gwf, iconvert=icelltype, ss=ss, sy=sy,
                                steady_state={3:True},)
        flopy.mf6.ModflowGwfic(gwf, strt=H1)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
        flopy.mf6.ModflowGwfrch(gwf, stress_period_data=rch_spd)

        head_filerecord = "{}.hds".format(sim_name)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            saverecord=[("HEAD", "LAST")],
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

# Create a water-table array

def get_water_table(h, bot):
    imask = (h > -1e30) & (h <= bot)
    h[imask] = -1e30
    return np.amax(h, axis=0)


# Function to plot the model results.


def plot_results(silent=True):
    verbose = not silent
    if verbose:
        verbosity_level = 1
    else:
        verbosity_level = 0

    if config.plotModel:
        fs = USGSFigure(figure_type="map", verbose=verbose)

        # load the newton model
        name = list(parameters.keys())[0]
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=sim_ws,
                                          verbosity_level=verbosity_level)
        gwf = sim.get_model(sim_name)
        bot = gwf.dis.botm.array
        xnode = gwf.modelgrid.xcellcenters[0, :]

        # create MODFLOW 6 head object
        file_name = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        hobj = flopy.utils.HeadFile(fpth)

        # get a list of times
        times = hobj.get_times()

        # load rewet model
        name = list(parameters.keys())[1]
        sim_ws = os.path.join(ws, name)

        # create MODFLOW 6 head object
        file_name = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        hobj1 = flopy.utils.HeadFile(fpth)

         # Create figure for simulation
        fig, axes = plt.subplots(ncols=1, nrows=4, sharex=True,
                                figsize=figure_size, constrained_layout=False)

        # plot the results
        for idx, ax in enumerate(axes):

            # extract heads and specific discharge for newton model
            head = hobj.get_data(totim=times[idx])
            head = get_water_table(head, bot)

            # extract heads and specific discharge for newton model
            head1 = hobj1.get_data(totim=times[idx])
            head1 = get_water_table(head1, bot)

            # calculate mean error
            diff = np.abs(head - head1)
            # print("max", diff.max(), np.argmax(diff))
            me = diff.sum() / float(ncol * nrow)
            me_text = "Mean absolute water-table error {:.3f} feet".format(me)

            ax.set_xlim(extents[:2])
            ax.set_ylim(extents[2:])
            mm = flopy.plot.PlotCrossSection(model=gwf, ax=ax,
                                             extent=extents, line={"row": 1})
            mm.plot_bc("CHD", color="cyan")
            mm.plot_grid(lw=0.5)
            ax.plot(xnode, head[0, :], lw=0.75, color="black", label="Newton-Raphson")
            ax.plot(xnode, head1[0, :], lw=0, marker="o", ms=4, mfc="none", mec="blue",
                    label="Rewetting")
            if idx == 0:
                ax.plot(-1000, -1000, lw=0, marker="s", ms=4, mec="0.5", mfc="none",
                        label="Model cell")
                ax.plot(-1000, -1000, lw=0, marker="s", ms=4, mec="0.5", mfc="cyan",
                        label="Constant head")
                fs.graph_legend(ax, loc="upper right", ncol=2,
                                frameon=True, facecolor="white", edgecolor="none")
            letter = chr(ord("@") + idx + 1)
            fs.heading(letter=letter, ax=ax)
            fs.add_text(ax, text=me_text, x=1, y=1.01, ha="right", bold=False)
            fs.remove_edge_ticks(ax)

            # set fake y-axis label
            ax.set_ylabel(" ")

        # set fake x-axis label
        ax.set_xlabel(" ")

        ax = fig.add_subplot(1, 1, 1, frameon=False)
        ax.tick_params(
            labelcolor="none", top="off", bottom="off", left="off", right="off"
        )
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 1])
        ax.set_xlabel("x-coordinate, in feet")
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        ax.set_ylabel("Water-table elevation above arbitrary datum, in meters")
        fs.remove_edge_ticks(ax)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..",
                "figures",
                "{}-01{}".format(sim_name, config.figure_ext),
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
    assert success, "could not run...{}".format(key)



# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(0, silent=False)


def test_02():
    simulation(1, silent=False)

def test_plot_results(silent=False):
    plot_results(silent=silent)


# nosetest end

if __name__ == "__main__":
    # ### MODFLOW-NWT Problem 2 Simulation
    #
    # Newton-Raphson.

    simulation(0)

    # Rewetting.

    simulation(1)

    # Plot results

    plot_results()