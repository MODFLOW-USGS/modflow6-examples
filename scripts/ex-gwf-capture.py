# ## Capture fraction analysis
#
# This problem is an example of a capture fraction analysis based on
# Leake and others (2010) using the model developed by Freyberg (1988) and
# the MODFLOW API. The MODFLOW API is used because the capture fraction
# for each cell can be calculated without regenerating the input files or
# running the model to completion (finalizing the time step), which writes
# output to the listing file. The capture fraction perturbation flux is
# added to the model using the API package, which adds the specified flux to
# the right-hand side of the system of equations.
#

# ### Capture Fraction Problem Setup
#
# Imports

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy
import modflowapi

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

sim_name = "ex-gwf-capture"

# Model units

length_units = "meters"
time_units = "seconds"

# Load the bottom, hydraulic conductivity, and idomain arrays

bottom = np.loadtxt(
    os.path.join("..", "data", sim_name, "bottom.txt"),
)
k11 = np.loadtxt(
    os.path.join("..", "data", sim_name, "hydraulic_conductivity.txt"),
)
idomain = np.loadtxt(
    os.path.join("..", "data", sim_name, "idomain.txt"),
    dtype=np.int32,
)


# Table Capture Fraction Model Parameters

nper = 1  # Number of periods
nlay = 1  # Number of layers
nrow = 40  # Number of rows
ncol = 20  # Number of columns
delr = 250.0  # Column width ($m$)
delc = 250.0  # Row width ($m$)
top = 35.0  # Top of the model ($m$)
icelltype = 1  # Cell conversion type
strt = 45.0  # Starting head ($m$)
recharge = 1.60000000E-09  # Recharge rate ($m/s$)
cf_q = -1e-3 # Perturbation flux ($m/s$)

# Static temporal data used by TDIS file

tdis_ds = (
    (1.0, 1.0, 1),
)


# ### Create Capture Fraction Model Boundary Conditions

# Well boundary conditions

wel_spd = {
    0: [
        [0, 8, 15, -0.00820000],
        [0, 10, 12, -0.00410000],
        [0, 19, 13, -0.00390000],
        [0, 25, 9, -8.30000000E-04],
        [0, 28, 5, -7.20000000E-04],
        [0, 33, 11, -0.00430000],
    ]
}

# Constant head boundary conditions

chd_spd = {
    0: [
        [0, 39, 5, 16.90000000],
        [0, 39, 6, 16.40000000],
        [0, 39, 7, 16.10000000],
        [0, 39, 8, 15.60000000],
        [0, 39, 9, 15.10000000],
        [0, 39, 10, 14.00000000],
        [0, 39, 11, 13.00000000],
        [0, 39, 12, 12.50000000],
        [0, 39, 13, 12.00000000],
        [0, 39, 14, 11.40000000],
    ]
}

# River boundary conditions

rbot = np.linspace(20., 10.25, num=nrow)
rstage = np.linspace(20.1, 11.25, num=nrow)
riv_spd = []
for idx, (s, b) in enumerate(zip(rstage, rbot)):
    riv_spd.append([0, idx, 14, s, 0.05, b])
riv_spd = {0: riv_spd}

# Solver parameters

nouter = 100
ninner = 25
hclose = 1e-9
rclose = 1e-3


# Create mapping array for the capture zone analysis
imap = idomain.copy()
for (_k, i, j, _h) in chd_spd[0]:
    imap[i, j] = 0


# ### Functions to build, write, run, and plot the MODFLOW 6 Capture Zone model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model():
    if config.buildModel:
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name=config.libmf6_exe,
        )
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_ds, time_units=time_units
        )
        flopy.mf6.ModflowIms(
            sim,
            linear_acceleration="BICGSTAB",
            outer_maximum=nouter,
            outer_dvclose=hclose * 10.,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord="{} strict".format(rclose),
        )
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=sim_name,
            newtonoptions="NEWTON UNDER_RELAXATION",
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
            botm=bottom,
            idomain=idomain,
        )

        flopy.mf6.ModflowGwfnpf(
            gwf,
            icelltype=icelltype,
            k=k11,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)
        flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd, pname="RIV-1")
        flopy.mf6.ModflowGwfwel(gwf, stress_period_data=wel_spd)
        flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
        flopy.mf6.ModflowGwfapi(
            gwf,
            maxbound=1,
            pname="CF-1",
        )
        flopy.mf6.ModflowGwfoc(
            gwf,
            printrecord=[("BUDGET", "ALL"),],
        )
        return sim
    else:
        return None


# Function to write MODFLOW 6 Capture Fraction model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to solve the system of equations to convergence

def solve_current(mobj):
    max_iter = mobj.get_value(mobj.get_var_address("MXITER", "SLN_1"))

    # convergence loop
    kiter = 0
    mobj.prepare_solve()

    while kiter < max_iter:
        has_converged = mobj.solve()
        kiter += 1
        if has_converged:
            break

    mobj.finalize_solve()


# Function to get the streamflow from memory

def get_streamflow(mobj):
    tag = mobj.get_var_address("SIMVALS", sim_name, "RIV-1")
    return mobj.get_value(tag).sum()


# Function to update the API package

def update_well(mobj, node, q=-1e-3):
    tag = mobj.get_var_address("NBOUND", sim_name, "CF-1")
    nbound = mobj.get_value(tag)
    if nbound[0] < 1:
        nbound[0] = 1
        mobj.set_value(tag, nbound)
    tag = mobj.get_var_address("NODELIST", sim_name, "CF-1")
    nodelist = mobj.get_value(tag)
    nodelist[0] = node + 1 # convert from zero-based to one-based node number
    mobj.set_value(tag, nodelist)
    tag = mobj.get_var_address("RHS", sim_name, "CF-1")
    rhs = mobj.get_value(tag)
    rhs[:] = -q
    mobj.set_value(tag, rhs)
    return

# Function to run the Capture Fraction model.
# True is returned if the model runs successfully
#

@config.timeit
def run_model():
    success = True
    if config.runModel:
        sim_ws = os.path.join(ws, sim_name)
        mf6 = modflowapi.ModflowApi(config.libmf6_exe, working_directory=sim_ws)
        mf6.initialize()
        # get dt and prepare for non-linear iterations
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)
        # run the base simulation
        solve_current(mf6)
        # get the base streamflow term
        qbase = get_streamflow(mf6)
        # create capture fraction array
        capture = np.zeros((nrow, ncol), dtype=float)
        # iterate through each active cell
        ireduced_node = -1
        for irow in range(nrow):
            for jcol in range(ncol):

                # skip inactive cells
                if imap[irow, jcol] < 1:
                    continue

                # increment node number
                ireduced_node += 1

                # update the api package with the well
                update_well(mf6, ireduced_node, q=cf_q)

                # solve with the updated well
                solve_current(mf6)

                # calculate the total streamflow
                rivcf = get_streamflow(mf6)

                # process the results
                cf = (rivcf - qbase) / abs(cf_q)

                # add the value to the capture array
                capture[irow, jcol] = cf

        # finalize libmf6
        mf6.finalize()

        # save the capture fraction array
        fpth = os.path.join(sim_ws, "capture.npz")
        np.savez_compressed(fpth, capture=capture)

    return success


# Function to plot the Capture Fraction model results with heads in each layer.
#

def plot_results(silent=True):
    if config.plotModel:
        verbose = not silent
        if silent:
            verbosity_level = 0
        else:
            verbosity_level = 1


        fs = USGSFigure(figure_type="map", verbose=verbose)
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation.load(
            sim_name=sim_name, sim_ws=sim_ws, verbosity_level=verbosity_level
        )
        gwf = sim.get_model(sim_name)

        # load the capture fraction data
        fpth = os.path.join(sim_ws, "capture.npz")
        capture = np.load(fpth)["capture"]

        # plot grid
        fig = plt.figure(figsize=(4, 3.75), constrained_layout=True)
        gs = mpl.gridspec.GridSpec(
            2,
            2,
            figure=fig,
            width_ratios=(4, 1),
            height_ratios=(1, 6),
        )

        ax = fig.add_subplot(gs[:, 0])
        ax.set_aspect("equal")

        # ax = fig.add_subplot(1, 1, 1)
        # ax.set_aspect("equal")
        mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
        cf = mm.plot_array(capture, vmin=0, vmax=1)
        mm.plot_grid(lw=0.5, color="0.5")
        mm.plot_bc('WEL')
        ax.axvline(x=14.5 * delc, lw=1.25, color="cyan")
        mm.plot_bc('CHD', color="green")
        mm.plot_ibound()
        ax.set_ylabel("y-coordinate, in feet")
        ax.set_xlabel("x-coordinate, in feet")
        fs.remove_edge_ticks(ax)

        ax = fig.add_subplot(gs[0, 1])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_color("none")
        ax.spines["bottom"].set_color("none")
        ax.spines["left"].set_color("none")
        ax.spines["right"].set_color("none")
        ax.patch.set_alpha(0.0)
        cbar = plt.colorbar(cf, ax=ax, orientation="horizontal")
        cbar.ax.set_xlabel("Streamflow capture fraction")
        ax.plot(
            -1000,
            -1000,
            "s",
            ms=5,
            color="green",
            mec="black",
            mew=0.5,
            label="Constant head",
        )
        ax.plot(
            -1000,
            -1000,
            color="cyan",
            lw=1.25,
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
            color="black",
            mec="black",
            mew=0.5,
            label="Inactive cell",
        )
        fs.graph_legend(
            ax,
            ncol=1,
            frameon=False,
            loc="upper center",
        )

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}-01{}".format(sim_name, config.figure_ext)
            )
            fig.savefig(fpth)


# Function that wraps all of the steps for the TWRI model
#
# 1. build_model,
# 2. write_model, and
# 3. run_model
# 4. plot_results.
#


def simulation(silent=True):

    sim = build_model()

    write_model(sim, silent=silent)

    success = run_model()

    assert success, "could not run...{}".format(sim_name)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(silent=False)


def test_plot():
    plot_results(silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Capture Zone Simulation
    #
    #  Capture zone examples using the MODFLOW API with the Freyberg (1988) model

    simulation()

    # Simulated streamflow capture fraction map for the Freyberg (1988) groundwater
    # flow model.

    plot_results()
