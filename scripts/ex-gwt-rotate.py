# ## Rotating Interface Problem
#
# Density driven groundwater flow
#
#


# ### Rotating Interface Problem Setup

# Imports

import os
import sys
import matplotlib.pyplot as plt
import flopy
import numpy as np

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# Import common functionality

import config
from figspecs import USGSFigure
from analytical import BakkerRotatingInterface

mf6exe = os.path.abspath(config.mf6_exe)
exe_name_mf = config.mf2005_exe
exe_name_mt = config.mt3dms_exe

# Set figure properties specific to this problem

figure_size = (6, 4)

# Base simulation and model name and workspace

ws = config.base_ws
example_name = "ex-gwt-rotate"

# Model units

length_units = "m"
time_units = "days"

# Table of model parameters

nper = 1  # Number of periods
nstp = 1000  # Number of time steps
perlen = 10000  # Simulation time length ($d$)
length = 300  # Length of box ($m$)
height = 40.0  # Height of box  ($m$)
nlay = 80  # Number of layers
nrow = 1  # Number of rows
ncol = 300  # Number of columns
system_length = 150.0  # Length of system ($m$)
delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
delv = 0.5  # Layer thickness
top = height / 2  # Top of the model ($m$)
hydraulic_conductivity = 2.0  # Hydraulic conductivity ($m d^{-1}$)
denseref = 1000.0  # Reference density
denseslp = 0.7  # Density and concentration slope
rho1 = 1000.0  # Density of zone 1 ($kg m^3$)
rho2 = 1012.5  # Density of zone 2 ($kg m^3$)
rho3 = 1025.0  # Density of zone 3 ($kg m^3$)
c1 = 0.0  # Concentration of zone 1 ($kg m^3$)
c2 = 17.5  # Concentration of zone 2 ($kg m^3$)
c3 = 35  # Concentration of zone 3 ($kg m^3$)
a1 = 40.0  # Interface extent for zone 1 and 2
a2 = 40.0  # Interface extent for zone 2 and 3
b = height / 2.0
x1 = 170.0  # X-midpoint location for zone 1 and 2 interface
x2 = 130.0  # X-midpoint location for zone 2 and 3 interface
porosity = 0.2  # Porosity (unitless)

botm = [top - k * delv for k in range(1, nlay + 1)]

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-8, 1e-8, 0.97


# ### Functions to build, write, run, and plot models
#
# MODFLOW 6 flopy GWF simulation object (sim) is returned
#


def get_cstrt(nlay, ncol, length, x1, x2, a1, a2, b, c1, c2, c3):
    cstrt = c1 * np.ones((nlay, ncol), dtype=float)
    from flopy.utils.gridintersect import GridIntersect
    from shapely.geometry import Polygon

    p3 = Polygon([(0, b), (x2 - a2, b), (x2 + a2, 0), (0, 0)])
    p2 = Polygon([(x2 - a2, b), (x1 - a1, b), (x1 + a1, 0), (x1 - a1, 0)])
    delc = b / nlay * np.ones(nlay)
    delr = length / ncol * np.ones(ncol)
    sgr = flopy.discretization.StructuredGrid(delc, delr)
    ix = GridIntersect(sgr, method="structured")
    for ival, p in [(c2, p2), (c3, p3)]:
        result = ix.intersect(p)
        for i, j in list(result["cellids"]):
            cstrt[i, j] = ival
    return cstrt


def build_model(sim_folder):
    print("Building model...{}".format(sim_folder))
    name = "flow"
    sim_ws = os.path.join(ws, sim_folder)
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=sim_ws,
        exe_name=config.mf6_exe,
        continue_=True,
    )
    tdis_ds = ((perlen, nstp, 1.0),)
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )
    sim.register_ims_package(ims, [gwf.name])
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
        save_specific_discharge=True,
        icelltype=0,
        k=hydraulic_conductivity,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    pd = [(0, denseslp, 0.0, "trans", "concentration")]
    flopy.mf6.ModflowGwfbuy(gwf, denseref=denseref, packagedata=pd)
    head_filerecord = "{}.hds".format(name)
    budget_filerecord = "{}.bud".format(name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwt.name),
    )
    sim.register_ims_package(imsgwt, [gwt.name])
    flopy.mf6.ModflowGwtdis(
        gwt,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    cstrt = get_cstrt(nlay, ncol, length, x1, x2, a1, a2, height, c1, c2, c3)
    flopy.mf6.ModflowGwtic(gwt, strt=cstrt)
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord="{}.cbc".format(gwt.name),
        concentration_filerecord="{}.ucn".format(gwt.name),
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )
    return sim


# Function to write model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)
    return


# Function to run the model
# True is returned if the model runs successfully


def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success = False
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)
    return success


# Function to plot the model results


def plot_velocity_profile(sim, idx):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = example_name
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model("flow")
    gwt = sim.get_model("trans")
    print("Creating velocity profile plot...")

    # find line of cells on left side of first interface
    cstrt = gwt.ic.strt.array
    cstrt = cstrt.reshape((nlay, ncol))
    interface_coords = []
    for k in range(nlay):
        crow = cstrt[k]
        j = (np.abs(crow - c2)).argmin() - 1
        interface_coords.append((k, j))

    # plot velocity
    xc, yc, zc = gwt.modelgrid.xyzcellcenters
    xg = []
    zg = []
    for k, j in interface_coords:
        x = xc[0, j]
        z = zc[k, 0, j]
        xg.append(x)
        zg.append(z)
    xg = np.array(xg)
    zg = np.array(zg)

    # set up plot
    fig = plt.figure(figsize=(4, 6))
    ax = fig.add_subplot(1, 1, 1)

    # plot analytical solution
    qx1, qz1 = BakkerRotatingInterface.get_w(
        xg, zg, hydraulic_conductivity, rho1, rho2, a1, b, x1
    )
    qx2, qz2 = BakkerRotatingInterface.get_w(
        xg, zg, hydraulic_conductivity, rho2, rho3, a2, b, x2
    )
    qx = qx1 + qx2
    qz = qz1 + qz2
    vh = qx + qz * a1 / b
    vh = vh / porosity
    ax.plot(vh, zg, "k-")

    # plot numerical results
    file_name = gwf.oc.budget_filerecord.get_data()[0][0]
    fpth = os.path.join(sim_ws, file_name)
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    kstpkper = bobj.get_kstpkper()
    spdis = bobj.get_data(text="DATA-SPDIS", kstpkper=kstpkper[0])[0]
    qxsim = spdis["qx"].reshape((nlay, ncol))
    qzsim = spdis["qz"].reshape((nlay, ncol))
    qx = []
    qz = []
    for k, j in interface_coords:
        qx.append(qxsim[k, j])
        qz.append(qzsim[k, j])
    qx = np.array(qx)
    qz = np.array(qz)
    vh = qx + qz * a1 / b
    vh = vh / porosity
    ax.plot(vh, zg, "bo", mfc="none")

    # configure plot and save
    ax.plot([0, 0], [-b, b], "k--", linewidth=0.5)
    ax.set_xlim(-0.1, 0.1)
    ax.set_ylim(-b, b)
    ax.set_ylabel("z location of left interface (m)")
    ax.set_xlabel("$v_h$ (m/d) of left interface at t=0")
    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-vh{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)


def plot_conc(sim, idx):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = example_name
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model("flow")
    gwt = sim.get_model("trans")

    # make initial conditions figure
    print("Creating initial conditions figure...")
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 0})
    pxs.plot_array(gwt.ic.strt.array, cmap="jet", vmin=c1, vmax=c3)
    pxs.plot_grid(linewidth=0.1)
    ax.set_ylabel("z position (m)")
    ax.set_xlabel("x position (m)")
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-bc{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    plt.close("all")

    # make results plot
    print("Making results plot...")
    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    # create MODFLOW 6 head object
    fpth = os.path.join(sim_ws, "trans.ucn")
    cobj = flopy.utils.HeadFile(fpth, text="concentration")
    times = cobj.get_times()
    times = np.array(times)

    # plot times in the original publication
    plot_times = [
        2000.0,
        10000.0,
    ]

    nplots = len(plot_times)
    for iplot in range(nplots):

        time_in_pub = plot_times[iplot]
        idx_conc = (np.abs(times - time_in_pub)).argmin()
        time_this_plot = times[idx_conc]
        conc = cobj.get_data(totim=time_this_plot)

        ax = fig.add_subplot(2, 1, iplot + 1)
        pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 0})
        pxs.plot_array(conc, cmap="jet", vmin=c1, vmax=c3)
        ax.set_xlim(0, length)
        ax.set_ylim(-height / 2.0, height / 2.0)
        ax.set_ylabel("z position (m)")
        ax.set_xlabel("x position (m)")
        ax.set_title("Time = {} days".format(time_this_plot))
    plt.tight_layout()

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-conc{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


def make_animated_gif(sim, idx):
    from matplotlib.animation import FuncAnimation, PillowWriter

    print("Creating animation...")
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = example_name
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model("flow")
    gwt = sim.get_model("trans")

    fpth = os.path.join(sim_ws, "trans.ucn")
    cobj = flopy.utils.HeadFile(fpth, text="concentration")
    times = cobj.get_times()
    times = np.array(times)
    conc = cobj.get_alldata()

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 0})
    pc = pxs.plot_array(conc[0], cmap="jet", vmin=c1, vmax=c3)

    def init():
        ax.set_xlim(0, length)
        ax.set_ylim(-height / 2, height / 2)
        ax.set_title("Time = {} seconds".format(times[0]))

    def update(i):
        pc.set_array(conc[i].flatten())
        ax.set_title("Time = {} days".format(times[i]))

    ani = FuncAnimation(
        fig, update, range(1, times.shape[0], 5), init_func=init
    )
    writer = PillowWriter(fps=50)
    fpth = os.path.join("..", "figures", "{}{}".format(sim_name, ".gif"))
    ani.save(fpth, writer=writer)
    return


def plot_results(sim, idx):
    if config.plotModel:
        plot_conc(sim, idx)
        plot_velocity_profile(sim, idx)
        if config.plotSave:
            make_animated_gif(sim, idx)
    return


# Function that wraps all of the steps for each scenario
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def scenario(idx, silent=True):
    sim = build_model(example_name)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(sim, idx)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Rotating Interface Problem

    # Plot showing MODFLOW 6 results

    scenario(0)
