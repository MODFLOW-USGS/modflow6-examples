# ## Curvilinear example
#
# This example, ex-gwf-curvilinear, shows how the MODFLOW 6 DISV Package
# can be used to simulate a curvilinear models.
#
# The example corresponds to Figure 3d (lower-right) in:
#    Romero, D. M., & Silver, S. E. (2006).
#    Grid cell distortion and MODFLOW's integrated finite difference
#    numerical solution. Groundwater, 44(6), 797-802.
#
# And the numerical result is compared against the analytical solution
# presented in Equation 5.4 of
#    Crank, J. (1975). The mathematics of diffusion.
#    Oxford. England: Clarendon.
# The equation is transformed here to use head instead of concentration

# Imports

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import flopy
from math import sqrt

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

from disv_curvilinear_builder import disv_curvilinear_builder


def analytical_model(r1, h1, r2, h2, r):
    # Analytical model from Equation 5.4 of
    #    Crank, J. (1975). The mathematics of diffusion.
    #    Oxford. England: Clarendon.
    num = h1 * np.log(r2 / r) + h2 * np.log(r / r1)
    den = np.log(r2 / r1)
    return num / den


# Set default figure properties

figure_size_grid = (6, 6)
figure_size_head = (7, 6)
figure_size_obsv = (6, 6)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-curve-90"

# Model units

length_units = "feet"
time_units = "days"

# Table Model Parameters

_ = "Steady-State"  # Simulation Type
nper = 1  # Number of periods
_ = 1  # Number of time steps

nlay = 1  # Number of layers
nradial = 16  # Number of radial direction cells (radial bands)
ncol = 18  # Number of columns in radial band (ncol)

_ = "0"  # Degree angle of column 1 boundary
_ = "90"  # Degree angle of column ncol boundary
_ = "5"  # Degree angle width of each column

r_inner = 4  # Model inner radius ($ft$)
r_outer = 20  # Model outer radius ($ft$)
r_width = 1  # Model radial band width ($ft$)

surface_elevation = 10.0  # Top of the model ($ft$)
model_base = 0.0  # Base of the model ($ft$)

Tran = 0.19  # Horizontal transmissivity ($ft^2/day$)
k11 = 0.019  # Horizontal hydraulic conductivity ($ft/day$)

bc0 = 10  # Inner Constant Head Boundary ($ft$)
_ = "3.334"  # Outer Constant Head Boundary ($ft$)

# Input specified in table as text
bc1 = bc0 / 3

angle_start = 0
angle_stop = 90
angle_step = 5

# Radius for each radial band.
#   First value is inner radius, the remaining are outer radii

radii = np.arange(r_inner, r_outer + r_width, r_width)

# Get the curvilinear model properties and vertices

curvlin = disv_curvilinear_builder(
    nlay,
    radii,
    angle_start,
    angle_stop,
    angle_step,
    surface_elevation=surface_elevation,
    layer_thickness=surface_elevation,
    single_center_cell=False,
)

# Constant head boundary condition
# Constant head is located along the innermost radial band (rad = 0)
# and outermost radial band (rad = nradial-1)

chd_inner = []
chd_outer = []
for lay in range(nlay):
    for node in curvlin.iter_radial_nodes(rad=0):
        chd_inner.append([(lay, node), bc0])
for lay in range(nlay):
    for node in curvlin.iter_radial_nodes(rad=nradial - 1):
        chd_outer.append([(lay, node), bc1])

chd_inner = {sp: chd_inner for sp in range(nper)}

chd_outer = {sp: chd_outer for sp in range(nper)}

# Static temporal data used by TDIS file
# Simulation is steady state so setup only a one day stress period.

tdis_ds = ((1.0, 1, 1),)

# Solver parameters

nouter = 500
ninner = 300
hclose = 1e-4
rclose = 1e-4


# ### Functions to build, write, run, and plot the MODFLOW 6 Curvilinear Model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(name):
    if config.buildModel:
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation(
            sim_name=name, sim_ws=sim_ws, exe_name="mf6"
        )
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_ds, time_units=time_units
        )
        flopy.mf6.ModflowIms(
            sim,
            print_option="summary",
            complexity="complex",
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
        )

        gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

        # **curvlin is an alias for **curvlin.disv_kw
        disv = flopy.mf6.ModflowGwfdisv(
            gwf, length_units=length_units, **curvlin
        )

        npf = flopy.mf6.ModflowGwfnpf(
            gwf,
            k=k11,
            k33=k11,
            save_flows=True,
            save_specific_discharge=True,
        )

        flopy.mf6.ModflowGwfsto(
            gwf,
            iconvert=0,
            steady_state=True,
            save_flows=True,
        )

        flopy.mf6.ModflowGwfic(gwf, strt=surface_elevation)

        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chd_inner,
            pname="CHD-INNER",
            filename=f"{sim_name}.inner.chd",
            save_flows=True,
        )
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chd_outer,
            pname="CHD-OUTER",
            filename=f"{sim_name}.outer.chd",
            save_flows=True,
        )

        flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=f"{name}.cbc",
            head_filerecord=f"{name}.hds",
            headprintrecord=[
                ("COLUMNS", nradial, "WIDTH", 15, "DIGITS", 6, "GENERAL")
            ],
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            filename=f"{name}.oc",
        )

        return sim
    return None


# Function to write model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the curvilinear model.
# True is returned if the model runs successfully.


@config.timeit
def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent, report=True)
        if not success:
            print("\n".join(buff))

    return success


# Function to plot the curvilinear model grid.


def plot_grid(sim, verbose=False):
    fs = USGSFigure(figure_type="map", verbose=verbose)
    gwf = sim.get_model(sim_name)

    fig = plt.figure(figsize=figure_size_grid)

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
    pmv.plot_grid()
    pmv.plot_bc(name="CHD-INNER", alpha=0.75, color="blue")
    pmv.plot_bc(name="CHD-OUTER", alpha=0.75, color="blue")
    ax.set_xlabel("x position (ft)")
    ax.set_ylabel("y position (ft)")
    for i, (x, y) in enumerate(
        zip(gwf.modelgrid.xcellcenters, gwf.modelgrid.ycellcenters)
    ):
        ax.text(
            x,
            y,
            f"{i + 1}",
            fontsize=6,
            horizontalalignment="center",
            verticalalignment="center",
        )
    v = gwf.disv.vertices.array
    ax.plot(v["xv"], v["yv"], "yo")
    for i in range(v.shape[0]):
        x, y = v["xv"][i], v["yv"][i]
        ax.text(
            x,
            y,
            f"{i + 1}",
            fontsize=5,
            color="red",
            horizontalalignment="center",
            verticalalignment="center",
        )

    fig.tight_layout()

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", f"{sim_name}-grid{config.figure_ext}"
        )
        fig.savefig(fpth)
    return


# Function to plot the curvilinear model results.


def plot_head(sim):
    fs = USGSFigure(figure_type="map", verbose=False)
    gwf = sim.get_model(sim_name)

    fig = plt.figure(figsize=figure_size_head)

    head = gwf.output.head().get_data()[:, 0, :]

    # create MODFLOW 6 cell-by-cell budget object
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
        gwf.output.budget().get_data(text="DATA-SPDIS", totim=1.0)[0],
        gwf,
    )

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
    cb = pmv.plot_array(head, cmap="jet", vmin=0.0, vmax=head.max())
    pmv.plot_vector(
        qx,
        qy,
        normalize=False,
        color="0.75",
    )
    cbar = plt.colorbar(cb, shrink=0.25)
    cbar.ax.set_xlabel(r"Head, ($ft$)")
    ax.set_xlabel("x position (ft)")
    ax.set_ylabel("y position (ft)")

    fig.tight_layout()

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", f"{sim_name}-head{config.figure_ext}"
        )
        fig.savefig(fpth)
    return


def plot_analytical(sim, verbose=False):
    gwf = sim.get_model(sim_name)

    head = gwf.output.head().get_data()[:, 0, :]

    col = ncol // 2 - 1  # Get head along middle of model

    head = [head[0, curvlin.get_node(rad, col)] for rad in range(nradial)]

    xrad = [0.5 * (radii[r - 1] + radii[r]) for r in range(1, nradial + 1)]

    analytical = [head[0]]
    r1 = xrad[0]
    r2 = xrad[-1]
    h1 = bc0
    h2 = bc1
    for rad in range(2, nradial):
        r = 0.5 * (radii[rad - 1] + radii[rad])
        analytical.append(analytical_model(r1, h1, r2, h2, r))
    analytical.append(head[-1])

    fs = USGSFigure(figure_type="graph", verbose=verbose)

    obs_fig = "obs-head"
    fig = plt.figure(figsize=figure_size_obsv)
    ax = fig.add_subplot()
    ax.set_xlabel("Radial distance (ft)")
    ax.set_ylabel("Head (ft)")
    ax.plot(xrad, head, "ob", label="MF6 Solution", markerfacecolor="none")
    ax.plot(xrad, analytical, "-b", label="Analytical Solution")

    fs.graph_legend(ax)

    fig.tight_layout()

    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-{}{}".format(sim_name, obs_fig, config.figure_ext),
        )
        fig.savefig(fpth)


# Function to plot the model results.


def plot_results(silent=True):
    if not config.plotModel:
        return

    if silent:
        verbosity_level = 0
    else:
        verbosity_level = 1

    sim_ws = os.path.join(ws, sim_name)
    sim = flopy.mf6.MFSimulation.load(
        sim_name=sim_name, sim_ws=sim_ws, verbosity_level=verbosity_level
    )

    verbose = not silent

    if config.plotModel:
        plot_grid(sim, verbose)
        plot_head(sim)
        plot_analytical(sim, verbose)
    return


def calculate_model_error():
    sim_ws = os.path.join(ws, sim_name)
    sim = flopy.mf6.MFSimulation.load(
        sim_name=sim_name, sim_ws=sim_ws, verbosity_level=0
    )

    gwf = sim.get_model(sim_name)

    head = gwf.output.head().get_data()[0, 0, :]

    xrad = [0.5 * (radii[r - 1] + radii[r]) for r in range(1, nradial + 1)]

    analytical = [head[0]]
    r1 = xrad[0]
    r2 = xrad[-1]
    h1 = bc0
    h2 = bc1
    for rad in range(2, nradial):
        r = 0.5 * (radii[rad - 1] + radii[rad])
        analytical.append(analytical_model(r1, h1, r2, h2, r))
    analytical.append(head[-1])

    dim = len(head)
    rel = 0.0
    sse = 0.0
    for rad in range(nradial):
        asol = analytical[rad]
        for node in curvlin.iter_radial_nodes(rad):
            diff = head[node] - asol
            rel += abs(diff / asol)
            sse += diff**2
    # for x, y in zip(head, analytical):
    #     mae += abs(x-y)
    #     sse += (x-y)**2
    rel /= dim
    rmse = sqrt(sse / dim)
    return rel, rmse


def check_model_error():
    if config.runModel:
        rel_error, rmse = calculate_model_error()
        assert rel_error < 0.001


# Function that wraps all of the steps for the curvilinear model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(silent=True):
    # key = list(parameters.keys())[idx]
    # params = parameters[key].copy()

    sim = build_model(sim_name)

    write_model(sim, silent=silent)

    success = run_model(sim, silent=silent)
    assert success, "could not run...{}".format(sim_name)


# nosetest - exclude block from this nosetest to the next nosetest
def test_and_plot():
    simulation(silent=False)
    plot_results(silent=False)
    return


# nosetest end


if __name__ == "__main__":
    # ### Curvilinear Example

    # MF6 Curvilinear Model
    simulation()

    # Solve analytical and plot results with MF6 results
    plot_results()

    check_model_error()
