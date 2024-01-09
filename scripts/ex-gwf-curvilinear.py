# ## Curvilinear example
#
# This example, ex-gwf-curvilinear, shows how the MODFLOW 6 DISV Package
# can be used to simulate a multipart curvilinear models.
#
# The example reproduces the hypothetical model grid presented in Figure 6 of
#    Romero, D. M., & Silver, S. E. (2006).
#    Grid cell distortion and MODFLOW's integrated finite difference
#    numerical solution. Groundwater, 44(6), 797-802.
#
# The hypothetical, curvilinear grid is built in three parts:
#    1) 180 to 270 degree curvilinear grid
#    2) 16 by 18 structured grid
#    3) 90 to 0 degree curvilinear grid
# that are merged, as 1-2-3, to make the final multipart curvilinear grid.

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
from flopy.plot.styles import styles

from DisvCurvilinearBuilder import DisvCurvilinearBuilder
from DisvStructuredGridBuilder import DisvStructuredGridBuilder
from DisvGridMerger import DisvGridMerger

# Set default figure properties

figure_size_grid_com = (6.5, 2.5)
figure_size_grid = (6.5, 3)
figure_size_head = (6.5, 2.5)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-curvilin"

# Model units

length_units = "feet"
time_units = "days"

# Table Model Parameters

_ = "Steady-State"  # Simulation Type
nper = 1  # Number of periods
_ = 1  # Number of time steps

nlay = 1  # Number of layers
_ = 864  # Number cells per layer

surface_elevation = 10.0  # Top of the model ($ft$)
model_base = 0.0  # Base of the model ($ft$)

Tran = 0.19  # Horizontal transmissivity ($ft^2/day$)
k11 = 0.019  # Horizontal hydraulic conductivity ($ft/day$)

bc0 = 10  # Left constant head boundary ($ft$)
_ = "3.334"  # Right constant head boundary ($ft$)

_ = " "  # --- Left Curvilinear Grid Properties ---

_ = "180"  # Degree angle of column 1 boundary
_ = "270"  # Degree angle of column ncol boundary
_ = "5"  # Degree angle width of each column

nradial1 = 16  # Number of radial direction cells (radial bands)
_ = 18  # Number of columns in radial band (ncol)

r_inner1 = 4  # Grid inner radius ($ft$)
r_outer1 = 20  # Grid outer radius ($ft$)
r_width1 = 1  # Radial band width ($ft$)

_ = " "  # --- Middle Structured Grid Properties ---
nrow = 16  # Number of rows
ncol = 18  # Number of columns
row_width = 1  # Row width ($ft$)
col_width = 1  # Column width ($ft$)

_ = " "  # --- Right Curvilinear Grid Properties ---

_ = "0"  # Degree angle of column 1 boundary
_ = "90"  # Degree angle of column ncol boundary
_ = "5"  # Degree angle width of each column

nradial2 = 16  # Number of radial direction cells (radial bands)
_ = 18  # Number of columns in radial band (ncol)

r_inner2 = 4  # Grid inner radius ($ft$)
r_outer2 = 20  # Grid outer radius ($ft$)
r_width2 = 1  # Grid radial band width ($ft$)

# Set up input that is not used in the table
# Left Curvilinear Model Angle and discretization
angle_start1 = 180
angle_stop1 = 270
angle_step1 = 5

# Right Curvilinear Model Angles
angle_start2 = 0
angle_stop2 = 90
angle_step2 = 5


# Right Curvilinear Model Boundary Condition
bc1 = bc0 / 3

# Radius for each radial band.
#   First value is inner radius, the remaining are outer radii

radii = np.arange(r_inner1, r_outer1 + r_width1, r_width1)

# Get the curvilinear model properties and vertices

# Left Curvilinear Model
curvlin1 = DisvCurvilinearBuilder(
    nlay,
    radii,
    angle_start1,
    angle_stop1,
    angle_step1,
    surface_elevation=surface_elevation,
    layer_thickness=surface_elevation,
    single_center_cell=False,
    origin_x=radii[-1],  # Shift to make merged image have (0, 0) for origin
    origin_y=radii[-1] + radii[0],
)

# Middle Structured Grid Model
rectgrid = DisvStructuredGridBuilder(
    nlay,
    nrow,
    ncol,
    row_width,
    col_width,
    surface_elevation,
    surface_elevation,
)

# Right Curvilinear Model
curvlin2 = DisvCurvilinearBuilder(
    nlay,
    radii,
    angle_start2,
    angle_stop2,
    angle_step2,
    surface_elevation=surface_elevation,
    layer_thickness=surface_elevation,
    single_center_cell=False,
)

# Combine the three models into one new vertex grid

grid_merger = DisvGridMerger()

grid_merger.add_grid("curvlin1", curvlin1)
grid_merger.add_grid("rectgrid", rectgrid)
grid_merger.add_grid("curvlin2", curvlin2)

# # Plot individual grids to find vertex connections
# grid_merger.plot_grid("curvlin1", show=False)
# grid_merger.plot_grid("rectgrid", show=False)
# grid_merger.plot_grid("curvlin2", show=False)
# plt.show()

# Setup vertex connections between model grids
grid_merger.set_vertex_connection("curvlin1", "rectgrid", 19 - 1, 1 - 1)
grid_merger.set_vertex_connection("rectgrid", "curvlin2", 19 - 1, 323 - 1)

# Merge grids into one single model grid
grid_merger.merge_grids()

# Shift first curvilinear grid for plotting against the orgin.
# (Note, grid_merger no longer needs curvlin1)
curvlin1.change_origin(0.0, 0.0)

# grid_merger.plot_grid(show=False, figsize=(23, 10))

# Constant head boundary condition
# Constant head is located along column 1 of curvlin1
# and column 1 of curvlin2
chd_left = []
chd_right = []
for lay in range(nlay):
    for cellid_old in curvlin1.iter_column_cellid(col=0):
        node = grid_merger.get_merged_cell2d("curvlin1", cellid_old)
        chd_left.append([(lay, node), bc0])

for lay in range(nlay):
    for cellid_old in curvlin2.iter_column_cellid(col=0):
        node = grid_merger.get_merged_cell2d("curvlin2", cellid_old)
        chd_right.append([(lay, node), bc1])

chd_left = {sp: chd_left for sp in range(nper)}

chd_right = {sp: chd_right for sp in range(nper)}

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

        disv = flopy.mf6.ModflowGwfdisv(
            gwf, length_units=length_units, **grid_merger.get_disv_kwargs()
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
            stress_period_data=chd_left,
            pname="CHD-LEFT",
            filename=f"{sim_name}.left.chd",
            save_flows=True,
        )
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chd_right,
            pname="CHD-RIGHT",
            filename=f"{sim_name}.right.chd",
            save_flows=True,
        )

        flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=f"{name}.cbc",
            head_filerecord=f"{name}.hds",
            headprintrecord=[
                (
                    "COLUMNS",
                    curvlin1.ncol + ncol + curvlin2.ncol,
                    "WIDTH",
                    15,
                    "DIGITS",
                    6,
                    "GENERAL",
                )
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
    with styles.USGSMap() as fs:
        gwf = sim.get_model(sim_name)

        fig = plt.figure(figsize=figure_size_grid)

        ax = fig.add_subplot(1, 1, 1, aspect="equal")
        pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
        pmv.plot_grid()
        pmv.plot_bc(name="CHD-LEFT", alpha=0.75, color="blue")
        pmv.plot_bc(name="CHD-RIGHT", alpha=0.75, color="blue")
        ax.set_xlabel("x position (ft)")
        ax.set_ylabel("y position (ft)")
        for i, (x, y) in enumerate(
            zip(gwf.modelgrid.xcellcenters, gwf.modelgrid.ycellcenters)
        ):
            ax.text(
                x,
                y,
                f"{i + 1}",
                fontsize=3,
                horizontalalignment="center",
                verticalalignment="center",
            )
        v = gwf.disv.vertices.array
        vert_size = 2
        ax.plot(v["xv"], v["yv"], "yo", markersize=vert_size)
        for i in range(v.shape[0]):
            x, y = v["xv"][i], v["yv"][i]
            ax.text(
                x,
                y,
                f"{i + 1}",
                fontsize=vert_size,
                color="red",
                horizontalalignment="center",
                verticalalignment="center",
            )

        fig.tight_layout()

        # Save components that made up the main grid
        fig2, ax2 = plt.subplots(
            1,
            3,
            figsize=figure_size_grid_com,
        )

        curvlin1.plot_grid(
            "Left Curvilinear Grid",
            ax_override=ax2[0],
            cell_dot=False,
            cell_num=False,
            vertex_dot=True,
            vertex_num=False,
            vertex_dot_size=3,
            vertex_dot_color="y",
        )

        rectgrid.plot_grid(
            "Center Rectangular Grid",
            ax_override=ax2[1],
            cell_dot=False,
            cell_num=False,
            vertex_dot=True,
            vertex_num=False,
            vertex_dot_size=3,
            vertex_dot_color="y",
        )

        curvlin2.plot_grid(
            "Right Curvilinear Grid",
            ax_override=ax2[2],
            cell_dot=False,
            cell_num=False,
            vertex_dot=True,
            vertex_num=False,
            vertex_dot_size=3,
            vertex_dot_color="y",
        )

        for ax_tmp in ax2:
            ax_tmp.set_xlabel("x position (ft)")
            ax_tmp.set_ylabel("y position (ft)")

        xshift, yshift = 0.0, 0.0
        for ax_tmp in ax2:
            xmin, xmax = ax_tmp.get_xlim()
            ymin, ymax = ax_tmp.get_ylim()
            if xshift < xmax - xmin:
                xshift = xmax - xmin
            if yshift < ymax - ymin:
                yshift = ymax - ymin

        for ax_tmp in ax2:
            xmin, xmax = ax_tmp.get_xlim()
            ymin, ymax = ax_tmp.get_ylim()
            ax_tmp.set_xlim(xmin, xmin + xshift)
            ax_tmp.set_ylim(ymin, ymin + yshift)

        ax2[0].annotate(
            "A",
            (-0.05, 1.05),
            xycoords="axes fraction",
            fontweight="black",
            fontsize="xx-large",
        )

        ax2[1].annotate(
            "B",
            (-0.05, 1.05),
            xycoords="axes fraction",
            fontweight="black",
            fontsize="xx-large",
        )

        ax2[2].annotate(
            "C",
            (-0.05, 1.05),
            xycoords="axes fraction",
            fontweight="black",
            fontsize="xx-large",
        )

        fig2.tight_layout()

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..",
                "figures",
                f"{sim_name}-grid{config.figure_ext}",
            )
            fig.savefig(fpth, dpi=600)

            fpth2 = os.path.join(
                "..",
                "figures",
                f"{sim_name}-grid-components{config.figure_ext}",
            )
            fig2.savefig(fpth2, dpi=300)


# Function to plot the curvilinear model results.


def plot_head(sim):
    with styles.USGSMap() as fs:
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
            fig.savefig(fpth, dpi=300)


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
    return


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
