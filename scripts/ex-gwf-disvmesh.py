# ## USG1DISV example
#
# Demonstration of a triangular mesh with the DISV Package to discretize a
# circular island with a radius of 1500 meters.  The model has 2 layers and
# uses 2778 vertices (NVERT) to delineate 5240 cells per layer (NCPL).
# General-head boundaries are assigned to model layer 1 for cells outside of
# a 1025 m radius circle.  Recharge is applied to the top of the model.
#

# ### USG1DISV Problem Setup
#
# Imports

import os
import pathlib as pl
from os import environ

import flopy
import flopy.utils.cvfdutil
import matplotlib.pyplot as plt
import numpy as np
from flopy.plot.styles import styles
from flopy.utils.geometry import get_polygon_area
from flopy.utils.gridintersect import GridIntersect
from modflow_devtools.misc import is_in_ci, timed
from shapely.geometry import Polygon

# Set default figure properties

figure_size = (5, 5)

# Simulation name and workspace
sim_name = "ex-gwf-disvmesh"
ws = pl.Path("../examples")
data_ws = pl.Path("../data")

# Configuration

writeModel = str(environ.get("WRITE", True)).lower() == "true"
runModel = str(environ.get("RUN", True)).lower() == "true"
plotModel = str(environ.get("PLOT", True)).lower() == "true"
plotSave = str(environ.get("SAVE", is_in_ci())).lower() == "true"
createGif = str(environ.get("GIF", False)).lower() == "true"

# Model units

length_units = "meters"
time_units = "days"

# Table USG1DISV Model Parameters

nper = 1  # Number of periods
nlay = 2  # Number of layers
top = 0.0  # Top of the model ($m$)
botm_str = "-20.0, -40.0"  # Layer bottom elevations ($m$)
strt = 0.0  # Starting head ($m$)
icelltype = 0  # Cell conversion type
k11 = 10.0  # Horizontal hydraulic conductivity ($m/d$)
k33 = 0.2  # Vertical hydraulic conductivity ($m/d$)
recharge = 4.0e-3  # Recharge rate ($m/d$)

# Static temporal data used by TDIS file
# Simulation has 1 steady stress period (1 day).

perlen = [1.0]
nstp = [1]
tsmult = [1.0]
tdis_ds = list(zip(perlen, nstp, tsmult))

# Parse strings into lists

botm = [float(value) for value in botm_str.split(",")]

# create the disv grid


def from_argus_export(fname):
    f = open(fname)
    line = f.readline()
    ll = line.split()
    ncells, nverts = ll[0:2]
    ncells = int(ncells)
    nverts = int(nverts)
    verts = np.empty((nverts, 2), dtype=float)

    # read the vertices
    f.readline()
    for ivert in range(nverts):
        line = f.readline()
        ll = line.split()
        c, iv, x, y = ll[0:4]
        verts[ivert, 0] = x
        verts[ivert, 1] = y

    # read the cell information and create iverts
    iverts = []
    for icell in range(ncells):
        line = f.readline()
        ll = line.split()
        ivlist = []
        for ic in ll[2:5]:
            ivlist.append(int(ic) - 1)
        if ivlist[0] != ivlist[-1]:
            ivlist.append(ivlist[0])
        ivlist.reverse()
        iverts.append(ivlist)

    # close file and return spatial reference
    f.close()
    return verts, iverts


# Load argus mesh and get disv grid properties

fname = os.path.join(data_ws, "ex-gwf-disvmesh", "argus.exp")
verts, iverts = from_argus_export(fname)
gridprops = flopy.utils.cvfdutil.get_disv_gridprops(verts, iverts)
cell_areas = []
for i in range(gridprops["ncpl"]):
    x = verts[iverts[i], 0]
    y = verts[iverts[i], 1]
    cell_verts = np.vstack((x, y)).transpose()
    cell_areas.append(get_polygon_area(cell_verts))

# Solver parameters

nouter = 50
ninner = 100
hclose = 1e-9
rclose = 1e-6

# ### Functions to build, write, run, and plot the MODFLOW 6 DISVMESH model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(sim_name):
    sim_ws = os.path.join(ws, sim_name)
    sim = flopy.mf6.MFSimulation(sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
    flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="bicgstab",
        outer_maximum=nouter,
        outer_dvclose=hclose,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=f"{rclose} strict",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True)
    flopy.mf6.ModflowGwfdisv(
        gwf,
        length_units=length_units,
        nlay=nlay,
        top=top,
        botm=botm,
        **gridprops,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        xt3doptions=True,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    theta = np.arange(0.0, 2 * np.pi, 0.2)
    radius = 1500.0
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    outer = [(x, y) for x, y in zip(x, y)]
    radius = 1025.0
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    hole = [(x, y) for x, y in zip(x, y)]
    p = Polygon(outer, holes=[hole])
    ix = GridIntersect(gwf.modelgrid, method="vertex", rtree=True)
    result = ix.intersect(p)
    ghb_cellids = np.array(result["cellids"], dtype=int)

    ghb_spd = []
    ghb_spd += [[0, i, 0.0, k33 * cell_areas[i] / 10.0] for i in ghb_cellids]
    ghb_spd = {0: ghb_spd}
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghb_spd,
        pname="GHB",
    )

    ncpl = gridprops["ncpl"]
    rchcells = np.array(list(range(ncpl)), dtype=int)
    rchcells[ghb_cellids] = -1
    rch_spd = [(0, rchcells[i], recharge) for i in range(ncpl) if rchcells[i] > 0]
    rch_spd = {0: rch_spd}
    flopy.mf6.ModflowGwfrch(gwf, stress_period_data=rch_spd, pname="RCH")

    head_filerecord = f"{sim_name}.hds"
    budget_filerecord = f"{sim_name}.cbc"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


# Function to write MODFLOW 6 DISVMESH model files


def write_model(sim, silent=True):
    if writeModel:
        sim.write_simulation(silent=silent)


# Function to run the FHB model.
# True is returned if the model runs successfully
#


@timed
def run_model(sim, silent=False):
    success = True
    if runModel:
        success, buff = sim.run_simulation(silent=silent, report=True)
        if not success:
            print(buff)
    return success


# Function to plot the DISVMESH model results.
#
def plot_grid(idx, sim):
    with styles.USGSMap() as fs:
        sim_ws = os.path.join(ws, sim_name)
        gwf = sim.get_model(sim_name)

        fig = plt.figure(figsize=figure_size)
        fig.tight_layout()

        ax = fig.add_subplot(1, 1, 1, aspect="equal")
        pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
        pmv.plot_grid(linewidth=1)
        pmv.plot_bc(name="GHB")
        ax.set_xlabel("x position (m)")
        ax.set_ylabel("y position (m)")

        # save figure
        if plotSave:
            fpth = os.path.join("..", "figures", f"{sim_name}-grid.png")
            fig.savefig(fpth)


def plot_head(idx, sim):
    with styles.USGSMap() as fs:
        sim_ws = os.path.join(ws, sim_name)
        gwf = sim.get_model(sim_name)

        fig = plt.figure(figsize=(7.5, 5))
        fig.tight_layout()

        head = gwf.output.head().get_data()[:, 0, :]

        # create MODFLOW 6 cell-by-cell budget object
        qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
            gwf.output.budget().get_data(text="DATA-SPDIS", totim=1.0)[0],
            gwf,
        )

        ax = fig.add_subplot(1, 2, 1, aspect="equal")
        pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
        cb = pmv.plot_array(head, cmap="jet", vmin=0.0, vmax=head.max())
        pmv.plot_vector(
            qx,
            qy,
            normalize=False,
            color="0.75",
        )
        cbar = plt.colorbar(cb, shrink=0.25)
        cbar.ax.set_xlabel(r"Head, ($m$)")
        ax.set_xlabel("x position (m)")
        ax.set_ylabel("y position (m)")
        styles.heading(ax, letter="A", heading="Layer 1")

        ax = fig.add_subplot(1, 2, 2, aspect="equal")
        pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=1)
        cb = pmv.plot_array(head, cmap="jet", vmin=0.0, vmax=head.max())
        pmv.plot_vector(
            qx,
            qy,
            normalize=False,
            color="0.75",
        )
        cbar = plt.colorbar(cb, shrink=0.25)
        cbar.ax.set_xlabel(r"Head, ($m$)")
        ax.set_xlabel("x position (m)")
        ax.set_ylabel("y position (m)")
        styles.heading(ax, letter="B", heading="Layer 2")

        # save figure
        if plotSave:
            fpth = os.path.join("..", "figures", f"{sim_name}-head.png")
            fig.savefig(fpth)


def plot_results(idx, sim, silent=True):
    if plotModel:
        if idx == 0:
            plot_grid(idx, sim)
        plot_head(idx, sim)


# Function that wraps all of the steps for the FHB model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(idx, silent=True):
    sim = build_model(sim_name)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(idx, sim, silent=silent)


# ### DISVMESH Simulation
#
# Model grid and simulated heads in the DISVMESH model

simulation(0)
