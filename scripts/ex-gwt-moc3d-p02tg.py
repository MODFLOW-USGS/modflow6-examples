# ## Three-Dimensional Steady Flow with Transport
#
# MOC3D Problem 2 simulated with a triangular grid
#
#


# ### Three-Dimensional Steady Flow with Transport Problem Setup

# Imports

import os
import sys
import matplotlib.pyplot as plt
import flopy
import flopy.utils.cvfdutil
import numpy as np

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# Import common functionality

import analytical
import config
from figspecs import USGSFigure

mf6exe = os.path.abspath(config.mf6_exe)
exe_name_mf = config.mf2005_exe
exe_name_mt = config.mt3dms_exe

# Set figure properties specific to this problem

figure_size = (6, 4)

# Base simulation and model name and workspace

ws = config.base_ws
example_name = "ex-gwt-moc3d-p02tg"

# Model units

length_units = "m"
time_units = "days"

# Table of model parameters

nper = 1  # Number of periods
nlay = 40  # Number of layers
nrow = 12  # Number of rows
ncol = 30  # Number of columns
delr = 3  # Column width ($m$)
delc = 0.5  # Row width ($m$)
delv = 0.05  # Layer thickness ($m$)
top = 0.0  # Top of the model ($m$)
bottom = -2.0  # Model bottom elevation ($m$)
velocity_x = 0.1  # Velocity in x-direction ($m d^{-1}$)
hydraulic_conductivity = 0.0125  # Hydraulic conductivity ($m d^{-1}$)
porosity = 0.25  # Porosity of mobile domain (unitless)
alpha_l = 0.6  # Longitudinal dispersivity ($m$)
alpha_th = 0.03  # Transverse horizontal dispersivity ($m$)
alpha_tv = 0.006  # Transverse vertical dispersivity ($m$)
total_time = 400.0  # Simulation time ($d$)
solute_mass_flux = 2.5  # Solute mass flux ($g d^{-1}$)
source_location = (1, 12, 8)  # Source location (layer, row, column)

botm = [-(k + 1) * delv for k in range(nlay)]
specific_discharge = velocity_x * porosity
source_location0 = tuple([idx - 1 for idx in source_location])

# ### Functions to build, write, run, and plot models
#
# MODFLOW 6 flopy GWF simulation object (sim) is returned
#


def grid_triangulator(itri, delr, delc):
    regular_grid = flopy.discretization.StructuredGrid(delc, delr)
    vertdict = {}
    icell = 0
    for i in range(nrow):
        for j in range(ncol):
            vs = regular_grid.get_cell_vertices(i, j)
            vs.reverse()
            if itri[i, j] == 0:
                vertdict[icell] = [vs[0], vs[3], vs[2], vs[1], vs[0]]
                icell += 1
            elif itri[i, j] == 1:
                vertdict[icell] = [vs[0], vs[3], vs[1], vs[0]]
                icell += 1
                vertdict[icell] = [vs[1], vs[3], vs[2], vs[1]]
                icell += 1
            elif itri[i, j] == 2:
                vertdict[icell] = [vs[0], vs[2], vs[1], vs[0]]
                icell += 1
                vertdict[icell] = [vs[0], vs[3], vs[2], vs[0]]
                icell += 1
            else:
                raise Exception("Unknown itri value: {}".format(itri[i, j]))
    verts, iverts = flopy.utils.cvfdutil.to_cvfd(vertdict)
    return verts, iverts


def cvfd_to_cell2d(verts, iverts):
    vertices = []
    for i in range(verts.shape[0]):
        x = verts[i, 0]
        y = verts[i, 1]
        vertices.append([i, x, y])
    cell2d = []
    for icell2d, vs in enumerate(iverts):
        points = [tuple(verts[ip]) for ip in vs]
        xc, yc = flopy.utils.cvfdutil.centroid_of_polygon(points)
        cell2d.append([icell2d, xc, yc, len(vs), *vs])
    return vertices, cell2d


def make_grid():
    # itri 0 means do not split the cell
    # itri 1 means split from upper right to lower left
    # itri 2 means split from upper left to lower right
    itri = np.zeros((nrow, ncol), dtype=np.int)
    itri[:, 1 : ncol - 1] = 2
    itri[source_location0[1], source_location0[2]] = 0
    delra = delr * np.ones(ncol, dtype=np.float)
    delca = delc * np.ones(nrow, dtype=np.float)
    verts, iverts = grid_triangulator(itri, delra, delca)
    vertices, cell2d = cvfd_to_cell2d(verts, iverts)

    # A grid array that has the cellnumber of the first triangular cell in
    # the original grid
    itricellnum = np.empty((nrow, ncol), dtype=np.int)
    icell = 0
    for i in range(nrow):
        for j in range(ncol):
            itricellnum[i, j] = icell
            if itri[i, j] != 0:
                icell += 2
            else:
                icell += 1
    return vertices, cell2d, itricellnum


def build_mf6gwf(sim_folder):
    print("Building mf6gwf model...{}".format(sim_folder))
    name = "flow"
    sim_ws = os.path.join(ws, sim_folder, "mf6gwf")
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=sim_ws, exe_name=config.mf6_exe
    )
    tdis_ds = ((total_time, 1, 1.0),)
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        inner_maximum=300,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    vertices, cell2d, itricellnum = make_grid()
    flopy.mf6.ModflowGwfdisv(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nvert=len(vertices),
        ncpl=len(cell2d),
        vertices=vertices,
        cell2d=cell2d,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=0,
        k=hydraulic_conductivity,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    chdspd = []
    welspd = []
    for k in range(nlay):
        for i in range(nrow):
            icpl = itricellnum[i, ncol - 1]
            rec = [(k, icpl), 0.0]
            chdspd.append(rec)
            icpl = itricellnum[i, 0]
            rec = [(k, icpl), specific_discharge * delc * delv]
            welspd.append(rec)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=welspd)
    head_filerecord = "{}.hds".format(name)
    budget_filerecord = "{}.bud".format(name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


# MODFLOW 6 flopy GWF simulation object (sim) is returned


def build_mf6gwt(sim_folder):
    print("Building mf6gwt model...{}".format(sim_folder))
    name = "trans"
    sim_ws = os.path.join(ws, sim_folder, "mf6gwt")
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=sim_ws, exe_name=config.mf6_exe
    )
    tdis_ds = ((total_time, 100, 1.0),)
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=0.01,
        inner_dvclose=0.01,
        under_relaxation="simple",
        under_relaxation_gamma=0.9,
        relaxation_factor=0.99,
        linear_acceleration="bicgstab",
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name, save_flows=True)
    vertices, cell2d, itricellnum = make_grid()
    flopy.mf6.ModflowGwfdisv(
        gwt,
        length_units=length_units,
        nlay=nlay,
        nvert=len(vertices),
        ncpl=len(cell2d),
        vertices=vertices,
        cell2d=cell2d,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtdsp(
        gwt, alh=alpha_l, ath1=alpha_th, ath2=alpha_tv,
    )
    pd = [
        ("GWFHEAD", "../mf6gwf/flow.hds".format(), None),
        ("GWFBUDGET", "../mf6gwf/flow.bud", None),
    ]
    flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd)
    sourcerecarray = [[]]
    icpl = itricellnum[source_location0[1], source_location0[2]]
    srcspd = [[(0, icpl), solute_mass_flux]]
    flopy.mf6.ModflowGwtsrc(gwt, stress_period_data=srcspd)
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    obs_data = {
        "{}.obs.csv".format(name): [
            ("SOURCELOC", "CONCENTRATION", source_location0),
        ],
    }
    obs_package = flopy.mf6.ModflowUtlobs(
        gwt, digits=10, print_input=True, continuous=obs_data
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord="{}.cbc".format(name),
        concentration_filerecord="{}.ucn".format(name),
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )
    return sim


def build_model(sim_name):
    sims = None
    if config.buildModel:
        sim_mf6gwf = build_mf6gwf(sim_name)
        sim_mf6gwt = build_mf6gwt(sim_name)
        sims = (sim_mf6gwf, sim_mf6gwt)
    return sims


# Function to write model files


def write_model(sims, silent=True):
    if config.writeModel:
        sim_mf6gwf, sim_mf6gwt = sims
        sim_mf6gwf.write_simulation(silent=silent)
        sim_mf6gwt.write_simulation(silent=silent)
    return


# Function to run the model
# True is returned if the model runs successfully


def run_model(sims, silent=True):
    success = True
    if config.runModel:
        success = False
        sim_mf6gwf, sim_mf6gwt = sims
        success, buff = sim_mf6gwf.run_simulation(silent=silent)
        if not success:
            print(buff)
        success, buff = sim_mf6gwt.run_simulation(silent=silent)
        if not success:
            print(buff)
    return success


# Function to plot the model results


def plot_analytical(ax, levels):
    n = porosity
    v = velocity_x
    al = alpha_l
    ath = alpha_th
    atv = alpha_tv
    c0 = 10.0
    xc = [22.5]
    yc = [0]
    zc = [0]
    q = [1.0]
    dx = v * al
    dy = v * ath
    dz = v * atv
    lam = 0.0
    x = np.arange(0 + delr / 2.0, ncol * delr + delr / 2.0, delr)
    y = np.arange(0 + delc / 2.0, nrow * delc + delc / 2.0, delc)
    x, y = np.meshgrid(x, y)
    z = 0
    t = 400.0
    c400 = analytical.Wexler3d().multiwell(
        x, y, z, t, v, xc, yc, zc, dx, dy, dz, n, q, lam, c0
    )
    cs = ax.contour(x, y, c400, levels=levels, colors="k")
    return cs


def plot_grid(sims):
    if config.plotModel:
        print("Plotting model results...")
        sim_mf6gwf, sim_mf6gwt = sims
        fs = USGSFigure(figure_type="map", verbose=False)

        sim_ws = sim_mf6gwt.simulation_data.mfpath.get_sim_path()
        fig, axs = plt.subplots(
            1, 1, figsize=figure_size, dpi=300, tight_layout=True
        )
        gwt = sim_mf6gwt.trans
        pmv = flopy.plot.PlotMapView(model=gwt, ax=axs)
        pmv.plot_grid()
        axs.set_xlabel("x position (m)")
        axs.set_ylabel("y position (m)")
        axs.set_aspect(4.0)

        # save figure
        if config.plotSave:
            sim_folder = os.path.split(sim_ws)[0]
            sim_folder = os.path.basename(sim_folder)
            fname = "{}-grid{}".format(sim_folder, config.figure_ext)
            fpth = os.path.join(ws, "..", "figures", fname)
            fig.savefig(fpth)


def plot_results(sims):
    if config.plotModel:
        print("Plotting model results...")
        sim_mf6gwf, sim_mf6gwt = sims
        fs = USGSFigure(figure_type="map", verbose=False)

        sim_ws = sim_mf6gwt.simulation_data.mfpath.get_sim_path()
        fname = os.path.join(sim_ws, "trans.ucn")
        conc = flopy.utils.HeadFile(
            fname, text="concentration", precision="double"
        )
        conc = conc.get_data()

        fig, axs = plt.subplots(
            1, 1, figsize=figure_size, dpi=300, tight_layout=True
        )

        gwt = sim_mf6gwt.trans
        pmv = flopy.plot.PlotMapView(model=gwt, ax=axs)
        # pmv.plot_array(conc, alpha=0.5)
        # pmv.plot_grid()
        levels = [1, 3, 10, 30, 100, 300]
        cs1 = plot_analytical(axs, levels)
        cs2 = pmv.contour_array(
            conc, colors="blue", linestyles="--", levels=levels
        )
        axs.set_xlabel("x position (m)")
        axs.set_ylabel("y position (m)")
        axs.set_aspect(4.0)

        labels = ["Analytical", "MODFLOW 6"]
        lines = [cs1.collections[0], cs2.collections[0]]
        axs.legend(lines, labels, loc="upper left")

        # save figure
        if config.plotSave:
            sim_folder = os.path.split(sim_ws)[0]
            sim_folder = os.path.basename(sim_folder)
            fname = "{}-map{}".format(sim_folder, config.figure_ext)
            fpth = os.path.join(ws, "..", "figures", fname)
            fig.savefig(fpth)


# Function that wraps all of the steps for each scenario
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def scenario(idx, silent=True):
    sims = build_model(example_name)
    write_model(sims, silent=silent)
    success = run_model(sims, silent=silent)
    if success:
        plot_grid(sims)
        plot_results(sims)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Model

    # Model run

    scenario(0)
