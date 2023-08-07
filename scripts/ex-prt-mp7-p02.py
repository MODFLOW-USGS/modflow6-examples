# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# ## Backward particle tracking through a steady-state flow field on an unstructured (DISV quadtree) grid
#
# Application of a MODFLOW 6 particle-tracking (PRT) model to solve example 2 from the MODPATH 7 documentation. This example problem demonstrates a steady-state MODFLOW 6 simulation using a quadpatch DISV grid.
#
# In part A, 16 particles are evenly distributed around the 4 side faces of a cell containing a well in a locally-refined region in the center of the grid, then tracked backwards to recharge locations at the water table.
#
# In part B, 100 particles are evenly distributed around the 4 side faces of the well's cell, with 16 additional particles distributed over the cell's top face, then particles are again tracked backwards to recharge locations at the water table.
#
# ### Problem setup
#
# First import dependencies.

# +
import os
import sys
from warnings import warn
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
from shapely.geometry import MultiPoint, LineString
import flopy
from flopy.utils.gridintersect import GridIntersect
import numpy as np
import pandas as pd


try:
    # append the common/ subdirectory to the system path
    # (assumes running one level down from project root)
    sys.path.append(os.path.join("..", "common"))
    import config

    buildModel = config.buildModel
    writeModel = config.writeModel
    runModel = config.runModel
    plotModel = config.plotModel
    plotSave = config.plotSave
    figure_ext = config.figure_ext
    base_ws = config.base_ws
    timeit = config.timeit
except:
    warn(f"Failed to import config")
    # default settings
    buildModel = True
    writeModel = True
    runModel = True
    plotModel = True
    plotSave = False
    figure_ext = ".png"
    base_ws = Path("../examples")

    def timeit(func):
        return func


try:
    from figspecs import USGSFigure

    styles = True
except:
    warn(f"Failed to import figure styles")
    styles = False
# -

# Define simulation and model name and workspace.

# +
sim_name = "mp7-p02"
example_name = "ex-prt-" + sim_name
gwf_name = sim_name + "-gwf"
prt_name = sim_name + "-prt"
mp7_name = sim_name + "-mp7"

sim_ws = base_ws / example_name
gwf_ws = sim_ws / "gwf"
prt_ws = sim_ws / "prt"
mp7_ws = sim_ws / "mp7"

# We only need 1 working directory for the GWF model
# since it can be reused for both subproblems A & B.
# The PRT and MP7 models each need separate working
# directories for each subproblem.
gwf_ws.mkdir(exist_ok=True, parents=True)

headfile = f"{gwf_name}.hds"
headfile_bkwd = f"{gwf_name}_bkwd.hds"
budgetfile = f"{gwf_name}.cbb"
budgetfile_bkwd = f"{gwf_name}_bkwd.cbb"
trackfile = f"{prt_name}.trk"
trackhdrfile = f"{prt_name}.trk.hdr"
trackcsvfile = f"{prt_name}.trk.csv"

# -

# Define model units.

length_units = "feet"
time_units = "days"

# Define time discretization parameters.

tdis_rc = [(1000.0, 1, 1.0)]
nper = len(tdis_rc)

# Define MODFLOW 6 flow model parameters.

# +
# Cell types by layer
icelltype = [1, 0, 0]

# Conductivities
k = [50.0, 0.01, 200.0]
k33 = [10.0, 0.01, 20.0]

# Well
# issue(flopyex): in the original flopy example, the well
#   coordinates were at the edge (not the center) of the cell,
#   but it apparently worked out anyway
wel_coords = [(4718.45, 5281.25)]
wel_q = [-150000.0]

# Recharge
rch = 0.005
rch_iface = 6
rch_iflowface = -1

# River
riv_h = 320.0
riv_z = 318.0
riv_c = 1.0e5
riv_iface = 6
riv_iflowface = -1
# -

# Initialize globally available lists of well and river cells, which will be filled when the MODFLOW 6 flow model is built.

welcells = []
rivcells = []

# ### Grid refinement
#
# [GRIDGEN](https://www.usgs.gov/software/gridgen-program-generating-unstructured-finite-volume-grids) can be used to create a quadpatch grid with a central refined region.
#
# The grid will have 3 refinement levels. First, create the top-level (base) grid discretization.

# +
Lx = 10000.0
Ly = 10500.0
nlay = 3
nrow = 21
ncol = 20
delr = Lx / ncol
delc = Ly / nrow
top = 400
botm = [220, 200, 0]

ms = flopy.modflow.Modflow()
dis = flopy.modflow.ModflowDis(
    ms,
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=delr,
    delc=delc,
    top=top,
    botm=botm,
)
# -

# Refine the grid.

# +
from flopy.utils.gridgen import Gridgen

# create Gridgen workspace
gridgen_ws = sim_ws / "gridgen"
gridgen_ws.mkdir(parents=True, exist_ok=True)

# create Gridgen object
g = Gridgen(ms.modelgrid, model_ws=gridgen_ws)

# add polygon for each refinement level
outer_polygon = [[(3500, 4000), (3500, 6500), (6000, 6500), (6000, 4000), (3500, 4000)]]
g.add_refinement_features([outer_polygon], "polygon", 1, range(nlay))
refshp0 = gridgen_ws / "rf0"

middle_polygon = [
    [(4000, 4500), (4000, 6000), (5500, 6000), (5500, 4500), (4000, 4500)]
]
g.add_refinement_features([middle_polygon], "polygon", 2, range(nlay))
refshp1 = gridgen_ws / "rf1"

inner_polygon = [[(4500, 5000), (4500, 5500), (5000, 5500), (5000, 5000), (4500, 5000)]]
g.add_refinement_features([inner_polygon], "polygon", 3, range(nlay))
refshp2 = gridgen_ws / "rf2"
# -

# Build the grid and plot it with refinement levels superimposed.

# +
g.build(verbose=False)
grid_props = g.get_gridprops_vertexgrid()
disv_props = g.get_gridprops_disv()
grid = flopy.discretization.VertexGrid(**grid_props)

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
pmv = flopy.plot.PlotMapView(model=ms)
grid.plot(ax=ax)

flopy.plot.plot_shapefile(refshp0, ax=ax, facecolor="green", alpha=0.3)
flopy.plot.plot_shapefile(refshp1, ax=ax, facecolor="green", alpha=0.5)
flopy.plot.plot_shapefile(str(refshp2), ax=ax, facecolor="green", alpha=0.7)
# -

# ### Model creation
#
# We are now ready to setup groundwater flow and particle tracking models.
#
# Define shared MODFLOW 6 PRT and MODPATH 7 particle-tracking model parameters.

# Porosity
porosity = 0.1

# Define particle release points for the PRT model. Note that release points for part A and B are not identical, with more particles released in part B.

# +
releasepts = []

# retrieve GRIDGEN-generated gridprops
ncpl = disv_props["ncpl"]
top = disv_props["top"]
botm = disv_props["botm"]
nvert = disv_props["nvert"]
vertices = disv_props["vertices"]
cell2d = disv_props["cell2d"]


def set_releasepts(part):
    from math import sqrt

    global releasepts

    welcellnode = welcells[0]
    xctr = cell2d[welcellnode][1]
    yctr = cell2d[welcellnode][2]
    vert0 = cell2d[welcellnode][4]
    vert1 = cell2d[welcellnode][5]
    x0 = vertices[vert0][1]
    y0 = vertices[vert0][2]
    x1 = vertices[vert1][1]
    y1 = vertices[vert1][2]
    dx = x1 - x0
    dy = y1 - y0
    delcell = sqrt(dx * dx + dy * dy)
    delcellhalf = 0.5 * delcell

    # Reinitialize for the current scenario
    releasepts = []

    if part == "a":
        # Example 2A

        zrpt = 0.5 * (botm[1][welcellnode] + botm[2][welcellnode])
        npart_on_side = 4
        delta = delcell / npart_on_side
        baseoffset = 0.5 * (npart_on_side + 1) * delta
        xbase = xctr - baseoffset
        ybase = yctr - baseoffset
        nrpt = -1
        for idiv in range(npart_on_side):
            i = idiv + 1
            xrpt = xbase + i * delta
            yrpt = yctr + delcellhalf
            nrpt += 1
            rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
            releasepts.append(rpt)
            yrpt = yctr - delcellhalf
            nrpt += 1
            rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
            releasepts.append(rpt)
            yrpt = ybase + i * delta
            xrpt = xctr + delcellhalf
            nrpt += 1
            rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
            releasepts.append(rpt)
            xrpt = xctr - delcellhalf
            nrpt += 1
            rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
            releasepts.append(rpt)

    else:
        # Example 2B

        # 4x4 array of particles on top of well cell
        zrpt = botm[1][welcellnode]
        npart_on_side = 4
        delta = delcell / npart_on_side
        baseoffset = 0.5 * (npart_on_side + 1) * delta
        xbase = xctr - baseoffset
        ybase = yctr - baseoffset
        nrpt = -1
        for idivx in range(npart_on_side):
            ix = idivx + 1
            xrpt = xbase + ix * delta
            for idivy in range(npart_on_side):
                iy = idivy + 1
                yrpt = ybase + iy * delta
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
        # 10x10 arrays of particles on the four sides of well cell
        zwel = 0.5 * (botm[1][welcellnode] + botm[2][welcellnode])
        npart_on_side = 10
        delcellz = botm[1][welcellnode] - botm[2][welcellnode]
        delta = delcell / npart_on_side
        deltaz = delcellz / npart_on_side
        baseoffset = 0.5 * (npart_on_side + 1) * delta
        baseoffsetz = 0.5 * (npart_on_side + 1) * deltaz
        xbase = xctr - baseoffset
        ybase = yctr - baseoffset
        zbase = zwel - baseoffsetz
        for idivz in range(npart_on_side):
            iz = idivz + 1
            zrpt = zbase + iz * deltaz
            for idiv in range(npart_on_side):
                i = idiv + 1
                xrpt = xbase + i * delta
                yrpt = yctr + delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
                yrpt = yctr - delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
                yrpt = ybase + i * delta
                xrpt = xctr + delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
                xrpt = xctr - delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)


# -

# Define particle release points for the MODPATH 7 model.


def get_particle_data(part):
    nodew = ncpl * 2 + welcells[0]

    if part == "a":
        # Example 2A
        pcoord = np.array(
            [
                [0.000, 0.125, 0.500],
                [0.000, 0.375, 0.500],
                [0.000, 0.625, 0.500],
                [0.000, 0.875, 0.500],
                [1.000, 0.125, 0.500],
                [1.000, 0.375, 0.500],
                [1.000, 0.625, 0.500],
                [1.000, 0.875, 0.500],
                [0.125, 0.000, 0.500],
                [0.375, 0.000, 0.500],
                [0.625, 0.000, 0.500],
                [0.875, 0.000, 0.500],
                [0.125, 1.000, 0.500],
                [0.375, 1.000, 0.500],
                [0.625, 1.000, 0.500],
                [0.875, 1.000, 0.500],
            ]
        )
        plocs = [nodew for i in range(pcoord.shape[0])]
        pgdata = flopy.modpath.ParticleData(
            plocs,
            structured=False,
            localx=pcoord[:, 0],
            localy=pcoord[:, 1],
            localz=pcoord[:, 2],
            drape=0,
        )

    else:
        # Example 2B
        facedata = flopy.modpath.FaceDataType(
            drape=0,
            verticaldivisions1=10,
            horizontaldivisions1=10,
            verticaldivisions2=10,
            horizontaldivisions2=10,
            verticaldivisions3=10,
            horizontaldivisions3=10,
            verticaldivisions4=10,
            horizontaldivisions4=10,
            rowdivisions5=0,
            columndivisions5=0,
            rowdivisions6=4,
            columndivisions6=4,
        )
        pgdata = flopy.modpath.NodeParticleData(subdivisiondata=facedata, nodes=nodew)

    return pgdata


# Define some variables used to plot model results.

# +
# colormap for boundary locations
cmapbd = mpl.colors.ListedColormap(
    [
        "r",
        "g",
    ]
)

# time series point colors by layer
colors = ["green", "orange", "red"]

# figure sizes
figure_size_solo = (8.0, 8.0)
figure_size_compare = (15, 6)


# -

# Define functions to build, write, run, and plot models.
#
# Below we use three distinct simulations:
#
#   1. MODFLOW 6 GWF (groundwater flow)
#   2. MODFLOW 6 PRT (particle-tracking)
#   3. MODPATH 7 particle-tracking


def build_mf6gwf():
    global welcells, rivcells

    print("Building GWF model")

    # ====================================
    # Create the MODFLOW 6 flow simulation
    # ====================================

    # Instantiate the MODFLOW 6 simulation object
    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name, exe_name="mf6", version="mf6", sim_ws=gwf_ws
    )

    # Instantiate the MODFLOW 6 temporal discretization package
    flopy.mf6.ModflowTdis(
        sim, pname="tdis", time_units="DAYS", perioddata=tdis_rc, nper=len(tdis_rc)
    )

    # -------------------
    # Build the GWF model
    # -------------------

    # Instantiate the MODFLOW 6 gwf (groundwater-flow) model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwf_name, model_nam_file="{}.nam".format(gwf_name)
    )
    gwf.name_file.save_flows = True

    # Instantiate the MODFLOW 6 gwf discretization package
    flopy.mf6.ModflowGwfdisv(
        gwf,
        length_units=length_units,
        **disv_props,
    )

    # GridIntersect object for setting up boundary conditions
    ix = GridIntersect(gwf.modelgrid, method="vertex", rtree=True)

    # Instantiate the MODFLOW 6 gwf initial conditions package
    flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=riv_h)

    # Instantiate the MODFLOW 6 gwf node property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=[("xt3d")],
        icelltype=icelltype,
        k=k,
        k33=k33,
        save_saturation=True,
        save_specific_discharge=True,
    )

    # Instantiate the MODFLOW 6 gwf recharge package
    flopy.mf6.ModflowGwfrcha(
        gwf,
        recharge=rch,
        auxiliary=["iface", "iflowface"],
        aux=[rch_iface, rch_iflowface],
    )

    # Instantiate the MODFLOW 6 gwf well package
    welcells = ix.intersects(MultiPoint(wel_coords))
    welcells = [icpl for (icpl,) in welcells]
    welspd = [[(2, icpl), wel_q[idx]] for idx, icpl in enumerate(welcells)]
    flopy.mf6.ModflowGwfwel(gwf, print_input=True, stress_period_data=welspd)

    # Instantiate the MODFLOW 6 gwf river package
    riverline = [(Lx - 1.0, Ly), (Lx - 1.0, 0.0)]
    rivcells = ix.intersects(LineString(riverline))
    rivcells = [icpl for (icpl,) in rivcells]
    rivspd = [
        [(0, icpl), riv_h, riv_c, riv_z, riv_iface, riv_iflowface] for icpl in rivcells
    ]
    flopy.mf6.ModflowGwfriv(
        gwf, stress_period_data=rivspd, auxiliary=[("iface", "iflowface")]
    )

    # Instantiate the MODFLOW 6 gwf output control package
    headfile = "{}.hds".format(gwf_name)
    head_record = [headfile]
    budgetfile = "{}.cbb".format(gwf_name)
    budget_record = [budgetfile]
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="oc",
        budget_filerecord=budget_record,
        head_filerecord=head_record,
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Create an iterative model solution (IMS) for the MODFLOW 6 gwf model
    ims = flopy.mf6.ModflowIms(
        sim,
        pname="ims",
        print_option="SUMMARY",
        complexity="SIMPLE",
        outer_dvclose=1.0e-5,
        outer_maximum=100,
        under_relaxation="NONE",
        inner_maximum=100,
        inner_dvclose=1.0e-6,
        rcloserecord=0.1,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=0.99,
    )
    sim.register_ims_package(ims, [gwf.name])

    return sim


def build_mf6prt(part):
    print(f"Building PRT model for {example_name}{part}")

    prpname = f"prp2{part}"
    prpfilename = f"{prt_name}_2{part}.prp"

    prt_ws_scen = prt_ws.parent / (prt_ws.name + part)
    prt_ws_scen.mkdir(exist_ok=True)

    # Instantiate the MODFLOW 6 simulation object
    simprt = flopy.mf6.MFSimulation(
        sim_name=prt_name, version="mf6", exe_name="mf6", sim_ws=prt_ws_scen
    )

    # Create a time discretization for backward tracking;
    # since there is only a single period with a single time step,
    # the time discretization is the same backward as forward
    tdis_bkwd = tdis_rc

    # Instantiate the MODFLOW 6 temporal discretization package
    flopy.mf6.ModflowTdis(
        simprt,
        pname="tdis",
        time_units="DAYS",
        perioddata=tdis_bkwd,
        nper=len(tdis_bkwd),
    )

    # -------------------
    # Build the PRT model
    # -------------------

    # Instantiate the MODFLOW 6 prt model
    prt = flopy.mf6.ModflowPrt(
        simprt, modelname=prt_name, model_nam_file="{}.nam".format(prt_name)
    )

    # Instantiate the MODFLOW 6 prt discretization package
    flopy.mf6.ModflowGwfdisv(
        prt,
        length_units=length_units,
        **disv_props,
    )

    # Instantiate the MODFLOW 6 prt model input package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # Set particle release point data according to the scenario
    set_releasepts(part)
    nreleasepts = len(releasepts)
    particle_data = get_particle_data(part)
    prp_pkg_data = list(particle_data.to_prp(grid))

    # Instantiate the MODFLOW 6 prt particle release point (prp) package
    pd = {0: ["FIRST"], 1: []}
    flopy.mf6.ModflowPrtprp(
        prt,
        pname=prpname,
        filename=prpfilename,
        nreleasepts=len(prp_pkg_data),
        packagedata=prp_pkg_data,
        perioddata=pd,
    )

    # Instantiate the MODFLOW 6 prt output control package
    budgetfile = "{}.bud".format(prt_name)
    trackfile = "{}.trk".format(prt_name)
    trackcsvfile = "{}.trk.csv".format(prt_name)
    budget_record = [budgetfile]
    track_record = [trackfile]
    trackcsv_record = [trackcsvfile]
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=budget_record,
        track_filerecord=track_record,
        trackcsv_filerecord=trackcsv_record,
        saverecord=[("BUDGET", "ALL")],
    )

    # Instantiate the MODFLOW 6 prt flow model interface
    # using "time-reversed" budget and head files
    pd = [
        ("GWFHEAD", headfile_bkwd),
        ("GWFBUDGET", budgetfile_bkwd),
    ]
    flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)

    # Create an explicit model solution (EMS) for the MODFLOW 6 prt model
    ems = flopy.mf6.ModflowEms(
        simprt,
        pname="ems",
        filename="{}.ems".format(prt_name),
    )
    simprt.register_solution_package(ems, [prt.name])

    return simprt


def build_mp7(part, gwf):
    print(f"Building mp7 model for {example_name}{part}")

    # Set parameters according to the scenario
    pgdata = get_particle_data(part)
    nm_mp7 = f"{mp7_name}{part}"
    pgname = f"BACKWARD2{part.upper()}"
    fpth = nm_mp7 + ".sloc"

    if part == "a":
        simtype = "combined"
        timepointdata = [500, 1000.0]
        pg = flopy.modpath.ParticleGroup(
            particlegroupname=pgname, particledata=pgdata, filename=fpth
        )
    else:
        pg = flopy.modpath.ParticleGroupNodeTemplate(
            particlegroupname=pgname, particledata=pgdata, filename=fpth
        )
        simtype = "combined"
        timepointdata = None

    # Instantiate the MODPATH 7 simulation object
    bfile = gwf_ws / budgetfile
    hfile = gwf_ws / headfile

    # create workspace for this scenario (part A or B)
    mp7_ws_scen = mp7_ws.parent / (mp7_ws.name + part)
    mp7_ws_scen.mkdir(exist_ok=True)

    mp7 = flopy.modpath.Modpath7(
        modelname=nm_mp7,
        flowmodel=gwf,
        exe_name="mp7",
        model_ws=mp7_ws_scen,
        budgetfilename=budgetfile,
        headfilename=headfile,
    )

    # Instantiate the MODPATH 7 basic data
    flopy.modpath.Modpath7Bas(mp7, porosity=porosity)

    # Instantiate the MODPATH 7 simulation data
    flopy.modpath.Modpath7Sim(
        mp7,
        simulationtype=simtype,
        trackingdirection="backward",
        weaksinkoption="pass_through",
        weaksourceoption="pass_through",
        referencetime=0.0,
        stoptimeoption="extend",
        timepointdata=timepointdata,
        particlegroups=pg,
    )

    # Return the simulation
    return mp7


def build_models(part):
    sim = build_mf6gwf()
    gwf = sim.get_model(gwf_name)
    simprt = build_mf6prt(part)
    mp7 = build_mp7(part, gwf)
    return sim, simprt, mp7


# Function to write simulation/model input files.


def write_models(sim, simprt, mp7, silent=True):
    sim.write_simulation(silent=silent)
    simprt.write_simulation(silent=silent)
    mp7.write_input()


# Because this problem tracks particles backwards, we need to reverse the head and budget files after running the groundwater flow model and before running the particle tracking model. Define functions to do this.

# +
import flopy.utils.binaryfile as bf


def reverse_budgetfile(fpth, rev_fpth, tdis):
    f = bf.CellBudgetFile(fpth, tdis=tdis)
    f.reverse(rev_fpth)


def reverse_headfile(fpth, rev_fpth, tdis):
    f = bf.HeadFile(fpth, tdis=tdis)
    f.reverse(rev_fpth)


# -

# Define a function to run all three simulations and their respective models.


@timeit
def run_models(part, sim, simprt, mp7, silent=True):
    # Run MODFLOW 6 flow simulation
    success, buff = sim.run_simulation(silent=silent, report=True)
    for line in buff:
        print(line)
    assert success

    # Process budget and head files for backward tracking
    prt_ws_scen = prt_ws.parent / (prt_ws.name + part)
    reverse_budgetfile(gwf_ws / budgetfile, prt_ws_scen / budgetfile_bkwd, sim.tdis)
    reverse_headfile(gwf_ws / headfile, prt_ws_scen / headfile_bkwd, sim.tdis)

    # Run MODFLOW 6 PRT particle-tracking simulation
    success, buff = simprt.run_simulation(silent=silent, report=True)
    for line in buff:
        print(line)
    assert success

    # Run MODPATH 7 simulation
    success, buff = mp7.run_model(silent=silent, report=True)
    for line in buff:
        print(line)
    assert success


# Define functions to load pathline and endpoint data from MODFLOW 6 PRT's particle track CSV files. Note that unlike MODPATH 7, MODFLOW 6 PRT does not make a distinction between pathline and endpoint output files &mdash; all pathline data is saved to track files, and endpoints are computed dynamically.

# +
from flopy.plot.plotutil import to_mp7_pathlines, to_mp7_endpoints


def load_mf6pathlines(p):
    pls = pd.read_csv(p)
    pls = to_mp7_pathlines(pls)
    return pls


def load_mf6endpoints(p):
    pls = pd.read_csv(p)
    eps = pls.sort_values("t").groupby(["imdl", "iprp", "irpt", "trelease"]).tail(1)
    return eps


# -

# Define a function for plotting flow model grid, boundary conditions, and head results.


# +
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


def plot_nodes_and_vertices(gwf, mg, ibd, ax):
    """
    Plot cell nodes and vertices (and IDs) on a zoomed inset
    """

    ax.set_aspect("equal")

    # set zoom area
    xmin, xmax = 2050, 4800
    ymin, ymax = 5200, 7550
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    # create map view plot
    pmv = flopy.plot.PlotMapView(gwf, ax=ax)
    v = pmv.plot_grid(lw=0.5, edgecolor="black")
    t = ax.set_title("Node and vertex indices (one-based)\n", fontsize=14)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    # plot vertices
    verts = mg.verts
    ax.plot(verts[:, 0], verts[:, 1], "bo")
    for i in range(ncpl):
        x, y = verts[i, 0], verts[i, 1]
        if xmin <= x <= xmax and ymin <= y <= ymax:
            ax.annotate(str(i + 1), verts[i, :], color="b")

    # plot nodes
    xc, yc = mg.get_xcellcenters_for_layer(0), mg.get_ycellcenters_for_layer(0)
    for i in range(ncpl):
        x, y = xc[i], yc[i]
        ax.plot(x, y, "o", color="grey")
        if xmin <= x <= xmax and ymin <= y <= ymax:
            ax.annotate(str(i + 1), (x, y), color="grey")

    # plot well
    ax.plot(wel_coords[0][0], wel_coords[0][1], "ro")

    # create legend
    ax.legend(
        handles=[
            mpl.patches.Patch(color="blue", label="vertex"),
            mpl.patches.Patch(color="grey", label="node"),
        ],
        loc="upper left",
    )


def plot_bc(sim, mg, ibd):
    if styles:
        fs = USGSFigure(figure_type="map", verbose=False)
    fig = plt.figure(figsize=figure_size_solo)
    fig.tight_layout()
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    gwf = sim.get_model(gwf_name)
    pmv = flopy.plot.PlotMapView(gwf, ax=ax)
    pmv.plot_bc("WEL", plotAll=True)
    pmv.plot_grid(lw=0.5)
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    pc = pmv.plot_array(ibd, cmap=cmapbd, edgecolor="gray")
    t = ax.set_title("Boundary conditions\n", fontsize=14)

    # create inset
    axins = ax.inset_axes([-0.76, 0.25, 0.7, 0.9])
    plot_nodes_and_vertices(gwf, mg, ibd, axins)
    ax.indicate_inset_zoom(axins)

    # create legend
    ax.legend(
        handles=[
            mpl.patches.Patch(color="red", label="Well"),
            mpl.patches.Patch(color="green", label="River"),
        ],
        loc="upper left",
    )

    # set title
    plt.suptitle(t="Example 2: model grid", fontsize=14, ha="center")

    plt.show()


def plot_head(mg):
    # Import simulated head values from the binary head file
    fname = gwf_ws / (gwf_name + ".hds")
    hdobj = flopy.utils.HeadFile(fname)
    head = hdobj.get_data()

    # Prepare a plot of heads in layer 3
    ilay = 2
    cint = 0.25
    if styles:
        fs = USGSFigure(figure_type="map", verbose=False)
    fig2 = plt.figure(figsize=figure_size_solo)
    fig2.tight_layout()
    ax = fig2.add_subplot(1, 1, 1, aspect="equal")
    mm = flopy.plot.PlotMapView(modelgrid=mg, ax=ax, layer=ilay)
    mm.plot_grid(lw=0.5)
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    pc = mm.plot_array(head[:, 0, :], cmap="jet", edgecolor="black")
    hmin = head[ilay, 0, :].min()
    hmax = head[ilay, 0, :].max()
    levels = np.arange(np.floor(hmin), np.ceil(hmax) + cint, cint)
    cs = mm.contour_array(head[:, 0, :], colors="white", levels=levels)
    plt.clabel(cs, fmt="%.1f", colors="white", fontsize=11)
    cb = plt.colorbar(pc, shrink=0.5)
    t = ax.set_title(
        "Example 2: Head in layer {}, hmin {:6.2f}, hmax {:6.2f}\n".format(
            ilay + 1, hmin, hmax
        ),
        fontsize=14,
    )
    plt.show()


def plot_gwf(sim, mg, ibd):
    plot_bc(sim, mg, ibd)
    plot_head(mg)


# -


# Define a function for plotting pathlines and time series points colored by layer.


def plot_pathlines_and_timeseries(ax, mg, ibd, pathlines, timeseries, plottitle):
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(modelgrid=mg, ax=ax)
    mm.plot_grid(lw=0.5)
    v = mm.plot_array(ibd, cmap=cmapbd, edgecolor="gray")
    mm.plot_pathline(pathlines, layer="all", colors=["blue"], lw=0.75)
    # issue(mf6): no way to output time series in prt
    if timeseries != None:
        for k in range(nlay - 1, -1, -1):
            mm.plot_timeseries(timeseries, layer=k, marker="o", lw=0, color=colors[k])
    ax.set_title(plottitle, fontsize=12)


# Define a utility function for plotting particle endpoints, colored by travel time to capture.


def plot_endpoints(ax, mg, ibd, endpointdata, plottitle):
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    mm = flopy.plot.PlotMapView(modelgrid=mg, ax=ax)
    mm.plot_grid(lw=0.5)
    v = mm.plot_array(ibd, cmap=cmapbd, edgecolor="gray")
    mm.plot_endpoint(endpointdata, direction="ending", colorbar=True, shrink=0.25)
    ax.set_title(plottitle, fontsize=12)


def plot_results(part, sim, simprt, mp7):
    # Get model grid
    fname = gwf_ws / (gwf_name + ".disv.grb")
    grd = flopy.mf6.utils.MfGrdFile(fname, verbose=False)
    mg = grd.modelgrid

    # Workspaces
    prt_ws_scen = prt_ws.parent / (prt_ws.name + part)
    mp7_ws_scen = mp7_ws.parent / (mp7_ws.name + part)

    # identify the boundary locations
    ibd = np.zeros((ncpl), dtype=int)
    ibd[np.array(welcells)] = 1
    ibd[np.array(rivcells)] = 2
    ibd = np.ma.masked_equal(ibd, 0)

    if part == "a":
        # =========================
        # MODFLOW 6 flow simulation
        # =========================

        plot_gwf(sim, mg, ibd)

        # ==========
        # Example 2A
        # ==========

        # Load MODFLOW 6 PRT pathlines
        mf6pathlines = load_mf6pathlines(prt_ws_scen / trackcsvfile)

        # Load MODPATH 7 pathlines
        plf = flopy.utils.PathlineFile(mp7_ws_scen / (mp7_name + f"{part}.mppth"))
        mp7pathlines = plf.get_alldata()

        # ------------------------------------------------
        # Pathlines and time series point colored by layer
        # ------------------------------------------------

        # Initialize plot
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=figure_size_compare)
        plt.suptitle(
            t="Example 2A: Pathlines and time series points tracked backward from well",
            fontsize=14,
        )
        axes = axes.flatten()

        # MODFLOW 6 PRT
        ax = axes[0]
        # issue(mf6): no way to output time series
        plot_pathlines_and_timeseries(ax, mg, ibd, mf6pathlines, None, "MODFLOW 6 PRT")

        # MODPATH 7
        ax = axes[1]
        plot_pathlines_and_timeseries(
            ax, mg, ibd, mp7pathlines, mp7pathlines, "MODPATH 7"
        )

        # Save figure
        if plotSave:
            fpth = os.path.join(
                "..", "figures", "{}-paths{}".format(sim_name, figure_ext)
            )
            fig.savefig(fpth)

    else:
        # ==========
        # Example 2B
        # ==========

        # Load MODFLOW 6 PRT endpoint data
        mf6endpoints = load_mf6endpoints(prt_ws_scen / trackcsvfile)

        # Load MODPATH 7 endpoint data
        epf = flopy.utils.EndpointFile(mp7_ws_scen / (mp7_name + f"{part}.mpend"))
        mp7endpoints = epf.get_alldata()

        # -------------------------------------------
        # Endpoints colored by travel time to capture
        # -------------------------------------------

        # Initialize plot
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=figure_size_compare)
        plt.suptitle(
            t="Example 2B: Endpoints of particle paths tracked backward from well",
            fontsize=14,
        )
        axes = axes.flatten()

        # MODFLOW 6 PRT
        ax = axes[0]
        plot_endpoints(ax, mg, ibd, mf6endpoints, "MODFLOW 6 PRT")

        # MODPATH 7
        ax = axes[1]
        plot_endpoints(ax, mg, ibd, mp7endpoints, "MODPATH 7")

        # issue(flopy): the center "dot" looks to be in a slightly different location
        # in the two plots -- true even if you plot the very same endpoint array in both plots

        # Save figure
        if plotSave:
            fpth = os.path.join(
                "..", "figures", "{}-endpts{}".format(sim_name, figure_ext)
            )
            fig.savefig(fpth)


# Define a function to wrap all of the steps for each scenario.
#
# 1. build model
# 2. write model
# 3. run model
# 4. plot results


def scenario(part, silent=True):
    sim, simprt, mp7 = build_models(part)
    write_models(sim, simprt, mp7, silent=silent)
    if runModel:
        run_models(part, sim, simprt, mp7, silent=silent)
        plot_results(part, sim, simprt, mp7)


# Run the scenario for problem 2A.

scenario("a", silent=False)

# Scenario for example 2B
scenario("b")
