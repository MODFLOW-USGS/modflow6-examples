# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# ## Particle tracking: steady state, unstructured quad-refined grid
#
# Application of a MODFLOW 6 particle-tracking (PRT) model and a MODPATH 7 (MP7) model to solve example 2 from the MODPATH 7 documentation.
#
# This example problem adapts the flow system in example 1, consisting of two aquifers separated by a low conductivity confining layer, with an unstructured grid with a quad-refined region around the central well. A river still runs along the grid's right-hand boundary.
#
# In part A, 16 particles are distributed evenly for release around the four horizontal faces of the well. To determine recharge points, particles are then tracked backwards to the water table.
#
# In part B, 100 particles are evenly distributed in a 10 x 10 square over the horizontal faces of the well, with 16 more release points on the well cell's top face. Particles are again tracked backwards to determine the well's capture area.
#

# ### Initial setup
#
# Import dependencies, define the example name and workspace, and read settings from environment variables.

# +
import pathlib as pl
from pprint import pformat

import flopy
from flopy.mf6 import MFSimulation
from shapely.geometry import MultiPoint, LineString
from flopy.utils.gridintersect import GridIntersect
import git
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flopy.plot.styles import styles
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from modflow_devtools.misc import get_env, timed

# Example name and workspace paths. If this example is running
# in the git repository, use the folder structure described in
# the README. Otherwise just use the current working directory.
sim_name = "mp7-p02"
gwf_name = sim_name + "-gwf"
prt_name = sim_name + "-prt"
mp7_name = sim_name + "-mp7"
try:
    root = pl.Path(git.Repo(".", search_parent_directories=True).working_dir)
except:
    root = None
workspace = root / "examples" if root else pl.Path.cwd()
figs_path = root / "figures" if root else pl.Path.cwd()
sim_ws = workspace / sim_name
gwf_ws = sim_ws / "gwf"
prt_ws = sim_ws / "prt"
mp7_ws = sim_ws / "mp7"
gwf_ws.mkdir(exist_ok=True, parents=True)
prt_ws.mkdir(exist_ok=True, parents=True)
mp7_ws.mkdir(exist_ok=True, parents=True)

# Define output file names
headfile = f"{gwf_name}.hds"
budgetfile = f"{gwf_name}.cbb"
headfile_bkwd = f"{gwf_name}_bkwd.hds"
budgetfile_bkwd = f"{gwf_name}_bkwd.cbb"
budgetfile_prt = f"{prt_name}.cbb"
trackfile_prt = f"{prt_name}.trk"
trackhdrfile_prt = f"{prt_name}.trk.hdr"
trackcsvfile_prt = f"{prt_name}.trk.csv"

# Settings from environment variables
write = get_env("WRITE", True)
run = get_env("RUN", True)
plot = get_env("PLOT", True)
plot_show = get_env("PLOT_SHOW", True)
plot_save = get_env("PLOT_SAVE", True)
# -

# ### Define parameters
#
# Define model units, parameters and other settings.

# +
# Model units
length_units = "feet"
time_units = "days"

# Model parameters
nper = 1  # Number of periods
nlay = 3  # Number of layers
nrow = 21  # Number of rows
ncol = 20  # Number of columns
delr = 500.0  # Column width ($ft$)
delc = 500.0  # Row width ($ft$)
top = 400.0  # Top of the model ($m$)
botm_str = "220.0, 200.0, 0.0"  # Layer bottom elevations ($ft$)
porosity = 0.1  # Soil porosity ($\%$)
rch = 0.005  # Recharge rate ($ft/s$)
kh = [50.0, 0.01, 200.0]  # Horizontal hydraulic conductivity ($ft/s$)
kv = [10.0, 0.01, 20.0]  # Vertical hydraulic conductivity ($ft/s$)
wel_q = -150000.0  # Well pumping rate
riv_h = 320.0  # River stage ($ft$)
riv_z = 317.0  # River bottom ($ft$)
riv_c = 1.0e5  # River conductance ($l^2/t$)

# Time discretization
nstp = 1
perlen = 1000.0
tsmult = 1.0
tdis_rc = [(perlen, nstp, tsmult)]

# parse parameter strings into tuples
botm = [float(value) for value in botm_str.split(",")]

# Cell types by layer
icelltype = [1, 0, 0]

# Conductivities
k = [50.0, 0.01, 200.0]
k33 = [10.0, 0.01, 20.0]

# Well
wel_coords = [(4718.45, 5281.25)]
wel_q = [-150000.0]
welcells = []

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
rivcells = []

# -

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
# -

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

    if part == "A":
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

    if part == "A":
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


# -

# ### Model setup
#
# Define functions to build models, write input files, and run the simulation.


# +


def build_gwf_sim():
    global welcells, rivcells

    # Instantiate the MODFLOW 6 simulation object
    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name, exe_name="mf6", version="mf6", sim_ws=gwf_ws
    )

    # Instantiate the MODFLOW 6 temporal discretization package
    flopy.mf6.ModflowTdis(
        sim, pname="tdis", time_units="DAYS", perioddata=tdis_rc, nper=len(tdis_rc)
    )

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


def build_prt_sim():
    global welcells, rivcells

    # Instantiate the MODFLOW 6 simulation object
    sim = flopy.mf6.MFSimulation(
        sim_name=prt_name, exe_name="mf6", version="mf6", sim_ws=prt_ws
    )

    # Instantiate the MODFLOW 6 temporal discretization package
    flopy.mf6.ModflowTdis(
        sim, pname="tdis", time_units="DAYS", perioddata=tdis_rc, nper=len(tdis_rc)
    )

    # Instantiate the MODFLOW 6 prt model
    prt = flopy.mf6.ModflowPrt(
        sim, modelname=prt_name, model_nam_file="{}.nam".format(prt_name)
    )

    # Instantiate the MODFLOW 6 prt discretization package
    flopy.mf6.ModflowGwfdisv(
        prt,
        length_units=length_units,
        **disv_props,
    )

    # Instantiate the MODFLOW 6 prt model input package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    def add_prp(part):
        prpname = f"prp2{part}"
        prpfilename = f"{prt_name}_2{part}.prp"

        # Set particle release point data according to the scenario
        set_releasepts(part)
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

    add_prp("A")
    add_prp("B")

    # Instantiate the MODFLOW 6 prt output control package
    budgetfile = "{}.bud".format(prt_name)
    trackfile = "{}.trk".format(prt_name)
    trackcsvfile = "{}.trk.csv".format(prt_name)
    budget_record = [budgetfile]
    track_record = [trackfile]
    trackcsv_record = [trackcsvfile]
    tracktimes = list(range(0, 72000, 1000))
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=budget_record,
        track_filerecord=track_record,
        trackcsv_filerecord=trackcsv_record,
        track_all=False,
        track_timesrecord=tracktimes,
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
        sim,
        pname="ems",
        filename="{}.ems".format(prt_name),
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_mp7_sim(gwf_model):
    # Create particle groups
    pga = flopy.modpath.ParticleGroup(
        particlegroupname="PG2A",
        particledata=get_particle_data("A"),
        filename=f"a.sloc",
    )
    pgb = flopy.modpath.ParticleGroupNodeTemplate(
        particlegroupname="PG2B",
        particledata=get_particle_data("B"),
        filename=f"b.sloc",
    )

    # Instantiate the MODPATH 7 simulation object
    mp7 = flopy.modpath.Modpath7(
        modelname=mp7_name,
        flowmodel=gwf_model,
        exe_name="mp7",
        model_ws=mp7_ws,
        budgetfilename=budgetfile,
        headfilename=headfile,
    )

    # Instantiate the MODPATH 7 basic data
    flopy.modpath.Modpath7Bas(mp7, porosity=porosity)

    # Instantiate the MODPATH 7 simulation data
    flopy.modpath.Modpath7Sim(
        mp7,
        simulationtype="combined",
        trackingdirection="backward",
        weaksinkoption="pass_through",
        weaksourceoption="pass_through",
        referencetime=0.0,
        stoptimeoption="extend",
        timepointdata=[500, 1000.0],
        particlegroups=[pga, pgb],
    )

    # Return the simulation
    return mp7


def build_models():
    gwf_sim = build_gwf_sim()
    prt_sim = build_prt_sim()
    mp7_sim = build_mp7_sim(gwf_sim.get_model(gwf_name))
    return gwf_sim, prt_sim, mp7_sim


def write_models(*sims, silent=True):
    for sim in sims:
        if isinstance(sim, MFSimulation):
            sim.write_simulation(silent=silent)
        else:
            sim.write_input()


@timed
def run_models(*sims, silent=True):
    for sim in sims:
        if isinstance(sim, MFSimulation):
            success, buff = sim.run_simulation(silent=silent, report=True)
        else:
            sim.run_model(silent=silent, report=True)
        assert success, pformat(buff)

        if "gwf" in sim.name:
            # Reverse budget and head files for backward tracking
            reverse_budgetfile(gwf_ws / budgetfile, prt_ws / budgetfile_bkwd, sim.tdis)
            reverse_headfile(gwf_ws / headfile, prt_ws / headfile_bkwd, sim.tdis)


# -

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


# Define a function to load pathline data from MODFLOW 6 PRT and MODPATH 7 pathline files.

# +


def get_mf6_pathlines(path):
    # load mf6 pathlines
    mf6pl = pd.read_csv(path)

    # index by particle group and particle ID
    mf6pl.set_index(["iprp", "irpt"], drop=False, inplace=True)

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mf6pl["subprob"] = mf6pl.apply(lambda row: "A" if row.iprp == 1 else "B", axis=1)

    # add release time and termination time columns
    mf6pl["t0"] = (
        mf6pl.groupby(level=["iprp", "irpt"])
        .apply(lambda x: x.t.min())
        .to_frame(name="t0")
        .t0
    )
    mf6pl["tt"] = (
        mf6pl.groupby(level=["iprp", "irpt"])
        .apply(lambda x: x.t.max())
        .to_frame(name="tt")
        .tt
    )

    # add markercolor column, color-coding by layer for plots
    mf6pl["mc"] = mf6pl.apply(
        lambda row: "green" if row.ilay == 1 else "yellow" if row.ilay == 2 else "red",
        axis=1,
    )

    return mf6pl


def get_mp7_pathlines(path, gwf_model):
    # load mp7 pathlines, letting flopy determine capture areas
    mp7plf = flopy.utils.PathlineFile(path)
    mp7pl = pd.DataFrame(
        mp7plf.get_destination_pathline_data(
            list(range(gwf_model.modelgrid.nnodes)), to_recarray=True
        )
    )

    # index by particle group and particle ID
    mp7pl.set_index(["particlegroup", "sequencenumber"], drop=False, inplace=True)

    # convert indices to 1-based (flopy converts them to 0-based, but PRT uses 1-based, so do the same for consistency)
    kijnames = [
        "k",
        "node",
        "particleid",
        "particlegroup",
        "particleidloc",
        "sequencenumber",
    ]

    for n in kijnames:
        mp7pl[n] += 1

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mp7pl["subprob"] = mp7pl.apply(
        lambda row: (
            "A" if row.particlegroup == 1 else "B" if row.particlegroup == 2 else pd.NA
        ),
        axis=1,
    )

    # add release time and termination time columns
    mp7pl["t0"] = (
        mp7pl.groupby(level=["particlegroup", "sequencenumber"])
        .apply(lambda x: x.time.min())
        .to_frame(name="t0")
        .t0
    )
    mp7pl["tt"] = (
        mp7pl.groupby(level=["particlegroup", "sequencenumber"])
        .apply(lambda x: x.time.max())
        .to_frame(name="tt")
        .tt
    )

    # add markercolor column, color-coding by layer for plots
    mp7pl["mc"] = mp7pl.apply(
        lambda row: "green" if row.k == 1 else "yellow" if row.k == 2 else "red", axis=1
    )

    return mp7pl


def get_mp7_timeseries(path, gwf_model):
    # load mp7 pathlines, letting flopy determine capture areas
    mp7tsf = flopy.utils.TimeseriesFile(path)
    mp7ts = pd.DataFrame(
        mp7tsf.get_destination_timeseries_data(list(range(gwf_model.modelgrid.nnodes)))
    )

    # index by particle group and particle ID
    mp7ts.set_index(["particlegroup", "particleid"], drop=False, inplace=True)

    # convert indices to 1-based (flopy converts them to 0-based, but PRT uses 1-based, so do the same for consistency)
    kijnames = [
        "k",
        "node",
        "particleid",
        "particlegroup",
        "particleidloc",
    ]

    for n in kijnames:
        mp7ts[n] += 1

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mp7ts["subprob"] = mp7ts.apply(
        lambda row: (
            "A" if row.particlegroup == 1 else "B" if row.particlegroup == 2 else pd.NA
        ),
        axis=1,
    )

    # add release time and termination time columns
    mp7ts["t0"] = (
        mp7ts.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.min())
        .to_frame(name="t0")
        .t0
    )
    mp7ts["tt"] = (
        mp7ts.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.max())
        .to_frame(name="tt")
        .tt
    )

    # add markercolor column, color-coding by layer for plots
    mp7ts["mc"] = mp7ts.apply(
        lambda row: "green" if row.k == 1 else "yellow" if row.k == 2 else "red", axis=1
    )

    return mp7ts


def get_mp7_endpoints(path, gwf_model):
    # load mp7 pathlines, letting flopy determine capture areas
    mp7epf = flopy.utils.EndpointFile(path)
    mp7ep = pd.DataFrame(
        mp7epf.get_destination_endpoint_data(list(range(gwf_model.modelgrid.nnodes)))
    )

    # index by particle group and particle ID
    mp7ep.set_index(["particlegroup", "particleid"], drop=False, inplace=True)

    # convert indices to 1-based (flopy converts them to 0-based, but PRT uses 1-based, so do the same for consistency)
    kijnames = [
        "k",
        "node",
        "particleid",
        "particlegroup",
        "particleidloc",
    ]

    for n in kijnames:
        mp7ep[n] += 1

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mp7ep["subprob"] = mp7ep.apply(
        lambda row: (
            "A" if row.particlegroup == 1 else "B" if row.particlegroup == 2 else pd.NA
        ),
        axis=1,
    )

    # add release time and termination time columns
    mp7ep["t0"] = (
        mp7ep.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.min())
        .to_frame(name="t0")
        .t0
    )
    mp7ep["tt"] = (
        mp7ep.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.max())
        .to_frame(name="tt")
        .tt
    )

    # add markercolor column, color-coding by layer for plots
    mp7ep["mc"] = mp7ep.apply(
        lambda row: "green" if row.k == 1 else "yellow" if row.k == 2 else "red", axis=1
    )

    return mp7ep


# -

# ### Plotting results
#
# Define functions to plot model results.


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


def plot_nodes_and_vertices(gwf, ax):
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
    pmv.plot_grid(lw=0.5, edgecolor="black", alpha=0.5)
    styles.heading(ax=ax, heading="Nodes and vertices (indices one-based)")
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    # plot vertices
    mg = gwf.modelgrid
    verts = mg.verts
    ax.plot(verts[:, 0], verts[:, 1], "bo", alpha=0.25, ms=2)
    for i in range(ncpl):
        x, y = verts[i, 0], verts[i, 1]
        if xmin <= x <= xmax and ymin <= y <= ymax:
            ax.annotate(str(i + 1), verts[i, :], color="b", alpha=0.5)

    # plot nodes
    xc, yc = mg.get_xcellcenters_for_layer(0), mg.get_ycellcenters_for_layer(0)
    for i in range(ncpl):
        x, y = xc[i], yc[i]
        ax.plot(x, y, "o", color="grey", alpha=0.25, ms=2)
        if xmin <= x <= xmax and ymin <= y <= ymax:
            ax.annotate(str(i + 1), (x, y), color="grey", alpha=0.5)

    # plot well
    ax.plot(wel_coords[0][0], wel_coords[0][1], "ro")

    # adjust left margin to compensate for inset
    plt.subplots_adjust(left=0.45)

    # create legend
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Vertex",
                markerfacecolor="blue",
                markersize=10,
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Node",
                markerfacecolor="grey",
                markersize=10,
            ),
        ],
        loc="upper left",
    )


def plot_bc(sim, ibd=None):
    with styles.USGSPlot():
        fig = plt.figure(figsize=figure_size_solo)
        fig.tight_layout()
        ax = fig.add_subplot(1, 1, 1, aspect="equal")
        styles.heading(ax=ax, heading="Boundary conditions")
        gwf = sim.get_model(gwf_name)
        pmv = flopy.plot.PlotMapView(gwf, ax=ax)
        pmv.plot_bc("WEL", plotAll=True)
        pmv.plot_grid(lw=0.5)
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)
        if ibd is not None:
            pmv.plot_array(ibd, cmap=cmapbd, edgecolor="gray")

        # create inset
        axins = ax.inset_axes([-0.76, 0.25, 0.7, 0.9])
        plot_nodes_and_vertices(gwf, axins)
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
        if plot_show:
            plt.show()
        if plot_save:
            plt.savefig(figs_path / "{}-bc".format(sim_name))


def plot_head(gwf, head):
    with styles.USGSPlot():
        fig = plt.figure(figsize=figure_size_solo)
        fig.tight_layout()
        ax = fig.add_subplot(1, 1, 1, aspect="equal")
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)
        ilay = 2
        cint = 0.25
        hmin = head[ilay, 0, :].min()
        hmax = head[ilay, 0, :].max()
        styles.heading(ax=ax, heading="Head, layer {}".format(ilay + 1))
        mm = flopy.plot.PlotMapView(gwf, ax=ax, layer=ilay)
        mm.plot_grid(lw=0.5)
        pc = mm.plot_array(head[:, 0, :], edgecolor="black")
        cb = plt.colorbar(pc, shrink=0.5, pad=0.1)
        cb.ax.set_xlabel(r"Head ($m$)")

        levels = np.arange(np.floor(hmin), np.ceil(hmax) + cint, cint)
        cs = mm.contour_array(head[:, 0, :], colors="white", levels=levels)
        plt.clabel(cs, fmt="%.1f", colors="white", fontsize=11)

        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / "{}-head".format(sim_name))


def plot_points(ax, gwf, data):
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(lw=0.5, alpha=0.5)
    return ax.scatter(
        data["x"],
        data["y"],
        color=data["mc"],
        s=3,
    )


def plot_tracks(
    ax, gwf, title=None, ibd=None, pathlines=None, timeseries=None, endpoints=None
):
    ax.set_aspect("equal")
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    if title is not None:
        ax.set_title(title, fontsize=12)

    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(lw=0.5)

    if ibd is not None:
        mm.plot_array(ibd, cmap=cmapbd, edgecolor="gray")

    pts = []
    if pathlines is not None:
        pts.append(
            mm.plot_pathline(
                pathlines, layer="all", colors=["blue"], lw=0.5, ms=1, alpha=0.25
            )
        )
        if timeseries is None and endpoints is None:
            plot_points(ax, gwf, pathlines[pathlines.ireason != 1])
    if timeseries is not None:
        for k in range(0, nlay - 1):
            pts.append(mm.plot_timeseries(timeseries, layer=k, lw=0, ms=1))
            plot_points(ax, gwf, timeseries)
    if endpoints is not None:
        pts.append(
            mm.plot_endpoint(endpoints, direction="ending", colorbar=False, shrink=0.25)
        )
    return pts[0] if len(pts) == 1 else pts


def plot_2a(gwf, mf6pl, mp7pl, mp7ts, ibd, title=None):
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=figure_size_compare)
    fig.tight_layout()
    axes = axes.flatten()

    plot_tracks(axes[0], gwf, ibd=ibd, pathlines=mf6pl[mf6pl.subprob == "A"])
    plot_tracks(axes[1], gwf, ibd=ibd, pathlines=mp7pl[mp7pl.subprob == "A"], timeseries=mp7ts[mp7ts.subprob == "A"])
    plot_tracks(axes[2], gwf, ibd=ibd, pathlines=mf6pl[mf6pl.subprob == "B"])
    plot_tracks(axes[3], gwf, ibd=ibd, pathlines=mp7pl[mp7pl.subprob == "B"], timeseries=mp7ts[mp7ts.subprob == "B"])

    if title is not None:
        styles.heading(axes[0], title)
    axes[0].set_ylabel("1A")
    axes[2].set_ylabel("1B")
    axes[2].set_xlabel("MODFLOW 6 PRT")
    axes[3].set_xlabel("MODPATH 7")

    plt.subplots_adjust(wspace=-0.7)

    if plot_show:
        plt.show()
    if plot_save:
        fig.savefig(figs_path / "{}-paths".format(sim_name))


def plot_2b(gwf, mf6endpoints, mp7endpoints, ibd, title=None):
    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=figure_size_compare)
    axes = axes.flatten()

    plot_tracks(axes[0], gwf, ibd=ibd, endpoints=mf6endpoints)
    pts = plot_tracks(axes[1], gwf, ibd=ibd, endpoints=mp7endpoints)

    if title is not None:
        styles.heading(axes[0], title)
    axes[0].set_xlabel("MODFLOW 6 PRT")
    axes[1].set_xlabel("MODPATH 7")

    cax = fig.add_axes([0.2, 0.085, 0.6, 0.02])
    cb = plt.colorbar(pts, cax=cax, orientation="horizontal")
    cb.set_label("Travel time (days)")

    plt.subplots_adjust(bottom=0.2, wspace=-0.2)

    if plot_show:
        plt.show()
    if plot_save:
        fig.savefig(figs_path / "{}-endpts".format(sim_name))


def plot_3d(gwf, pathlines, endpoints=None, title=None):
    from flopy.export.vtk import Vtk
    import pyvista as pv

    pv.set_plot_theme("document")

    axes = pv.Axes(show_actor=False, actor_scale=2.0, line_width=5)

    vert_exag = 10
    pathlines = pathlines.to_records(index=False)
    pathlines["z"] = pathlines["z"] * vert_exag
    
    vtk = Vtk(model=gwf, binary=False, vertical_exageration=vert_exag, smooth=False)
    vtk.add_model(gwf)
    gwf_mesh = vtk.to_pyvista()
    pls_mesh = pv.PolyData(np.array(tuple(map(tuple, pathlines[["x", "y", "z"]]))))
    riv_mesh = pv.Box(
        bounds=[
            gwf.modelgrid.extent[1] - delc,
            gwf.modelgrid.extent[1],
            gwf.modelgrid.extent[2],
            gwf.modelgrid.extent[3],
            220 * vert_exag,
            gwf.output.head().get_data()[(0, 0, ncol - 1)] * vert_exag,
        ]
    )
    wel_mesh = pv.Box(bounds=(4700, 4800, 5200, 5300, 0, 200 * vert_exag))
    bed_mesh = pv.Box(
        bounds=[
            gwf.modelgrid.extent[0],
            gwf.modelgrid.extent[1],
            gwf.modelgrid.extent[2],
            gwf.modelgrid.extent[3],
            200 * vert_exag,
            220 * vert_exag,
        ]
    )
    gwf_mesh.rotate_z(-30, point=axes.origin, inplace=True)
    gwf_mesh.rotate_y(-10, point=axes.origin, inplace=True)
    gwf_mesh.rotate_x(10, point=axes.origin, inplace=True)
    pls_mesh.rotate_z(-30, point=axes.origin, inplace=True)
    pls_mesh.rotate_y(-10, point=axes.origin, inplace=True)
    pls_mesh.rotate_x(10, point=axes.origin, inplace=True)
    riv_mesh.rotate_z(-30, point=axes.origin, inplace=True)
    riv_mesh.rotate_y(-10, point=axes.origin, inplace=True)
    riv_mesh.rotate_x(10, point=axes.origin, inplace=True)
    wel_mesh.rotate_z(-30, point=axes.origin, inplace=True)
    wel_mesh.rotate_y(-10, point=axes.origin, inplace=True)
    wel_mesh.rotate_x(10, point=axes.origin, inplace=True)
    bed_mesh.rotate_z(-30, point=axes.origin, inplace=True)
    bed_mesh.rotate_y(-10, point=axes.origin, inplace=True)
    bed_mesh.rotate_x(10, point=axes.origin, inplace=True)

    p = pv.Plotter(window_size=[500, 500])
    if title is not None:
        p.add_title(title, font_size=5)
    p.add_mesh(gwf_mesh, opacity=0.025, style="wireframe")
    p.add_mesh(
        pls_mesh,
        scalars=pathlines.k.ravel(),
        cmap=["green", "gold", "red"],
        point_size=2,
    )
    if endpoints is not None:
        endpoints = endpoints.to_records(index=False)
        endpoints["z"] = endpoints["z"] * vert_exag
        eps_mesh = pv.PolyData(np.array(tuple(map(tuple, endpoints[["x", "y", "z"]]))))
        eps_mesh.rotate_z(-30, point=axes.origin, inplace=True)
        eps_mesh.rotate_y(-10, point=axes.origin, inplace=True)
        eps_mesh.rotate_x(10, point=axes.origin, inplace=True)
        p.add_mesh(
            eps_mesh,
            scalars=endpoints.k.ravel(),
            cmap=["green", "gold", "red"],
            point_size=2,
        )
    p.add_mesh(riv_mesh, color="teal", opacity=0.2)
    p.add_mesh(wel_mesh, color="red", opacity=0.3)
    p.add_mesh(bed_mesh, color="tan", opacity=0.1)
    p.remove_scalar_bar()
    p.add_legend(
        labels=[("Layer 1", "green"), ("Layer 2", "gold"), ("Layer 3", "red")],
        bcolor="white",
        face="r",
        size=(0.1, 0.1),
    )

    if plot_save:
        p.save_graphic(figs_path / f"{sim_name}-paths-3d.pdf", raster=False)
    if plot_show:
        p.camera.zoom(1.7)
        p.show()


def load_head():
    head_file = flopy.utils.HeadFile(gwf_ws / (gwf_name + ".hds"))
    return head_file.get_data()


def get_ibound():
    ibd = np.zeros((ncpl), dtype=int)
    ibd[np.array(welcells)] = 1
    ibd[np.array(rivcells)] = 2
    return np.ma.masked_equal(ibd, 0)


def plot_results(gwf_sim):
    gwf_model = gwf_sim.get_model(gwf_name)
    ibound = get_ibound()
    plot_bc(gwf_sim, ibound)
    plot_head(gwf_model, load_head())

    mf6pl = get_mf6_pathlines(prt_ws / trackcsvfile_prt)
    mf6ep = mf6ep = mf6pl[mf6pl.ireason == 3]  # termination event
    mp7pl = get_mp7_pathlines(mp7_ws / f"{mp7_name}.mppth", gwf_model)
    mp7ts = get_mp7_timeseries(mp7_ws / f"{mp7_name}.timeseries", gwf_model)
    mp7ep = get_mp7_endpoints(mp7_ws / f"{mp7_name}.mpend", gwf_model)
    plot_2a(
        gwf_model,
        mf6pl=mf6pl,
        mp7pl=mp7pl,
        mp7ts=mp7ts,
        title="Pathlines, 1000-day time points colored by layer",
        ibd=ibound,
    )
    plot_3d(
        gwf_model,
        pathlines=mp7pl[mp7pl.subprob == "A"],
        endpoints=mp7ep,
        title="Path points (2A) and recharge points (2B),\ncolored by layer"
    )
    plot_2b(
        gwf_model,
        mf6endpoints=mf6ep[mf6ep.subprob == "B"],
        mp7endpoints=mp7ep[mp7ep.subprob == "B"],
        title="Recharge points (2B), colored by travel time",
        ibd=ibound,
    )


# -

# ### Running the example
#
# Define and invoke a function to run the example scenario, then plot results.


# +
def scenario(silent=False):
    gwf_sim, prt_sim, mp7_sim = build_models()
    if write:
        write_models(gwf_sim, prt_sim, mp7_sim, silent=silent)
    if run:
        run_models(gwf_sim, prt_sim, mp7_sim, silent=silent)
    if plot:
        plot_results(gwf_sim)  # , prt_sim, mp7_sim)


# We are now ready to run the example problem. Subproblems 2A and 2B are solved by a single MODFLOW 6 run and a single MODPATH 7 run, so they are included under one "scenario". Each of the two subproblems is represented by its own particle release package (for MODFLOW 6) or particle group (for MODPATH 7).
scenario()
# -
