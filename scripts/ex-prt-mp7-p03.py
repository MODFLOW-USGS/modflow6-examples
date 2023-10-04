# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: venv
#     language: python
#     name: python3
# ---

# ## Particle tracking through a transient flow system
#
# Application of a MODFLOW 6 particle-tracking (PRT) model to solve example 3 from the MODPATH 7 documentation.
#
# This is a transient MODFLOW 6 simulation with a flow system similar to examples 1 and 2.
# Starting at 90,000 days, a 25 particles are released from 4 cells in the grid's top left corner.
# This occurs every 20 days for 200 days. 1000 particles are released in total.
# The model's 2 wells begin to pump at a constant rate 100,000 days into the simulation,
# and continue to pump at a constant rate for the remainder of the simulation.
#
# The original problem specifies three stress periods:
#
# | Stress period | Type         | Time steps | Length (days) | Multiplier |
# |:--------------|:-------------|:-----------|:--------------|:-----------|
# | 1             | steady-state | 1          | 100000        | 1          |
# | 2             | transient    | 10         | 36500         | 1.5        |
# | 3             | steady-state | 1          | 100000        | 1          |
#
# The original solution achieves particle release timing by setting MODPATH 7's reference time option to 2 and value to 0.9. Since PRT does not yet support fractional time steps, we modify the problems' time discretization to 5 stress periods, with a finer time discretization in the 2nd for particle release:
#
# | Stress period | Type         | Time steps | Length (days) | Multiplier |
# |:--------------|:-------------|:-----------|:--------------|:-----------|
# | 1             | steady-state | 1          | 90000         | 1          |
# | 2             | steady-state | 10         | 200           | 1          |
# | 3             | steady-state | 1          | 9800          | 1          |
# | 4             | transient    | 10         | 36500         | 1.5        |
# | 5             | steady-state | 1          | 100000        | 1          |
#

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

perioddata = [
    # perlen, nstp, tsmult
    (90000, 1, 1),
    (200, 10, 1),
    (9800, 1, 1),
    (36500, 10, 1.5),
    (100000, 1, 1),
]

#
# First import dependencies.

# +
import os
import sys
import matplotlib.pyplot as plt
import flopy
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.append(os.path.join("..", "common"))
import config
# -

# Setup model name and workspace variables.

# +
ws = config.base_ws
sim_name = "mp7-p03"
example_name = "ex-prt-" + sim_name
sim_ws = Path(ws) / example_name

nm_mf6 = sim_name
nm_prt = sim_name + "_prt"

headfile = f"{sim_name}.hds"
budgetfile = f"{sim_name}.cbb"
budgetfile_prt = f"{nm_prt}.cbb"
trackcsvfile_prt = f"{nm_prt}.trk.csv"
pathlinefile_mp7 = f"{sim_name}_mp.mppth"
endpointfile_mp7 = f"{sim_name}_mp.mpend"
# -

# Define some shared variables, starting with discretization and other flow model data.

nlay, nrow, ncol = 3, 21, 20
delr = delc = 500.0
top = 350.0
botm = [220.0, 200.0, 0.0]
laytyp = [1, 0, 0]
kh = [50.0, 0.01, 200.0]
kv = [10.0, 0.01, 20.0]
rchv = 0.005
riv_h = 320.0
riv_z = 317.0
riv_c = 1.0e5

# Define well data. Although this notebook will refer to layer/row/column indices starting at 1, indices in FloPy (and more generally in Python) are zero-based. A negative discharge indicates pumping, while a positive value indicates injection.

wells = [
    # layer, row, col, discharge
    (0, 10, 9, -75000),
    (2, 12, 4, -100000),
]

# Define the drain location.

drain = (0, 14, (9, 20))

# Configure locations for particle tracking to terminate. We have three explicitly defined termination zones:
#
# - 2: the well in layer 1, at row 11, column 10
# - 3: the well in layer 3, at row 13, column 5
# - 4: the drain in layer 1, running through row 15 from column 10-20
# - 5: the river in layer 1, running through column 20
#
# MODFLOW 6 reserves zone number 1 to indicate that particles may move freely within the zone.

# +
zone_maps = []

# zone 1 is the default (non-terminating)
def fill_zone_1():
    return np.ones((nrow, ncol), dtype=np.int32)

# zone map for layer 1
za = fill_zone_1()
za[wells[0][1:3]] = 2 # well
za[drain[1], drain[2][0] : drain[2][1]] = 4 # drain
za[:, ncol - 1] = 5 # river
zone_maps.append(za)

# zone map for layer 2 (use constant to apply zone 1 to all cells)
zone_maps.append(1)

# zone map for layer 3
za = fill_zone_1()
za[wells[1][1:3]] = 3 # well
zone_maps.append(za)
# -

# TODO: plot zone map
#
# Define particles to track. Particles are released from the top of a 2x2 square of cells in the upper left of the midel grid's top layer.

# +

rel_minl = rel_maxl = 0
rel_minr = 2
rel_maxr = 3
rel_minc = 2
rel_maxc = 3
celldata = flopy.modpath.CellDataType(
    drape=0,
    rowcelldivisions=5,
    columncelldivisions=5,
    layercelldivisions=1,
)
lrcregions = [[rel_minl, rel_minr, rel_minc, rel_maxl, rel_maxr, rel_maxc]]
lrcpd = flopy.modpath.LRCParticleData(
    subdivisiondata=celldata,
    lrcregions=lrcregions,
)
pg = flopy.modpath.ParticleGroupLRCTemplate(
    particlegroupname="PG1",
    particledata=lrcpd,
    filename=f"{sim_name}.pg1.sloc",
    releasedata=(10, 0, 20)
)
pgs = [pg]
defaultiface = {"RECHARGE": 6, "ET": 6}
# -

# Define a function to create the MODFLOW 6 simulation, which will include a groundwater flow (GWF) model and a particle tracking (PRT) model.

def build_mf6_model():
    print("Building MODFLOW 6 model")

    # simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name,
        exe_name="mf6",
        version="mf6",
        sim_ws=sim_ws
    )

    # temporal discretization
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=len(perioddata),
        perioddata=perioddata
    )

    # groundwater flow (gwf) model
    model_nam_file = f"{sim_name}.nam"
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=sim_name, model_nam_file=model_nam_file, save_flows=True
    )

    # iterative model solver (ims) package
    ims = flopy.mf6.modflow.mfims.ModflowIms(
        sim,
        pname="ims",
        complexity="SIMPLE",
    )

    # grid discretization
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        length_units="FEET",
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic", strt=320)

    # node property flow
    npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf, pname="npf", icelltype=laytyp, k=kh, k33=kv, save_flows=True, save_specific_discharge=True, save_saturation=True
    )

    # recharge
    rch = flopy.mf6.modflow.mfgwfrcha.ModflowGwfrcha(gwf, recharge=rchv)

    # storage
    sto = flopy.mf6.modflow.mfgwfsto.ModflowGwfsto(
        gwf, save_flows=True, iconvert=1, ss=0.0001, sy=0.1,
        steady_state={0: True, 4: True},
        transient={3: True},
    )

    # wells
    def no_flow(w):
        return w[0], w[1], w[2], 0

    nf_wells = [no_flow(w) for w in wells] 
    wel = flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(
        gwf,
        maxbound=2,
        stress_period_data={0: nf_wells, 1: nf_wells, 2: nf_wells, 3: wells, 4: wells},
    )

    # river
    riv_iface = 6
    riv_iflowface = -1
    rd = [[(0, i, ncol - 1), riv_h, riv_c, riv_z, riv_iface, riv_iflowface] for i in range(nrow)]
    flopy.mf6.modflow.mfgwfriv.ModflowGwfriv( 
        gwf,
        auxiliary=["iface", "iflowface"],
        stress_period_data={0: rd}
    )

    # drain (set auxiliary IFACE var to 6 for top of cell)
    drn_iface = 6
    drn_iflowface = -1
    dd = [
        [drain[0], drain[1], i + drain[2][0], 322.5, 100000.0, drn_iface, drn_iflowface]
        for i in range(drain[2][1] - drain[2][0])
    ]
    drn = flopy.mf6.modflow.mfgwfdrn.ModflowGwfdrn(
        gwf,
        auxiliary=["iface", "iflowface"],
        stress_period_data={0: dd})

    # output control
    oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
        gwf,
        pname="oc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        head_filerecord=[headfile],
        budget_filerecord=[budgetfile],
    )

    # Instantiate the MODFLOW 6 prt model
    prt = flopy.mf6.ModflowPrt(
        sim, modelname=nm_prt, model_nam_file="{}.nam".format(nm_prt)
    )

    # Instantiate the MODFLOW 6 prt discretization package
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt, pname="dis",
        nlay=nlay, nrow=nrow, ncol=ncol,
        length_units="FEET",
        delr=delr, delc=delc,
        top=top, botm=botm,
    )

    # Instantiate the MODFLOW 6 prt model input package
    porosity = 0.1
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    releasepts = list(lrcpd.to_prp(gwf.modelgrid))
    flopy.mf6.ModflowPrtprp(
        prt, pname="prp1", filename="{}_1.prp".format(nm_prt),
        nreleasepts=len(releasepts), packagedata=releasepts,
        perioddata={
            0: [],
            1: ["STEPS"] + [str(i + 1) for i in range(10)],
            2: [],
            3: [],
            4: []
        },
    )

    # Instantiate the MODFLOW 6 prt output control package
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[budgetfile_prt],
        trackcsv_filerecord=[trackcsvfile_prt],
        saverecord=[("BUDGET", "ALL")],
    )

    # Instantiate the MODFLOW 6 prt flow model interface
    flopy.mf6.ModflowPrtfmi(prt, packagedata=[
        ("GWFHEAD", headfile),
        ("GWFBUDGET", budgetfile)
    ])

    # Create the MODFLOW 6 gwf-prt model exchange
    flopy.mf6.ModflowGwfprt(
        sim, exgtype="GWF6-PRT6",
        exgmnamea=nm_mf6, exgmnameb=nm_prt,
        filename="{}.gwfprt".format(nm_mf6),
    )

    # Create an explicit model solution (EMS) for the MODFLOW 6 prt model
    ems = flopy.mf6.ModflowEms(
        sim, pname="ems",
        filename=f"{nm_prt}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


# Define a function create the MODPATH 7 model & simulation.

def build_mp7_model(gwf):
    print("Building MODPATH 7 model...")

    mp = flopy.modpath.Modpath7(
        modelname=f"{sim_name}_mp",
        flowmodel=gwf,
        exe_name="mp7",
        model_ws=sim_ws,
    )
    mpbas = flopy.modpath.Modpath7Bas(mp, porosity=0.1, defaultiface=defaultiface)
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="combined",
        trackingdirection="forward",
        weaksinkoption="pass_through",
        weaksourceoption="pass_through",
        budgetoutputoption="summary",
        referencetime=[1, 0, 0],
        timepointdata=[30, 2000.0],
        zonedataoption="on",
        zones=zone_maps,
        particlegroups=pgs,
    )

    return mp


# Define a function to build both simulations.

def build_models():
    mf6sim = build_mf6_model()
    gwf = mf6sim.get_model(nm_mf6)
    mp7sim = build_mp7_model(gwf)
    return mf6sim, mp7sim


# Define a function to run the MODFLOW 6 and MOPATH 7 models/simulations.

def run_models(mf6sim, mp7sim):
    mf6sim.write_simulation()
    success, buff = mf6sim.run_simulation(silent=True, report=True)
    for line in buff:
        print(line)
    assert success, "Failed to run MODFLOW 6"

    mp7sim.write_input()
    success, buff = mp7sim.run_model(silent=True, report=True)
    for line in buff:
        print(line)
    assert success, "Failed to run MODPATH 7"


# Define a function to determine termination zones.

def get_term_zone_nns(grid):
    wel_locs = [w[0:3] for w in wells]
    riv_locs = [(0, i, 19) for i in range(20)]
    drn_locs = [(drain[0], drain[1], d) for d in range(drain[2][0], drain[2][1])]
    wel_nids = grid.get_node(wel_locs)
    riv_nids = grid.get_node(riv_locs)
    drn_nids = grid.get_node(drn_locs)
    return wel_nids, drn_nids, riv_nids


# Define functions to load pathline data from the budget file created by the PRT model.

# +
from flopy.plot.plotutil import to_mp7_pathlines, to_mp7_endpoints

def load_mf6_pathlines(p):
    # read pathlines into dataframe from csv
    pls = pd.read_csv(p)

    # assign a unique particle index column incrementing an integer
    # for each unique combination of irpt, iprp, imdl, and trelease
    id_cols = ["imdl", "iprp", "irpt", "trelease"]
    pls = pls.sort_values(id_cols)
    particles = pls.groupby(id_cols)
    pls["particleid"] = particles.ngroup()

    # select endpoints
    eps = pls.sort_values("t").groupby(id_cols).tail(1)

    # add a terminating zone column to pathlines dataframe
    def set_term_zone(row):
        id = row["particleid"]
        ep = eps[eps["particleid"] == id].iloc[0]
        return ep["izone"]
    pls["izone_term"] = pls.apply(lambda r: set_term_zone(r), axis=1)

    # return a dict keyed on capture zone
    return {
        "well": to_mp7_pathlines(pls[(pls["izone_term"] == 2) | (pls["izone_term"] == 3)]),
        "drain": to_mp7_pathlines(pls[pls["izone_term"] == 4]),
        "river": to_mp7_pathlines(pls[pls["izone_term"] == 5]),
    }

def load_mf6_endpoints(p):
    # read pathlines into dataframe from csv
    pls = pd.read_csv(p)

    # assign a unique particle index column incrementing an integer
    # for each unique combination of irpt, iprp, imdl, and trelease
    id_cols = ["imdl", "iprp", "irpt", "trelease"]
    pls = pls.sort_values(id_cols)
    particles = pls.groupby(id_cols)
    pls["particleid"] = particles.ngroup()

    # select endpoints
    eps = pls.sort_values("t").groupby(id_cols).tail(1)

    # add a terminating zone column to pathlines dataframe
    def set_term_zone(row):
        id = row["particleid"]
        ep = eps[eps["particleid"] == id].iloc[0]
        return ep["izone"]
    pls["izone_term"] = pls.apply(lambda r: set_term_zone(r), axis=1)

    # return a dict keyed on capture zone
    return {
        "all": to_mp7_endpoints(pls),
        "well": to_mp7_endpoints(pls[(pls["izone_term"] == 2) | (pls["izone_term"] == 3)]),
        "drain": to_mp7_endpoints(pls[pls["izone_term"] == 4]),
        "river": to_mp7_endpoints(pls[pls["izone_term"] == 5]),
    }
# -

# Define functions to load pathline and endpoint data from the MODPATH 7 models' output files.

# +
from flopy.utils.modpathfile import EndpointFile

def load_mp7_pathlines(pth, grid):
    wel_nids, drn_nids, riv_nids = get_term_zone_nns(grid)
    mp7pathlines = {}
    p = flopy.utils.PathlineFile(pth)
    mp7pathlines["well"] = p.get_destination_pathline_data(wel_nids, to_recarray=True)
    mp7pathlines["drain"] = p.get_destination_pathline_data(drn_nids, to_recarray=True)
    mp7pathlines["river"] = p.get_destination_pathline_data(riv_nids, to_recarray=True)
    return mp7pathlines

def load_mp7_endpoints(pth, grid):
    wel_nids, drn_nids, riv_nids = get_term_zone_nns(grid)
    mp7endpoints = {}
    p = EndpointFile(pth)
    # endpts = pthobj.get_alldata()
    mp7endpoints["well"] = p.get_destination_endpoint_data(wel_nids)
    mp7endpoints["drain"] = p.get_destination_endpoint_data(drn_nids)
    mp7endpoints["river"] = p.get_destination_endpoint_data(riv_nids)
    return mp7endpoints


# -

# Define functions to plot the heads and pathline results, with pathlines colored by capture zone.

# +
from flopy.export.vtk import Vtk
import pyvista as pv

def plot_startpoints(ax, gwf, sp, title):
    ax.set_aspect("equal")
    ax.set_title(title)
    mv = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mv.plot_grid(lw=0.5)
    mv.plot_bc("WEL", plotAll=True)
    mv.plot_bc("DRN")
    mv.plot_bc("RIV")
    mv.plot_endpoint(sp, label="start points", color="yellow")
    ax.set_xlim([1000, 2000])
    ax.set_ylim([9000, 10000])

def plot_endpoints(ax, gwf, ep, title):
    ax.set_aspect("equal")
    ax.set_title(title)

    mv = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mv.plot_grid(lw=0.5)
    mv.plot_bc("DRN", alpha=0.2)
    mv.plot_bc("RIV", alpha=0.2)
    mv.plot_bc("WEL", alpha=0.2, plotAll=True)

    if ep["well"] is not None:
        mv.plot_endpoint(ep["well"], label="end points", color="red")
    if ep["drain"] is not None:
        mv.plot_endpoint(ep["drain"], label="end points", color="green")
    if ep["river"] is not None:
        mv.plot_endpoint(ep["river"], label="end points", color="blue")

    # inset axes, zoom into particle endpoint clusters, first well..
    axins1 = ax.inset_axes([-0.82, 0.25, 0.5, 0.5])
    mvins1 = flopy.plot.PlotMapView(model=gwf, ax=axins1)
    mvins1.plot_grid(lw=0.5)

    if ep["well"] is not None:
        mvins1.plot_endpoint(ep["well"], label="end points", color="red")
    if ep["drain"] is not None:
        mvins1.plot_endpoint(ep["drain"], label="end points", color="green")
    if ep["river"] is not None:
        mvins1.plot_endpoint(ep["river"], label="end points", color="blue")

    axins1.set_xlim([2000, 3000])
    axins1.set_ylim([3900, 4600])
    ax.indicate_inset_zoom(axins1)

    # ..then river
    axins2 = ax.inset_axes([1.12, 0.25, 0.5, 0.5])
    mvins2 = flopy.plot.PlotMapView(model=gwf, ax=axins2)
    mvins2.plot_grid(lw=0.5)

    if ep["well"] is not None:
        mvins2.plot_endpoint(ep["well"], label="end points", color="red")
    if ep["drain"] is not None:
        mvins2.plot_endpoint(ep["drain"], label="end points", color="green")
    if ep["river"] is not None:
        mvins2.plot_endpoint(ep["river"], label="end points", color="blue")

    axins2.set_xlim([9000, 10000])
    axins2.set_ylim([2500, 3900])
    ax.indicate_inset_zoom(axins2)

def plot_pathlines(ax, gwf, hd, pl, title):
    ax.set_aspect("equal")
    ax.set_title(title)
    mv = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mv.plot_grid(lw=0.5)
    mv.plot_bc("WEL", plotAll=True)
    mv.plot_bc("DRN")
    mv.plot_bc("RIV")
    hd = mv.plot_array(hd, alpha=0.2)
    cb = plt.colorbar(hd, shrink=0.5, ax=ax)
    cb.set_label("Head")
    mv.plot_pathline(
        pl['well'],
        layer="all",
        colors=["red"],
        label="captured by well",
    )
    mv.plot_pathline(
        pl['drain'],
        layer="all",
        colors=["green"],
        label="captured by drain",
    )
    mv.plot_pathline(
        pl['river'],
        layer="all",
        colors=["blue"],
        label="captured by river",
    )
    mv.ax.legend()

    # inset axes, zoom into particle release locations....
    axins = ax.inset_axes([-0.82, 0.5, 0.5, 0.5])
    mvins = flopy.plot.PlotMapView(model=gwf, ax=axins)
    mvins.plot_grid(lw=0.5)
    mvins.plot_pathline(
        pl['well'],
        layer="all",
        colors=["red"],
        label="captured by well",
    )
    mvins.plot_pathline(
        pl['drain'],
        layer="all",
        colors=["green"],
        label="captured by drain",
    )
    mvins.plot_pathline(
        pl['river'],
        layer="all",
        colors=["blue"],
        label="captured by river",
    )
    axins.set_xlim([1000, 2000])
    axins.set_ylim([8500, 9500])
    ax.indicate_inset_zoom(axins)

def plot_pathlines_xc(ax, gwf, hd, pl, line, title):
    # ax.set_aspect("equal")
    ax.set_title(title)
    mv = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line=line)
    mv.plot_grid(lw=0.5)
    mv.plot_bc("WEL")
    mv.plot_bc("DRN")
    mv.plot_bc("RIV")
    hd = mv.plot_array(hd, alpha=0.2)
    cb = plt.colorbar(hd, shrink=0.5, ax=ax)
    cb.set_label("Head")
    mv.plot_pathline(
        pl['well'],
        colors=["red"],
        label="captured by well",
    )
    mv.plot_pathline(
        pl['drain'],
        colors=["green"],
        label="captured by drain",
    )
    mv.plot_pathline(
        pl['river'],
        colors=["blue"],
        label="captured by river",
    )
    mv.ax.legend()

def plot_results_2d(gwf, heads, ep1, ep2, pl1, pl2, title1="MODFLOW 6 PRT", title2="MODPATH 7"):
    fig, axes = plt.subplots(ncols=2, nrows=5, figsize=(18, 18))
    plt.suptitle(
        t="Example 3: Structured transient forward-tracking model with multiple release times",
        fontsize=14,
        y=0.92,
    )

    # map view pathlines
    plot_pathlines(axes[0][0], gwf, heads, pl1, f"{title1} map view")
    plot_pathlines(axes[0][1], gwf, heads, pl2, f"{title2} map view")

    # map view end points
    plot_endpoints(axes[1][0], gwf, ep1, f"{title1} end points")
    plot_endpoints(axes[1][1], gwf, ep2, f"{title2} end points")

    # cross sections
    xc_row = 11
    xc_line = {"row": xc_row}
    plot_pathlines_xc(axes[2][0], gwf, heads, pl1, xc_line, f"{title1} cross section, row {xc_row}")
    plot_pathlines_xc(axes[2][1], gwf, heads, pl2, xc_line, f"{title2} cross section, row {xc_row}")
    xc_row = 12
    xc_line = {"row": xc_row}
    plot_pathlines_xc(axes[3][0], gwf, heads, pl1, xc_line, f"{title1} cross section, row {xc_row}")
    plot_pathlines_xc(axes[3][1], gwf, heads, pl2, xc_line, f"{title2} cross section, row {xc_row}")
    xc_row = 14
    xc_line = {"row": xc_row}
    plot_pathlines_xc(axes[4][0], gwf, heads, pl1, xc_line, f"{title1} cross section, row {xc_row}")
    plot_pathlines_xc(axes[4][1], gwf, heads, pl2, xc_line, f"{title2} cross section, row {xc_row}")
    
    plt.show()

def plot_results_3d(gwf, heads, path1, path2, title1="MODFLOW 6 PRT", title2="MODPATH 7"):
    # read the VTK files into PyVista meshes
    grid = pv.read(path1.parent / f"{path1.stem}_000000.vtk")
    mf6_pls = pv.read(path1.parent / f"{path1.stem}_pathline.vtk")
    mp7_pls = pv.read(path2.parent / f"{path2.stem}_pathline.vtk")

    # rotate and scale the plot
    axes = pv.Axes(show_actor=True, actor_scale=2.0, line_width=5)
    grid.rotate_z(160, point=axes.origin, inplace=True)
    mf6_pls.rotate_z(160, point=axes.origin, inplace=True)
    mp7_pls.rotate_z(160, point=axes.origin, inplace=True)

    # check grid properties
    assert grid.n_cells == gwf.modelgrid.nnodes
    print("Grid has", grid.n_cells, "cells")
    print("Grid has", grid.n_arrays, "arrays")

    # set PyVista theme and backend
    pv.set_plot_theme("document")
    pv.set_jupyter_backend('trame')
    # pv.set_jupyter_backend('static')

    # create the plot and add the grid and pathline meshes
    p = pv.Plotter(shape=(1, 2))

    # todo: when PRT updated to output capture zone,
    # add scalars and use them to color pathline pts

    # add subplot for MODFLOW 6 PRT
    p.subplot(0, 0)
    p.add_text(title1, font_size=20)
    p.add_mesh(grid, opacity=0.05)
    p.add_mesh(mf6_pls)

    # add subplot for MODPATH 7
    p.subplot(0, 1)
    p.add_text(title2, font_size=20)
    p.add_mesh(grid, opacity=0.05)
    p.add_mesh(mp7_pls)

    # zoom in and show the plot
    p.camera.zoom(2.4)
    p.show()


# -

# Define a function to wrap the entire scenario.

# +
import flopy.utils.binaryfile as bf

def scenario():
    # build models
    mf6sim, mp7sim = build_models()
    gwf = mf6sim.get_model(nm_mf6)
    grid = gwf.modelgrid

    # run models
    if not runModel:
        return

    run_models(mf6sim, mp7sim)

    # load results
    hds = flopy.utils.HeadFile(sim_ws / headfile).get_data()
    mf6_pl = load_mf6_pathlines(sim_ws / trackcsvfile_prt)
    mf6_ep = load_mf6_endpoints(sim_ws / trackcsvfile_prt)
    mp7_pl = load_mp7_pathlines(sim_ws / pathlinefile_mp7, grid)
    mp7_ep = load_mp7_endpoints(sim_ws / endpointfile_mp7, grid)

    # export pathline results to VTK files
    vtk = Vtk(model=gwf, binary=True, vertical_exageration=True, smooth=False)
    vtk.add_model(gwf)
    vtk.add_pathline_points(mf6_pl["well"] + mf6_pl["river"])
    mf6_vtk_path = mf6sim.sim_path / f"{sim_name}_mf6.vtk"
    vtk.write(mf6_vtk_path)
    vtk = Vtk(model=gwf, binary=True, vertical_exageration=True, smooth=False)
    vtk.add_model(gwf)
    vtk.add_pathline_points([mp7_pl["well"], mp7_pl["river"]])
    mp7_vtk_path = mf6sim.sim_path / f"{sim_name}_mp7.vtk"
    vtk.write(mp7_vtk_path)

    # plot results
    plot_results_2d(gwf, hds, mf6_ep, mp7_ep, mf6_pl, mp7_pl)
    plot_results_3d(gwf, hds, mf6_vtk_path, mp7_vtk_path)


# -

# Run the MODPATH 7 example problem 3 scenario.

scenario()

# +
# wdir = Path("/Users/wes/dev/modpath-v7/examples/ex03_mf6")
# 
# mf6sim, mp7sim = build_models()
# gwf = mf6sim.get_model(nm_mf6)
# grid = gwf.modelgrid
# 
# # load results
# hds = flopy.utils.HeadFile(sim_ws / headfile).get_data()
# cbb = bf.CellBudgetFile(sim_ws / budgetfile_prt)
# 
# mf6_ep = load_mf6_endpoints(cbb, gwf)
# mp7_ep = load_mp7_endpoints(sim_ws / endpointfile_mp7, grid)
# mf6_pl = load_mf6_pathlines(cbb, gwf)
# mp7_pl = load_mp7_pathlines(sim_ws / pathlinefile_mp7, grid)
# 
# ref_ep = load_mp7_endpoints(wdir / "ex03a_mf6.endpoint7", grid)
# ref_pl = load_mp7_pathlines(wdir / "ex03a_mf6.pathline7", grid)
# 
# plot_results_2d(gwf, hds, ref_ep, mp7_ep, ref_pl, mp7_pl, "MODPATH 7 (reference)", "MODPATH 7")
