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

# ## Particle tracking: steady state, structured grid
#
# Application of a MODFLOW 6 particle-tracking (PRT) model and a MODPATH 7 (MP7) model to solve example 1 from the MODPATH 7 documentation.
#
# This example problem involves a flow system consisting of two aquifers separated by a low conductivity confining layer, modeled on a nearly square structured grid with uniform square cells. There is a single well in the center of the grid, with a river running along the grid's right-hand boundary.
#
# In part A, 21 particles are released at the water table in layer 1, all along the grid's third column, and tracked until discharge locations are reached. Some of the particles discharge to the well, while some discharge to the river.
#
# In part B, 9 particles are released from points evenly distributed over the top faces of all cells in layer 1, then the capture zones of the river and the central well are computed.
#

# ### Initial setup
#
# Import dependencies, define the example name and workspace, and read settings from environment variables.

# +
import pathlib as pl
from pprint import pformat

import flopy
import git
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flopy.plot.styles import styles
import matplotlib as mpl
from modflow_devtools.misc import get_env, timed

# Example name and workspace paths. If this example is running
# in the git repository, use the folder structure described in
# the README. Otherwise just use the current working directory.
sim_name = "mp7-p01"
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
mf6_ws = sim_ws / "mf6"
mp7_ws = sim_ws / "mp7"
mf6_ws.mkdir(exist_ok=True, parents=True)
mp7_ws.mkdir(exist_ok=True, parents=True)

# Define output file names
headfile = f"{gwf_name}.hds"
budgetfile = f"{gwf_name}.cbb"
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
perlen = 1.0
tsmult = 1.0

# parse parameter strings into tuples
botm = [float(value) for value in botm_str.split(",")]
# -

# Define data for the MODFLOW 6 well, river, and recharge packages.

# +
# Well package
wel_loc = (2, 10, 9)
wd = [(wel_loc, wel_q)]

# River package
riv_iface = 6
riv_iflowface = -1
rd = []
for i in range(nrow):
    rd.append([(0, i, ncol - 1), riv_h, riv_c, riv_z, riv_iface, riv_iflowface])

# Recharge package
rch_iface = 6
rch_iflowface = -1
# -

# Define shared MODFLOW 6 PRT and MODPATH 7 particle-tracking model parameters.

# +
# Zones
def_zone = 1
wel_zone = 2
zones_lay1 = def_zone
zones_lay2 = def_zone
zones_lay3 = np.full((nrow, ncol), def_zone, dtype=np.int32)
zones_lay3[wel_loc[1:]] = wel_zone

# Starting location template size for example 1:
# 3x3 array of particles in each cell in layer 1
sloc_tmpl_size = 3
# -

# Define MODPATH 7 model parameters.

# +
# Zones
zones = [zones_lay1, zones_lay2, zones_lay3]

# Default iface
defaultiface = {"RCH": 6, "EVT": 6}
# -

# Define particle release point data for the MODFLOW 6 PRT particle-tracking model.

# +
# Example 1A
releasepts = {}
releasepts["1A"] = []
zrpt = top
k = 0
j = 2
for i in range(nrow):
    nrpt = i
    xrpt = (j + 0.5) * delr
    yrpt = (nrow - i - 0.5) * delc
    rpt = [nrpt, k, i, j, xrpt, yrpt, zrpt]
    releasepts["1A"].append(rpt)

# Example 1B
releasepts["1B"] = []
ndivc = sloc_tmpl_size
ndivr = sloc_tmpl_size
deldivc = delc / ndivc
deldivr = delr / ndivr
k = 0
zrpt = top
nrpt = -1
for i in range(nrow):
    y0 = (nrow - i - 1) * delc
    for j in range(ncol):
        x0 = j * delr
        for idiv in range(ndivc):
            dy = (idiv + 0.5) * deldivc
            yrpt = y0 + dy
            for jdiv in range(ndivr):
                dx = (jdiv + 0.5) * deldivr
                xrpt = x0 + dx
                nrpt += 1
                rpt = [nrpt, k, i, j, xrpt, yrpt, zrpt]
                releasepts["1B"].append(rpt)
# -

# Define particle data for the MODPATH 7 model

# +
# Example 1A
plocs = []
pids = []
for idx in range(nrow):
    plocs.append((0, idx, 2))
    pids.append(idx)
# issue(flopyex): in the flopy example this notebook is based on,
#  localz is not set to 1.0 like in the MODPATH examples doc,
#  so it defaults to 0.5, but it shouldn't really matter because
#  the particle gets placed at the water table anyway
part0 = flopy.modpath.ParticleData(
    plocs, drape=0, structured=True, particleids=pids, localz=1.0
)

# Example 1B
divs = sloc_tmpl_size
locs1b = [[0, 0, 0, 0, nrow - 1, ncol - 1]]
sd = flopy.modpath.CellDataType(
    drape=0, columncelldivisions=1, rowcelldivisions=1, layercelldivisions=1
)
sd = flopy.modpath.FaceDataType(
    drape=0,
    verticaldivisions1=0,
    horizontaldivisions1=0,
    verticaldivisions2=0,
    horizontaldivisions2=0,
    verticaldivisions3=0,
    horizontaldivisions3=0,
    verticaldivisions4=0,
    horizontaldivisions4=0,
    rowdivisions5=0,
    columndivisions5=0,
    rowdivisions6=divs,
    columndivisions6=divs,
)
p = flopy.modpath.LRCParticleData(subdivisiondata=[sd], lrcregions=[locs1b])
# -

# Define well and river cell numbers, used to extract and plot model results later.

# Get well and river cell numbers
nodes = {}
k, i, j = wel_loc
nodes["well"] = ncol * (nrow * k + i) + j
nodes["river"] = []
for rivspec in rd:
    k, i, j = rivspec[0]
    node = ncol * (nrow * k + i) + j
    nodes["river"].append(node)


# ### Model setup
#
# Define functions to build models, write input files, and run the simulation.


# +
def build_model(example_name):
    print("Building models...{}".format(example_name))

    # Instantiate the MODFLOW 6 simulation object
    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name, exe_name="mf6", version="mf6", sim_ws=mf6_ws
    )

    # Instantiate the MODFLOW 6 temporal discretization package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, tsmult)],
    )

    # Instantiate the MODFLOW 6 gwf (groundwater-flow) model
    model_nam_file = "{}.nam".format(gwf_name)
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwf_name, model_nam_file=model_nam_file, save_flows=True
    )

    # Instantiate the MODFLOW 6 gwf discretization package
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
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

    # Instantiate the MODFLOW 6 gwf initial conditions package
    flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic", strt=top)

    # Instantiate the MODFLOW 6 gwf node property flow package
    flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        pname="npf",
        icelltype=[1, 0, 0],
        k=kh,
        k33=kv,
        save_saturation=True,
        save_specific_discharge=True,
    )

    # Instantiate the MODFLOW 6 gwf recharge package
    flopy.mf6.modflow.mfgwfrcha.ModflowGwfrcha(
        gwf,
        recharge=rch,
        auxiliary=["iface", "iflowface"],
        aux=[rch_iface, rch_iflowface],
    )

    # Instantiate the MODFLOW 6 gwf well package
    flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(
        gwf, maxbound=1, stress_period_data={0: wd}
    )

    # Instantiate the MODFLOW 6 gwf river package
    flopy.mf6.modflow.mfgwfriv.ModflowGwfriv(
        gwf, auxiliary=["iface", "iflowface"], stress_period_data={0: rd}
    )

    # Instantiate the MODFLOW 6 gwf output control package
    head_record = [headfile]
    budget_record = [budgetfile]
    saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]
    flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
        gwf,
        pname="oc",
        saverecord=saverecord,
        head_filerecord=head_record,
        budget_filerecord=budget_record,
    )

    # Instantiate the MODFLOW 6 prt model
    prt = flopy.mf6.ModflowPrt(
        sim, modelname=prt_name, model_nam_file="{}.nam".format(prt_name)
    )

    # Instantiate the MODFLOW 6 prt discretization package
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
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

    # Instantiate the MODFLOW 6 prt model input package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # Instantiate the MODFLOW 6 prt particle release point (prp) package for example 1A
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1a",
        filename="{}_1a.prp".format(prt_name),
        nreleasepts=len(releasepts["1A"]),
        packagedata=releasepts["1A"],
        perioddata={
            0: ["FIRST"],
        },
    )

    # Instantiate the MODFLOW 6 prt particle release point (prp) package for example 1B
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1b",
        filename="{}_1b.prp".format(prt_name),
        nreleasepts=len(releasepts["1B"]),
        packagedata=releasepts["1B"],
        perioddata={
            0: ["FIRST"],
        },
    )

    # Instantiate the MODFLOW 6 prt output control package
    budget_record = [budgetfile_prt]
    track_record = [trackfile_prt]
    trackcsv_record = [trackcsvfile_prt]
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=budget_record,
        track_filerecord=track_record,
        trackcsv_filerecord=trackcsv_record,
        saverecord=[("BUDGET", "ALL")],
    )

    # Instantiate the MODFLOW 6 prt flow model interface
    flopy.mf6.ModflowPrtfmi(prt)

    # Create the MODFLOW 6 gwf-prt model exchange
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwf_name,
        exgmnameb=prt_name,
        filename="{}.gwfprt".format(sim_name),
    )

    # Create an iterative model solution (IMS) for the MODFLOW 6 gwf model
    ims = flopy.mf6.ModflowIms(
        sim,
        pname="ims",
        complexity="SIMPLE",
        outer_dvclose=1e-6,
        inner_dvclose=1e-6,
        rcloserecord=1e-6,
    )

    # Create an explicit model solution (EMS) for the MODFLOW 6 prt model
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename="{}.ems".format(prt_name),
    )
    sim.register_solution_package(ems, [prt.name])

    # Instantiate the MODPATH 7 simulation object
    mp7 = flopy.modpath.Modpath7(
        modelname=mp7_name,
        flowmodel=gwf,
        exe_name="mp7",
        model_ws=mp7_ws,
        budgetfilename=budgetfile,
        headfilename=headfile,
    )

    # Instantiate the MODPATH 7 basic data
    flopy.modpath.Modpath7Bas(mp7, porosity=porosity, defaultiface=defaultiface)

    # Instantiate the MODPATH 7 particle groups
    pg1a = flopy.modpath.ParticleGroup(
        particlegroupname="PG1A", particledata=part0, filename=sim_name + "a.sloc"
    )
    pg1b = flopy.modpath.ParticleGroupLRCTemplate(
        particlegroupname="PG1B", particledata=p, filename=sim_name + "b.sloc"
    )

    # Instantiate the MODPATH 7 simulation data
    mpsim = flopy.modpath.Modpath7Sim(
        mp7,
        simulationtype="combined",
        trackingdirection="forward",
        weaksinkoption="pass_through",
        weaksourceoption="pass_through",
        budgetoutputoption="summary",
        referencetime=[0, 0, 0.0],
        stoptimeoption="extend",
        timepointdata=[500, 1000.0],
        zonedataoption="on",
        zones=zones,
        particlegroups=[pg1a, pg1b],
    )

    return sim, mp7


def write_model(sim, mp7, silent=True):
    sim.write_simulation(silent=silent)
    mp7.write_input()


@timed
def run_model(sim, mp7, silent=True):
    # Run MODFLOW 6 simulation
    success, buff = sim.run_simulation(silent=silent, report=True)
    assert success, pformat(buff)

    # Run MODPATH 7 simulation
    success, buff = mp7.run_model(silent=silent, report=True)
    assert success, pformat(buff)


# -


# Define functions to load pathline and endpoint data from MODFLOW 6 PRT's particle track CSV files. Note that unlike MODPATH 7, MODFLOW 6 PRT does not make a distinction between pathline and endpoint output files &mdash; all pathline data is saved to track files, and endpoints are computed dynamically.

# +
from flopy.plot.plotutil import to_mp7_endpoints, to_mp7_pathlines


def load_mf6pathlines(p):
    pls = pd.read_csv(p)
    mp7_pls = to_mp7_pathlines(pls)
    mp7_pls_a = mp7_pls[mp7_pls["particlegroup"] == 1]
    mp7_pls_b = mp7_pls[mp7_pls["particlegroup"] == 2]
    mp7_eps = to_mp7_endpoints(pls)
    mp7_eps_a = mp7_eps[mp7_eps["particlegroup"] == 1]
    mp7_eps_b = mp7_eps[mp7_eps["particlegroup"] == 2]

    # determine which particles ended up in which capture area
    # (for now, use x coordinate of endpoint, todo: use izone)
    wel_pids_a = mp7_eps_a[mp7_eps_a["x"] <= 9000]["particleid"].unique()
    riv_pids_a = mp7_eps_a[mp7_eps_a["x"] > 9000]["particleid"].unique()
    wel_pids_b = mp7_eps_b[mp7_eps_b["x"] <= 9000]["particleid"].unique()
    riv_pids_b = mp7_eps_b[mp7_eps_b["x"] > 9000]["particleid"].unique()

    # return a dict keyed by subproblem (particle group) and capture area
    return {
        ("1A", "well"): mp7_pls_a[mp7_pls_a["particleid"].isin(wel_pids_a)],
        ("1A", "river"): mp7_pls_a[mp7_pls_a["particleid"].isin(riv_pids_a)],
        ("1A", "all"): mp7_pls_a,
        ("1B", "well"): mp7_pls_b[(mp7_pls_b["particleid"].isin(wel_pids_b))],
        ("1B", "river"): mp7_pls_b[(mp7_pls_b["particleid"].isin(riv_pids_b))],
        ("1B", "all"): mp7_pls_b,
    }


def load_mf6endpoints(p):
    pls = pd.read_csv(p)
    mp7_eps = to_mp7_endpoints(pls)
    mp7_eps_a = mp7_eps[mp7_eps["particlegroup"] == 1]
    mp7_eps_b = mp7_eps[mp7_eps["particlegroup"] == 2]

    # determine which particles ended up in which capture area
    # (for now, use x coordinate of endpoint, todo: use izone)
    wel_pids_a = mp7_eps_a[mp7_eps_a["x"] <= 9000]["particleid"].unique()
    riv_pids_a = mp7_eps_a[mp7_eps_a["x"] > 9000]["particleid"].unique()
    wel_pids_b = mp7_eps_b[mp7_eps_b["x"] <= 9000]["particleid"].unique()
    riv_pids_b = mp7_eps_b[mp7_eps_b["x"] > 9000]["particleid"].unique()

    # return a dict keyed by subproblem (particle group) and capture area
    return {
        ("1A", "well"): mp7_eps_a[mp7_eps_a["particleid"].isin(wel_pids_a)],
        ("1A", "river"): mp7_eps_a[mp7_eps_a["particleid"].isin(riv_pids_a)],
        ("1A", "all"): mp7_eps_a,
        ("1B", "well"): mp7_eps_b[mp7_eps_b["particleid"].isin(wel_pids_b)],
        ("1B", "river"): mp7_eps_b[mp7_eps_b["particleid"].isin(riv_pids_b)],
        ("1B", "all"): mp7_eps_b,
    }


def load_mp7pathlines(plf):
    mp7pathlines = {}
    for dest in ["well", "river"]:
        pdata = plf.get_destination_pathline_data(nodes[dest], to_recarray=True)
        for pgidx, subprob in enumerate(["1A", "1B"]):
            mp7pathlines[(subprob, dest)] = np.array(
                [point for point in pdata if point["particlegroup"] == pgidx],
                dtype=pdata.dtype,
            )

    return mp7pathlines


def load_mp7endpointdata(epf):
    mp7endpointdata = {}
    for dest in ["well", "river"]:
        epd = epf.get_destination_endpoint_data(dest_cells=nodes[dest])
        for pgidx, subprob in enumerate(["1A", "1B"]):
            mp7endpointdata[(subprob, dest)] = np.array(
                [point for point in epd if point["particlegroup"] == pgidx],
                dtype=epd.dtype,
            )
    for subprob in ["1A", "1B"]:
        mp7endpointdata[(subprob, "all")] = np.concatenate(
            (mp7endpointdata[(subprob, "well")], mp7endpointdata[(subprob, "river")])
        )

    return mp7endpointdata


# -

# ### Plotting results
#
# Define functions to plot model results.


# +
# Pathline and starting point colors by capture destination
colordest = dict.fromkeys(["well", "river"], "blue")
colordest["well"] = "red"


def plot_pathlines(ax, gwf, pathlines, subprob, plottitle=None, legend=False, **kwargs):
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(lw=0.5)
    if plottitle:
        ax.set_title(plottitle, fontsize=12)

    for dest in ["well", "river"]:
        label = None
        if "colordest" in kwargs:
            # Use colordest to color pathlines by capture destination
            color = kwargs["colordest"][dest]
            label = "captured by " + dest
        elif "color" in kwargs:
            # Use the specified color for all pathlines
            color = kwargs["color"]
        else:
            # Use a default color for all pathlines
            color = "blue"
        mm.plot_pathline(
            pathlines[(subprob, dest)], layer="all", colors=[color], label=label
        )
        if legend and label != None:
            ax.legend(loc="lower right", bbox_to_anchor=(0.3, -0.3))


def plot_endpoints(fig, ax, gwf, pointdata, subprob, plottitle, starting, legend=False, **kwargs):
    ax.set_title(plottitle, fontsize=12)
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(lw=0.5)
    startingpoint_markersize = 5
    if "colordest" in kwargs:
        for dest in ["well", "river"]:
            color = kwargs["colordest"][dest]
            label = "captured by " + dest
            if "PRT" in plottitle:
                ax.scatter(
                    pointdata[(subprob, dest)]["x0" if starting else "x"],
                    pointdata[(subprob, dest)]["y0" if starting else "y"],
                    s=startingpoint_markersize
                    ** 2,  # todo: why does flopy plot_endpoint square s?
                    color=color,
                    label=label,
                )
            else:
                mm.plot_endpoint(
                    pointdata[(subprob, dest)],
                    direction="starting" if starting else "ending",
                    s=startingpoint_markersize,
                    color=color,
                    label=label,
                )
            if legend:
                ax.legend(loc="lower right", bbox_to_anchor=(0.3, -0.3))
    else:
        if "PRT" in plottitle:
            pts = ax.scatter(
                pointdata[(subprob, "all")]["x0" if starting else "x"],
                pointdata[(subprob, "all")]["y0" if starting else "y"],
                s=startingpoint_markersize
                ** 2,  # todo: why does flopy plot_endpoint square s?
                c=pointdata[(subprob, "all")]["time"],
            )
            cax = fig.add_axes([0.05, 0.2, 0.9, 0.01])
            cb = plt.colorbar(pts, cax=cax, orientation='horizontal')
            cb.set_label("Travel time")
        else:
            mm.plot_endpoint(
                pointdata[(subprob, "all")],
                direction="starting" if starting else "ending",
                s=startingpoint_markersize,
            )


def plot_grid(gwf):
    fig = plt.figure(figsize=(8, 8))
    fig.tight_layout()
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(gwf, ax=ax)
    cb = pmv.plot_array(gwf.output.head().get_data(), alpha=0.2)
    cbar = plt.colorbar(cb, shrink=0.5)
    cbar.ax.set_xlabel(r"Head, ($m$)")
    pmv.plot_bc("WEL", plotAll=True)
    pmv.plot_bc("RIV", plotAll=True)
    pmv.plot_grid(lw=0.5)
    ax.legend(
        handles=[
            mpl.patches.Patch(color="red", label="Well"),
            mpl.patches.Patch(color="teal", label="River"),
        ],
        loc="upper left",
    )
    plt.suptitle(t="Example 1 head and boundary conditions", fontsize=14, ha="center")
    if plot_show:
        plt.show()
    if plot_save:
        fig.savefig(figs_path / f"{sim_name}-grid.png")


def plot_paths(gwf, mf6pl, mp7pl):
    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))
        plt.suptitle(
            t="Example 1 pathlines, colored by destination",
            fontsize=14,
        )
        axes = axes.flatten()
        plot_pathlines(axes[0], gwf, mf6pl, "1A", "MODFLOW 6 PRT", colordest=colordest)
        plot_pathlines(axes[1], gwf, mp7pl, "1A", "MODPATH 7", colordest=colordest)
        plot_pathlines(axes[2], gwf, mf6pl, "1B", legend=True, colordest=colordest)
        plot_pathlines(axes[3], gwf, mp7pl, "1B", colordest=colordest)
        plt.subplots_adjust(hspace=0.1)
        axes[0].set_ylabel("1A")
        axes[2].set_ylabel("1B")
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-paths.png")


def plot_release(gwf, mf6ep, mp7ep):
    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8, 8), constrained_layout=True)
        plt.suptitle(
            t="Example 1B release points, colored by destination",
            fontsize=14,
            y=0.8
        )
        plt.tight_layout()
        axes = axes.flatten()
        plot_endpoints(
            fig, axes[0], gwf, mf6ep, "1B", "MODFLOW 6 PRT", colordest=colordest, starting=True, legend=True
        )
        plot_endpoints(
            fig, axes[1], gwf, mp7ep, "1B", "MODPATH 7", colordest=colordest, starting=True
        )
        
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-rel-capt.png")

    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8, 8))
        plt.suptitle(
            t="Example 1B release points, colored by travel time",
            fontsize=14,
            y=0.8
        )
        plt.tight_layout()
        axes = axes.flatten()
        plot_endpoints(fig, axes[0], gwf, mf6ep, "1B", "MODFLOW 6 PRT", starting=True, legend=True)
        plot_endpoints(fig, axes[1], gwf, mp7ep, "1B", "MODPATH 7", starting=True)
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-rel-tt.png")


def plot_termination(gwf, mf6ep, mp7ep):
    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8, 8))
        plt.suptitle(
            t="Example 1B terminating points, colored by destination",
            fontsize=14,
            y=0.8
        )
        axes = axes.flatten()
        plot_endpoints(
            fig,
            axes[0],
            gwf,
            mf6ep,
            "1B",
            "MODFLOW 6 PRT",
            colordest=colordest,
            starting=False,
            legend=True
        )
        plot_endpoints(
            fig, axes[1], gwf, mp7ep, "1B", "MODPATH 7", colordest=colordest, starting=False
        )
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-term-capt.png")

    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8, 8))
        plt.suptitle(
            t="Example 1B terminating points, colored by travel time",
            fontsize=14,
            y=0.8
        )
        axes = axes.flatten()
        plot_endpoints(fig, axes[0], gwf, mf6ep, "1B", "MODFLOW 6 PRT", starting=False, legend=True)
        plot_endpoints(fig, axes[1], gwf, mp7ep, "1B", "MODPATH 7", starting=False)
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-term-tt.png")


def plot_results(sim):
    mf6pathlines = load_mf6pathlines(mf6_ws / trackcsvfile_prt)
    mf6endpoints = load_mf6endpoints(mf6_ws / trackcsvfile_prt)
    mp7pathlines = load_mp7pathlines(
        flopy.utils.PathlineFile(mp7_ws / f"{mp7_name}.mppth")
    )
    mp7endpoints = load_mp7endpointdata(
        flopy.utils.EndpointFile(mp7_ws / f"{mp7_name}.mpend")
    )
    gwf = sim.get_model(gwf_name)
    plot_grid(gwf)
    plot_paths(gwf, mf6pathlines, mp7pathlines)
    plot_release(gwf, mf6endpoints, mp7endpoints)
    plot_termination(gwf, mf6endpoints, mp7endpoints)


# -

# ### Running the example
#
# Define a function to run the example scenarios and plot results.


def scenario(silent=True):
    sim, mp7 = build_model(sim_name)
    if write:
        write_model(sim, mp7, silent=silent)
    if run:
        run_model(sim, mp7, silent=silent)
    if plot:
        plot_results(sim)


# We are now ready to run the example problem. Subproblems 1A and 1B are solved by a single MODFLOW 6 run and a single MODPATH 7 run, so they are included under one "scenario". Each of the two subproblems is represented by its own particle release package (for MODFLOW 6) or particle group (for MODPATH 7).
scenario(silent=False)
