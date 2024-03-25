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
from matplotlib.lines import Line2D
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
nodes["well"] = [ncol * (nrow * k + i) + j]
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

    # Instantiate the MODFLOW 6 prt model input package.
    # Assign a different zone number to active cells, well cells, and river cells.
    # This makes it easier to determine where particles terminate.
    izone = np.zeros((nlay, nrow, ncol), dtype=int)
    for l, r, c in gwf.modelgrid.get_lrc(nodes["well"]):
        izone[l, r, c] = 1
    for l, r, c in gwf.modelgrid.get_lrc(nodes["river"]):
        izone[l, r, c] = 2
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity, izone=izone)

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
        zones=izone,
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


# Define a function to load pathline data from MODFLOW 6 PRT and MODPATH 7 pathline files.

# +


def get_pathlines(mf6_path, mp7_path):
    # load mf6 pathlines
    mf6pl = pd.read_csv(mf6_path)

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

    # determine which particles ended up in which capture area
    # (for now, use x coordinate of endpoint, todo: use izone)
    mf6pl["destzone"] = mf6pl[mf6pl.istatus > 1].izone
    mf6pl["dest"] = mf6pl.apply(
        lambda row: (
            "well" if row.destzone == 1 else "river" if row.destzone == 2 else pd.NA
        ),
        axis=1,
    )

    # load mp7 pathlines, letting flopy determine capture areas
    mp7plf = flopy.utils.PathlineFile(mp7_path)
    mp7pl_wel = pd.DataFrame(
        mp7plf.get_destination_pathline_data(nodes["well"], to_recarray=True)
    )
    mp7pl_wel["destzone"] = 1
    mp7pl_wel["dest"] = "well"
    mp7pl_riv = pd.DataFrame(
        mp7plf.get_destination_pathline_data(nodes["river"], to_recarray=True)
    )
    mp7pl_riv["destzone"] = 2
    mp7pl_riv["dest"] = "river"
    mp7pl = pd.concat([mp7pl_wel, mp7pl_riv])

    # index by particle group and particle ID
    mp7pl.set_index(["particlegroup", "sequencenumber"], drop=False, inplace=True)

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mp7pl["subprob"] = mp7pl.apply(
        lambda row: (
            "A" if row.particlegroup == 0 else "B" if row.particlegroup == 1 else pd.NA
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

    return mf6pl, mp7pl


# -

# ### Plotting results
#
# Define functions to plot model results.


# +
# Pathline and starting point colors by capture destination
colordest = {"well": "red", "river": "blue"}


def plot_lines(ax, gwf, data, **kwargs):
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(lw=0.5, alpha=0.5)
    for dest in ["well", "river"]:
        label = None
        if "colordest" in kwargs:
            color = kwargs["colordest"][dest]
            label = "captured by " + dest
        elif "color" in kwargs:
            color = kwargs["color"]
        else:
            color = "blue"
        mm.plot_pathline(
            data[data.dest == dest],
            layer="all",
            colors=[color],
            label=label,
            linewidth=0.5,
        )


def plot_points(ax, gwf, data, **kwargs):
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(lw=0.5, alpha=0.5)
    if "colordest" in kwargs:
        pts = []
        for dest in ["well", "river"]:
            color = kwargs["colordest"][dest]
            label = "captured by " + dest
            pdata = data[data.dest == dest]
            pts.append(
                ax.scatter(
                    pdata["x"],
                    pdata["y"],
                    s=3,
                    color=color,
                    label=label,
                )
            )
        return pts
    else:
        return ax.scatter(data["x"], data["y"], s=1, c=data["tt"])


def plot_grid(gwf, head=False, title=None, idx=0):
    with styles.USGSPlot():
        fig = plt.figure(figsize=(6, 6))
        fig.tight_layout()
        ax = fig.add_subplot(1, 1, 1, aspect="equal")
        pmv = flopy.plot.PlotMapView(gwf, ax=ax)
        if head:
            cb = pmv.plot_array(gwf.output.head().get_data(), alpha=0.25)
            cbar = plt.colorbar(cb, shrink=0.5, pad=0.1)
            cbar.ax.set_xlabel(r"Head ($m$)")
        pmv.plot_bc("WEL", plotAll=True)
        pmv.plot_bc("RIV", plotAll=True)
        pmv.plot_grid(lw=0.5, alpha=0.5)
        ax.legend(
            handles=[
                mpl.patches.Patch(color="red", label="Well"),
                mpl.patches.Patch(color="teal", label="River"),
            ],
            loc="upper left",
        )
        if title:
            styles.heading(ax, heading=f" {title}", idx=idx)
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-grid.png")


def plot_pathlines(gwf, mf6pl, mp7pl, idx):
    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(6, 6))
        styles.heading(
            axes[0][0], heading=" Pathlines, colored by destination", idx=idx
        )
        plot_lines(axes[0][0], gwf, mf6pl[mf6pl.subprob == "A"], colordest=colordest)
        plot_lines(axes[1][0], gwf, mf6pl[mf6pl.subprob == "B"], colordest=colordest)
        plot_lines(axes[0][1], gwf, mp7pl[mp7pl.subprob == "A"], colordest=colordest)
        plot_lines(axes[1][1], gwf, mp7pl[mp7pl.subprob == "B"], colordest=colordest)
        plt.subplots_adjust(hspace=0.1)
        axes[0][0].set_ylabel("1A")
        axes[1][0].set_ylabel("1B")
        axes[1][0].set_xlabel("MODFLOW 6 PRT")
        axes[1][1].set_xlabel("MODPATH 7")
        styles.graph_legend(
            handles=[
                Line2D([0], [0], color="red", lw=4, label="Well"),
                Line2D([0], [0], color="blue", lw=4, label="River"),
            ],
            bbox_to_anchor=(0.77, 0.09),
            bbox_transform=fig.transFigure,
            ncol=2,
        )
        plt.subplots_adjust(bottom=0.15)
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-paths.png")


def plot_release(gwf, mf6pl, mp7pl, idx, color="destination"):
    mf6sp = mf6pl[mf6pl.ireason == 0]  # release event
    mp7sp = mp7pl[mp7pl.time == 0]
    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(6, 6))
        axes = axes.flatten()
        styles.heading(axes[0], heading=f" Release points, colored by {color}", idx=idx)
        kwargs = {}
        if color == "destination":
            kwargs["colordest"] = colordest
        plot_points(axes[0], gwf, mf6sp[mf6sp.subprob == "A"], **kwargs)
        plot_points(axes[1], gwf, mp7sp[mp7sp.subprob == "A"], **kwargs)
        plot_points(axes[2], gwf, mf6sp[mf6sp.subprob == "B"], **kwargs)
        pts = plot_points(axes[3], gwf, mp7sp[mp7sp.subprob == "B"], **kwargs)
        plt.subplots_adjust(hspace=0.1)
        axes[0].set_ylabel("1A")
        axes[2].set_ylabel("1B")
        axes[2].set_xlabel("MODFLOW 6 PRT")
        axes[3].set_xlabel("MODPATH 7")
        if color == "destination":
            styles.graph_legend(
                handles=[
                    Line2D([0], [0], color="red", lw=4, label="Well"),
                    Line2D([0], [0], color="blue", lw=4, label="River"),
                ],
                bbox_to_anchor=(0.77, 0.09),
                bbox_transform=fig.transFigure,
                ncol=2,
            )
        else:
            cax = fig.add_axes([0.05, 0.06, 0.9, 0.01])
            cb = plt.colorbar(pts, cax=cax, orientation="horizontal")
            cb.set_label("Travel time")
        plt.subplots_adjust(bottom=0.15)
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-rel-{color}.png")


def plot_termination(gwf, mf6pl, mp7pl, idx, color="destination"):
    mf6ep = mf6pl[mf6pl.ireason == 3]  # termination event
    mp7ep = (
        mp7pl.sort_values("time")
        .groupby(level=["particlegroup", "sequencenumber"])
        .tail(1)
    )
    with styles.USGSPlot():
        fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(6, 6))
        axes = axes.flatten()
        styles.heading(
            axes[0], heading=f" Terminating points, colored by {color}", idx=idx
        )
        kwargs = {}
        if color == "destination":
            kwargs["colordest"] = colordest
        plot_points(axes[0], gwf, mf6ep[mf6ep.subprob == "A"], **kwargs)
        plot_points(axes[1], gwf, mp7ep[mp7ep.subprob == "A"], **kwargs)
        plot_points(axes[2], gwf, mf6ep[mf6ep.subprob == "B"], **kwargs)
        pts = plot_points(axes[3], gwf, mp7ep[mp7ep.subprob == "B"], **kwargs)
        axes[0].set_ylabel("1A")
        axes[2].set_ylabel("1B")
        axes[2].set_xlabel("MODFLOW 6 PRT")
        axes[3].set_xlabel("MODPATH 7")
        if color == "destination":
            styles.graph_legend(
                handles=[
                    Line2D([0], [0], color="red", lw=4, label="Well"),
                    Line2D([0], [0], color="blue", lw=4, label="River"),
                ],
                bbox_to_anchor=(0.77, 0.09),
                bbox_transform=fig.transFigure,
                ncol=2,
            )
        else:
            cax = fig.add_axes([0.05, 0.06, 0.9, 0.01])
            cb = plt.colorbar(pts, cax=cax, orientation="horizontal")
            cb.set_label("Travel time")
        plt.subplots_adjust(bottom=0.15)
        if plot_show:
            plt.show()
        if plot_save:
            fig.savefig(figs_path / f"{sim_name}-term-{color}.png")


def plot_all(sim):
    # get gwf model
    gwf = sim.get_model(gwf_name)

    # get pathlines
    mf6pathlines, mp7pathlines = get_pathlines(
        mf6_ws / trackcsvfile_prt,
        mp7_ws / f"{mp7_name}.mppth",
    )

    # plot the results
    plot_grid(gwf, title="Boundary Conditions", idx=0)
    plot_grid(gwf, head=True, title="Simulated Head", idx=1)
    plot_pathlines(gwf, mf6pathlines, mp7pathlines, idx=2)
    plot_release(gwf, mf6pathlines, mp7pathlines, idx=3, color="destination")
    plot_release(gwf, mf6pathlines, mp7pathlines, idx=4, color="travel-time")
    plot_termination(gwf, mf6pathlines, mp7pathlines, idx=5, color="destination")


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
        plot_all(sim)


# We are now ready to run the example problem. Subproblems 1A and 1B are solved by a single MODFLOW 6 run and a single MODPATH 7 run, so they are included under one "scenario". Each of the two subproblems is represented by its own particle release package (for MODFLOW 6) or particle group (for MODPATH 7).
scenario(silent=False)
