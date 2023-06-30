# ## LGRV example
#
# These are the models described in Vilhelmsen et al. (2012).  The parent
# model is 9 layers, the child model is 25 layers.
#

# ### LGRV Problem Setup
#
# Imports

import os
import sys

import flopy
import flopy.utils.lgrutil
import matplotlib.pyplot as plt
import numpy as np

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set default figure properties

figure_size = (5, 4)

# Base simulation and data workspace

ws = config.base_ws
data_ws = os.path.join(config.data_ws, "ex-gwf-lgrv")

# Model units

length_units = "meters"
time_units = "seconds"

# Scenario parameters

parameters = {
    "ex-gwf-lgrv-gr": {"configuration": "Refined"},
    "ex-gwf-lgrv-gc": {"configuration": "Coarse"},
    "ex-gwf-lgrv-lgr": {"configuration": "LGR"},
}

# Table LGRV Model Parameters

nper = 1  # Number of periods
nlay = 25  # Number of layers in refined model
nrow = 183  # Number of rows in refined model
ncol = 147  # Number of columns in refined model
nlaygc = 9  # Number of layers in coarsened model
nrowcg = 61  # Number of rows in coarsened model
ncolcg = 49  # Number of columns in coarsened model
delr = 35.0  # Column width ($m$) in refined model
delc = 25.0  # Row width ($m$) in refined model
delv = 5.0  # Layer thickness ($m$) in refined model
delrgc = 105.0  # Column width ($m$) in coarsened model
delcgc = 75.0  # Row width ($m$) in coarsened model
delvgc = 15.0  # Layer thickness ($m$) in coarsened model
top_str = "variable"  # Top of the model ($m$)
botm_str = "30 to -90"  # Layer bottom elevations ($m$)
icelltype = 0  # Cell conversion type
recharge = 1.1098e-09  # Recharge rate ($m/s$)
k11_str = "5.e-07, 1.e-06, 5.e-05"  # Horizontal hydraulic conductivity ($m/s$)

# Static temporal data used by TDIS file
# Simulation has 1 steady stress period (1 day)

perlen = [1.0]
nstp = [1]
tsmult = [1.0]
tdis_ds = list(zip(perlen, nstp, tsmult))

# load data files and process into arrays

fname = os.path.join(data_ws, "top.dat")
top = np.loadtxt(fname)
ikzone = np.empty((nlay, nrow, ncol), dtype=float)
for k in range(nlay):
    fname = os.path.join(data_ws, f"ikzone{k + 1}.dat")
    ikzone[k, :, :] = np.loadtxt(fname)
fname = os.path.join(data_ws, "riv.dat")
dt = [
    ("k", int),
    ("i", int),
    ("j", int),
    ("stage", float),
    ("conductance", float),
    ("rbot", float),
]
rivdat = np.loadtxt(fname, dtype=dt)
rivdat["k"] -= 1
rivdat["i"] -= 1
rivdat["j"] -= 1
riv_spd = [[(k, i, j), stage, cond, rbot] for k, i, j, stage, cond, rbot in rivdat]

botm = [30 - k * delv for k in range(nlay)]
botm = np.array(botm)
k11_values = [float(value) for value in k11_str.split(",")]
k11 = np.zeros((nlay, nrow, ncol), dtype=float)
for i, kval in enumerate(k11_values):
    k11 = np.where(ikzone == i + 1, kval, k11)

# Define model extent and child model extent

xmin = 0
xmax = ncol * delr
ymin = 0.0
ymax = nrow * delc
model_domain = [xmin, xmax, ymin, ymax]
child_domain = [
    xmin + 15 * 3 * delr,
    xmin + 41 * 3 * delr,
    ymax - 49 * 3 * delc,
    ymax - 19 * 3 * delc,
]

# Solver parameters

nouter = 50
ninner = 100
hclose = 1e-6
rclose = 100.0

# ### Functions to build, write, run, and plot the MODFLOW 6 LGRV model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def coarsen_shape(icoarsen, nrow, ncol):
    nrowc = int(np.ceil(nrow / icoarsen))
    ncolc = int(np.ceil(ncol / icoarsen))
    return (nrowc, ncolc)


def create_resampling_labels(a, icoarsen):
    nrow, ncol = a.shape
    labels = np.zeros((nrow, ncol), dtype=int)
    nodec = 0
    for ic in range(0, nrow, icoarsen):
        for jc in range(0, ncol, icoarsen):
            labels[ic : ic + icoarsen, jc : jc + icoarsen] = nodec
            nodec += 1
    return labels


def array_resampler(a, icoarsen, method):
    import scipy.ndimage as ndimage

    assert method in ["mean", "minimum", "maximum", "sum"]
    nrow, ncol = a.shape
    nrowc, ncolc = coarsen_shape(icoarsen, nrow, ncol)
    labels = create_resampling_labels(a, icoarsen)
    idx = np.array(range(nrowc * ncolc))
    if method == "mean":
        ar = ndimage.mean(a, labels=labels, index=idx)
    elif method == "minimum":
        ar = ndimage.minimum(a, labels=labels, index=idx)
    elif method == "maximum":
        ar = ndimage.maximum(a, labels=labels, index=idx)
    elif method == "sum":
        ar = ndimage.sum(a, labels=labels, index=idx)
    return ar.reshape((nrowc, ncolc))


def riv_resample(icoarsen, nrow, ncol, rivdat, idomain, rowcolspan):
    stage_grid = np.zeros((nrow, ncol), dtype=float)
    cond_grid = np.zeros((nrow, ncol), dtype=float)
    rbot_grid = np.zeros((nrow, ncol), dtype=float)
    count_grid = np.zeros((nrow, ncol), dtype=int)
    for k, i, j, stage, cond, rbot in rivdat:
        stage_grid[i, j] = stage
        cond_grid[i, j] = cond
        rbot_grid[i, j] = rbot
        count_grid[i, j] += 1
    stagec_grid = array_resampler(stage_grid, icoarsen, "sum")
    condc_grid = array_resampler(cond_grid, icoarsen, "sum")
    rbotc_grid = array_resampler(rbot_grid, icoarsen, "sum")
    countc_grid = array_resampler(count_grid, icoarsen, "sum")
    stagec_grid = np.divide(stagec_grid, countc_grid)
    rbotc_grid = np.divide(rbotc_grid, countc_grid)
    if rowcolspan is not None:
        istart, istop, jstart, jstop = rowcolspan
        stagec_grid = stagec_grid[istart:istop, jstart:jstop]
        condc_grid = condc_grid[istart:istop, jstart:jstop]
        rbotc_grid = rbotc_grid[istart:istop, jstart:jstop]
        countc_grid = countc_grid[istart:istop, jstart:jstop]
    rows, cols = np.where(condc_grid > 0.0)
    rivdatc = []
    for i, j in zip(rows, cols):
        k = 0
        if idomain[k, i, j] == 1:
            rivdatc.append(
                [
                    (k, i, j),
                    stagec_grid[i, j],
                    condc_grid[i, j],
                    rbotc_grid[i, j],
                ]
            )
    return rivdatc


def build_lgr_model(sim_name):
    sim_ws = os.path.join(ws, sim_name)
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6"
    )
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
    flopy.mf6.ModflowIms(
        sim,
        outer_maximum=nouter,
        outer_dvclose=hclose,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=f"{rclose} strict",
    )

    # parent model with coarse grid
    icoarsen = 3
    ncppl = [1, 3, 3, 3, 3, 3, 3, 3, 3]
    sim = build_parent_model(sim, sim_name, icoarsen=icoarsen, ncppl=ncppl)
    gwf = sim.get_model("parent")

    # child model with fine grid
    sim = build_child_model(sim, sim_name)
    gwfc = sim.get_model("child")

    # use flopy lgr utility to wire up connections between parent and child
    nlayp = len(ncppl)
    nrowp = gwf.dis.nrow.get_data()
    ncolp = gwf.dis.ncol.get_data()
    delrp = gwf.dis.delr.array
    delcp = gwf.dis.delc.array
    topp = gwf.dis.top.array
    botmp = gwf.dis.botm.array
    idomainp = gwf.dis.idomain.array
    lgr = flopy.utils.lgrutil.Lgr(
        nlayp,
        nrowp,
        ncolp,
        delrp,
        delcp,
        topp,
        botmp,
        idomainp,
        ncpp=icoarsen,
        ncppl=ncppl,
    )

    # swap out lgr child top and botm with
    topc = gwfc.dis.top.array
    botmc = gwfc.dis.botm.array
    lgr.top = topc
    lgr.botm = botmc
    exgdata = lgr.get_exchange_data(angldegx=True, cdist=True)
    flopy.mf6.ModflowGwfgwf(
        sim,
        nexg=len(exgdata),
        exgtype="GWF6-GWF6",
        exgmnamea="parent",
        exgmnameb="child",
        exchangedata=exgdata,
        auxiliary=["angldegx", "cdist"],
    )

    return sim


def build_parent_model(sim, sim_name, icoarsen, ncppl):
    xminp, xmaxp, yminp, ymaxp = model_domain
    xminc, xmaxc, yminc, ymaxc = child_domain
    delcp = delc * icoarsen
    delrp = delr * icoarsen
    istart = int((ymaxp - ymaxc) / delcp)
    istop = int((ymaxp - yminc) / delcp)
    jstart = int((xminc - xminp) / delrp)
    jstop = int((xmaxc - xminp) / delrp)
    nrowp, ncolp = coarsen_shape(icoarsen, nrow, ncol)
    nlayp = len(ncppl)
    idomain = np.ones((nlayp, nrowp, ncolp), dtype=int)
    idomain[:, istart:istop, jstart:jstop] = 0
    sim = build_model(
        sim_name,
        icoarsen=icoarsen,
        ncppl=ncppl,
        idomain=idomain,
        sim=sim,
        modelname="parent",
    )
    return sim


def build_child_model(sim, sim_name):
    icoarsen = 1
    xminp, xmaxp, yminp, ymaxp = model_domain
    xminc, xmaxc, yminc, ymaxc = child_domain
    delcp = delc * icoarsen
    delrp = delr * icoarsen
    istart = int((ymaxp - ymaxc) / delcp)
    istop = int((ymaxp - yminc) / delcp)
    jstart = int((xminc - xminp) / delrp)
    jstop = int((xmaxc - xminp) / delrp)
    nrowp, ncolp = coarsen_shape(icoarsen, nrow, ncol)
    sim = build_model(
        sim_name,
        rowcolspan=[istart, istop, jstart, jstop],
        sim=sim,
        modelname="child",
        xorigin=xminc,
        yorigin=yminc,
    )
    return sim


def build_model(
    sim_name,
    icoarsen=1,
    ncppl=None,
    rowcolspan=None,
    idomain=None,
    sim=None,
    modelname=None,
    xorigin=None,
    yorigin=None,
):
    if config.buildModel:
        if sim is None:
            sim_ws = os.path.join(ws, sim_name)
            sim = flopy.mf6.MFSimulation(
                sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6"
            )
            flopy.mf6.ModflowTdis(
                sim, nper=nper, perioddata=tdis_ds, time_units=time_units
            )
            flopy.mf6.ModflowIms(
                sim,
                outer_maximum=nouter,
                outer_dvclose=hclose,
                inner_maximum=ninner,
                inner_dvclose=hclose,
                rcloserecord=f"{rclose} strict",
            )
        if modelname is None:
            modelname = sim_name
        gwf = flopy.mf6.ModflowGwf(sim, modelname=modelname, save_flows=True)

        if ncppl is not None:
            nlayc = len(ncppl)
            layer_index = [ncppl[0] - 1]
            for iln in ncppl[1:]:
                last = layer_index[-1]
                layer_index.append(iln + last)
        else:
            nlayc = nlay
            layer_index = list(range(nlayc))
        nrowc, ncolc = coarsen_shape(icoarsen, nrow, ncol)
        delrc = delr * icoarsen
        delcc = delc * icoarsen
        topc = array_resampler(top, icoarsen, "mean")
        if rowcolspan is not None:
            istart, istop, jstart, jstop = rowcolspan
            nrowc = istop - istart
            ncolc = jstop - jstart
        else:
            istart = 0
            istop = nrow
            jstart = 0
            jstop = ncol
        if idomain is None:
            idomain = 1
        topc = topc[istart:istop, jstart:jstop]
        flopy.mf6.ModflowGwfdis(
            gwf,
            length_units=length_units,
            nlay=nlayc,
            nrow=nrowc,
            ncol=ncolc,
            delr=delrc,
            delc=delcc,
            top=topc,
            botm=botm[layer_index],
            idomain=idomain,
            xorigin=xorigin,
            yorigin=yorigin,
        )
        idomain = gwf.dis.idomain.array
        k11c = []
        for k in range(nlayc):
            ilay = layer_index[k]
            a = array_resampler(k11[ilay], icoarsen, "maximum")
            k11c.append(a[istart:istop, jstart:jstop])
        flopy.mf6.ModflowGwfnpf(
            gwf,
            k33overk=True,
            icelltype=icelltype,
            k=k11c,
            save_specific_discharge=True,
            k33=1.0,
        )
        strt = nlayc * [topc]
        flopy.mf6.ModflowGwfic(gwf, strt=strt)

        rivdatc = riv_resample(icoarsen, nrow, ncol, rivdat, idomain, rowcolspan)
        riv_spd = {0: rivdatc}
        flopy.mf6.ModflowGwfriv(
            gwf,
            stress_period_data=riv_spd,
            pname="RIV",
        )
        flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge, pname="RCH")
        head_filerecord = f"{modelname}.hds"
        budget_filerecord = f"{modelname}.cbc"
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )
        return sim
    return None


# Function to write MODFLOW 6 LGRV model files


def write_model(sim, silent=True):
    if config.writeModel:
        print(f"Writing simulation {sim.name}")
        sim.write_simulation(silent=silent)


# Function to run the LGRV model.
# True is returned if the model runs successfully
#


@config.timeit
def run_model(sim, silent=False):
    success = True
    if config.runModel:
        print(f"Running simulation {sim.name}")
        success, buff = sim.run_simulation(silent=silent, report=True)
        if not success:
            print(buff)
    return success


# Function to plot the LGRV model results.
#
def plot_grid(sim):
    print(f"Plotting grid for {sim.name}...")
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_ws = sim.simulation_data.mfpath.get_sim_path()
    sim_name = sim.name
    gwf = sim.get_model("parent")
    gwfc = None
    if "child" in list(sim.model_names):
        gwfc = sim.get_model("child")

    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
    # pmv.plot_grid()
    idomain = gwf.dis.idomain.array
    tp = np.ma.masked_where(idomain[0] == 0, gwf.dis.top.array)
    vmin = tp.min()
    vmax = tp.max()
    if gwfc is not None:
        tpc = gwfc.dis.top.array
        vmin = min(vmin, tpc.min())
        vmax = max(vmax, tpc.max())

    cb = pmv.plot_array(tp, cmap="jet", alpha=0.25, vmin=vmin, vmax=vmax)
    pmv.plot_bc(name="RIV")
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    cbar = plt.colorbar(cb, shrink=0.5)
    cbar.ax.set_xlabel(r"Top, ($m$)")
    if gwfc is not None:
        pmv = flopy.plot.PlotMapView(model=gwfc, ax=ax, layer=0)
        _ = pmv.plot_array(
            tpc,
            cmap="jet",
            alpha=0.25,
            masked_values=[1e30],
            vmin=vmin,
            vmax=vmax,
        )
        pmv.plot_bc(name="RIV")
    if gwfc is not None:
        xmin, xmax, ymin, ymax = child_domain
        ax.plot(
            [xmin, xmax, xmax, xmin, xmin],
            [ymin, ymin, ymax, ymax, ymin],
            "k--",
        )
    xmin, xmax, ymin, ymax = model_domain
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # save figure
    if config.plotSave:
        fpth = os.path.join("..", "figures", f"{sim_name}-grid{config.figure_ext}")
        fig.savefig(fpth)
    return


def plot_xsect(sim):
    print(f"Plotting cross section for {sim.name}...")
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_ws = sim.simulation_data.mfpath.get_sim_path()
    sim_name = sim.name
    gwf = sim.get_model("parent")

    fig = plt.figure(figsize=(5, 2.5))
    fig.tight_layout()

    ax = fig.add_subplot(1, 1, 1)
    irow, icol = gwf.modelgrid.intersect(3000.0, 3000.0)
    pmv = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"column": icol})
    pmv.plot_grid(linewidth=0.5)
    hyc = np.log(gwf.npf.k.array)
    cb = pmv.plot_array(hyc, cmap="jet", alpha=0.25)
    ax.set_xlabel("y position (m)")
    ax.set_ylabel("z position (m)")
    cbar = plt.colorbar(cb, shrink=0.5)
    cbar.ax.set_xlabel(r"K, ($m/s$)")

    # save figure
    if config.plotSave:
        fpth = os.path.join("..", "figures", f"{sim_name}-xsect{config.figure_ext}")
        fig.savefig(fpth)
    return


def plot_heads(sim):
    print(f"Plotting results for {sim.name} ...")
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_ws = sim.simulation_data.mfpath.get_sim_path()
    sim_name = sim.name
    gwf = sim.get_model("parent")
    modelname = gwf.name
    gwfc = None
    if "child" in list(sim.model_names):
        gwfc = sim.get_model("child")

    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    print("  Loading heads...")
    layer = 0
    head = gwf.output.head().get_data()
    head = np.ma.masked_where(head > 1e29, head)
    vmin = head[layer].min()
    vmax = head[layer].max()
    if gwfc is not None:
        headc = gwfc.output.head().get_data()
        vmin = min(vmin, headc.min())
        vmax = max(vmax, headc.max())

    print("  Making figure...")
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
    cb = pmv.plot_array(head, cmap="jet", masked_values=[1e30], vmin=vmin, vmax=vmax)
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    cbar = plt.colorbar(cb, shrink=0.5)
    cbar.ax.set_xlabel(r"Head, ($m$)")
    if gwfc is not None:
        pmv = flopy.plot.PlotMapView(model=gwfc, ax=ax, layer=0)
        cb = pmv.plot_array(
            headc, cmap="jet", masked_values=[1e30], vmin=vmin, vmax=vmax
        )
        xmin, xmax, ymin, ymax = child_domain
        ax.plot(
            [xmin, xmax, xmax, xmin, xmin],
            [ymin, ymin, ymax, ymax, ymin],
            "k--",
        )
    xmin, xmax, ymin, ymax = model_domain
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # save figure
    if config.plotSave:
        fpth = os.path.join("..", "figures", f"{sim_name}-head{config.figure_ext}")
        fig.savefig(fpth)
    return


def plot_results(sim, silent=True):
    if config.plotModel:
        plot_grid(sim)
        plot_xsect(sim)
        plot_heads(sim)
    return


# Function that wraps all of the steps for the LGRV model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(idx, silent=True):
    key = list(parameters.keys())[idx]
    params = parameters[key].copy()
    if params["configuration"] == "Refined":
        sim = build_model(key, modelname="parent")
    elif params["configuration"] == "Coarse":
        ncppl = [1, 3, 3, 3, 3, 3, 3, 3, 3]
        sim = build_model(key, icoarsen=3, ncppl=ncppl, modelname="parent")
    elif params["configuration"] == "LGR":
        sim = build_lgr_model(key)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(sim, silent=silent)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(0, silent=False)


def test_02():
    simulation(1, silent=False)


def test_03():
    simulation(2, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### LGRV Simulation
    #
    # Global Refined Model

    simulation(0)

    # Global Coarse Model

    simulation(1)

    # Locally Refined Grid Model

    simulation(2)
