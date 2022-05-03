# ## USG-ex1 with GWF-GWF Exchange and XT3D
#
# This example shows how the MODFLOW 6 GWF-GWF Exchange can be used to simulate
# a nested grid problem. The example corresponds to the first example
# described in the MODFLOW-USG documentation. Instead of the ghost node feature,
# we use the XT3D option in the NPF package to improve the accuracy at the
# interface between the models.
#
# The problem is run for three different scenarios:
#
# 1. without XT3D enabled in the NPF package
# 2. with XT3D enabled in both models
# 3. with XT3D enabled in both models and at the interface
# 4. with XT3D enabled _only_ at the interface between the models
#
# ### Setup
#
# Imports

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import flopy
from flopy.utils.lgrutil import Lgr

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

# Set default figure properties

figure_size = (6, 6)

# Base simulation and model name and workspace

ws = config.base_ws

# Model units

length_units = "meters"
time_units = "days"

# Scenario parameters

parameters = {
    "ex-gwf-u1gwfgwf-s1": {"xt3d": "none", },
    "ex-gwf-u1gwfgwf-s2": {"xt3d": "models", },
    "ex-gwf-u1gwfgwf-s3": {"xt3d": "models_exchange", },
    "ex-gwf-u1gwfgwf-s4": {"xt3d": "exchange", },
}

# Table with Model Parameters

nper = 1  # Number of periods
nlay = 1  # Number of layers
top = 0.0  # Top of the model ($m$)
botm = -100.0  # Layer bottom elevations ($m$)
strt = 0.0  # Starting head ($m$)
h_left = 1.0  # Constant head boundary LEFT ($m$)
h_right = 0.0  # Constant head boundary RIGHT ($m$)
icelltype = 0  # Cell conversion type
k11 = 1.0  # Horizontal hydraulic conductivity ($m/d$)

# Static temporal data used by TDIS file
# Simulation has 1 steady stress period (1 day)
# with 1 time step

perlen = [1.0]
nstp = [1]
tsmult = [1.0, 1.0, 1.0]
tdis_ds = list(zip(perlen, nstp, tsmult))

# Coarse model grid

nlay = 1
nrow = ncol = 7
delr = 100.0
delc = 100.0
tp = 0.0
bt = -100.0
idomain = np.ones((nlay, nrow, ncol))
idomain[:, 2:5, 2:5] = 0
gwfname_outer = "outer"

# Refined model grid

rfct = 3
nrow_inner = ncol_inner = 9
delr_inner = 100.0 / rfct
delc_inner = 100.0 / rfct
idomain_inner = np.ones((nlay, nrow_inner, ncol_inner))
xorigin = 200.0
yorigin = 200.0
gwfname_inner = "inner"

# Solver parameters

nouter = 50
ninner = 100
hclose = 1e-9
rclose = 1e-6


# ### Functions to build, write, run, and plot the model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(sim_name, xt3d):
    xt3d_global = (xt3d == "models") or (xt3d == "models_exchange")
    xt3d_exg = (xt3d == "models_exchange") or (xt3d == "exchange")

    if config.buildModel:
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
        )
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_ds, time_units=time_units
        )
        flopy.mf6.ModflowIms(
            sim,
            linear_acceleration="bicgstab",
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord="{} strict".format(rclose),
        )

        # The coarse, outer model
        gwf_outer = flopy.mf6.ModflowGwf(sim, modelname=gwfname_outer, save_flows=True)
        flopy.mf6.ModflowGwfdis(
            gwf_outer,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            idomain=idomain,
            top=top,
            botm=botm,
        )
        flopy.mf6.ModflowGwfnpf(
            gwf_outer,
            icelltype=icelltype,
            k=k11,
            save_specific_discharge=True,
            xt3doptions=xt3d_global,
        )
        flopy.mf6.ModflowGwfic(gwf_outer, strt=strt)

        chd_spd = []
        chd_spd += [[0, i, 1.0] for i in [0, 7, 14, 18, 22, 26, 33]]
        chd_spd = {0: chd_spd}

        # constant head boundary LEFT
        left_chd = [
            [(ilay, irow, 0), h_left]
            for ilay in range(nlay)
            for irow in range(nrow)
        ]
        chd_spd = {0: left_chd}
        flopy.mf6.ModflowGwfchd(
            gwf_outer,
            stress_period_data=chd_spd,
            pname="CHD-LEFT",
            filename="{}.left.chd".format(gwfname_outer),
        )

        # constant head boundary RIGHT
        right_chd = [
            [(ilay, irow, ncol - 1), h_right]
            for ilay in range(nlay)
            for irow in range(nrow)
        ]
        chd_spd = {0: right_chd}
        flopy.mf6.ModflowGwfchd(
            gwf_outer,
            stress_period_data=chd_spd,
            pname="CHD-RIGHT",
            filename="{}.right.chd".format(gwfname_outer),
        )

        head_filerecord = "{}.hds".format(gwfname_outer)
        budget_filerecord = "{}.cbc".format(gwfname_outer)
        flopy.mf6.ModflowGwfoc(
            gwf_outer,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        # the refined, inner model
        gwf_inner = flopy.mf6.ModflowGwf(sim, modelname=gwfname_inner, save_flows=True)
        flopy.mf6.ModflowGwfdis(
            gwf_inner,
            nlay=nlay,
            nrow=nrow_inner,
            ncol=ncol_inner,
            delr=delr_inner,
            delc=delc_inner,
            top=top,
            botm=botm,
            xorigin=xorigin,
            yorigin=yorigin,
            length_units=length_units,
        )
        flopy.mf6.ModflowGwfic(gwf_inner, strt=strt)
        flopy.mf6.ModflowGwfnpf(
            gwf_inner,
            save_specific_discharge=True,
            xt3doptions=xt3d_global,
            save_flows=True,
            icelltype=icelltype,
            k=k11,
        )

        head_filerecord = "{}.hds".format(gwfname_inner)
        budget_filerecord = "{}.cbc".format(gwfname_inner)
        flopy.mf6.ModflowGwfoc(
            gwf_inner,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        # Use Lgr to get the exchange data
        nrowp = gwf_outer.dis.nrow.get_data()
        ncolp = gwf_outer.dis.ncol.get_data()
        delrp = gwf_outer.dis.delr.array
        delcp = gwf_outer.dis.delc.array
        topp = gwf_outer.dis.top.array
        botmp = gwf_outer.dis.botm.array
        idomainp = gwf_outer.dis.idomain.array

        lgr = Lgr(
            nlay,
            nrowp,
            ncolp,
            delrp,
            delcp,
            topp,
            botmp,
            idomainp,
            ncpp=rfct,
            ncppl=1,
        )

        exgdata = lgr.get_exchange_data(angldegx=True, cdist=True)
        for exg in exgdata:
            l = exg
            angle = l[-2]
            if angle == 0:
                bname = "left"
            elif angle == 90.0:
                bname = "bottom"
            elif angle == 180.0:
                bname = "right"
            elif angle == 270.0:
                bname = "top"
            l.append(bname)

        # group exchanges based on boundname
        exgdata.sort(key=lambda x: x[-3])

        flopy.mf6.ModflowGwfgwf(
            sim,
            exgtype="GWF6-GWF6",
            nexg=len(exgdata),
            exgmnamea=gwfname_outer,
            exgmnameb=gwfname_inner,
            exchangedata=exgdata,
            xt3d=xt3d_exg,
            print_input=True,
            print_flows=True,
            save_flows=True,
            boundnames=True,
            auxiliary=["ANGLDEGX", "CDIST"],
        )

        return sim
    return None


# Function to write model files

def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the model.
# True is returned if the model runs successfully
#

@config.timeit
def run_model(sim, silent=False):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent, report=True)
        if not success:
            print(buff)
    return success


# Functions to plot model results.
#

def plot_grid(idx, sim):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = list(parameters.keys())[idx]
    gwf_outer = sim.get_model(gwfname_outer)
    gwf_inner = sim.get_model(gwfname_inner)

    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf_outer, ax=ax, layer=0)
    pmv_inner = flopy.plot.PlotMapView(model=gwf_inner, ax=ax, layer=0)

    pmv.plot_grid()
    pmv_inner.plot_grid()

    pmv.plot_bc(name="CHD-LEFT", alpha=0.75)
    pmv.plot_bc(name="CHD-RIGHT", alpha=0.75)
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-grid{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


def plot_head(idx, sim):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = list(parameters.keys())[idx]
    gwf_outer = sim.get_model(gwfname_outer)
    gwf_inner = sim.get_model(gwfname_inner)

    fig = plt.figure(figsize=(15, 10))
    fig.tight_layout()

    head = gwf_outer.output.head().get_data()[0]
    head_inner = gwf_inner.output.head().get_data()[0]
    head[head == 1e+30] = np.nan
    head_inner[head_inner == 1e+30] = np.nan

    # create MODFLOW 6 cell-by-cell budget objects
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
        gwf_outer.output.budget().get_data(text="DATA-SPDIS", totim=1.0)[0],
        gwf_outer,
    )
    qx_inner, qy_inner, qz_inner = flopy.utils.postprocessing.get_specific_discharge(
        gwf_inner.output.budget().get_data(text="DATA-SPDIS", totim=1.0)[0],
        gwf_inner,
    )

    # create plot with head values and spdis
    ax = fig.add_subplot(1, 2, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf_outer, ax=ax, layer=0)
    pmv_inner = flopy.plot.PlotMapView(model=gwf_inner, ax=ax, layer=0, extent=pmv.extent)
    cb = pmv.plot_array(head, cmap="jet", vmin=0.0, vmax=1.0)
    cb = pmv_inner.plot_array(head_inner, cmap="jet", vmin=0.0, vmax=1.0)
    pmv.plot_grid()
    pmv_inner.plot_grid()
    pmv.plot_vector(
        qx, qy, normalize=False, color="0.75",
    )
    pmv_inner.plot_vector(
        qx_inner, qy_inner, normalize=False, color="0.75",
    )
    cbar = plt.colorbar(cb, shrink=0.25)
    cbar.ax.set_xlabel(r"Head, ($m$)")
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    fs.heading(ax, letter="A", heading="Simulated Head")

    # create plot with error in head
    ax = fig.add_subplot(1, 2, 2, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf_outer, ax=ax, layer=0)
    pmv_inner = flopy.plot.PlotMapView(model=gwf_inner, ax=ax, layer=0, extent=pmv.extent)
    pmv.plot_grid()
    pmv_inner.plot_grid()
    x = np.array(gwf_outer.modelgrid.xcellcenters) - 50.0
    x_inner = np.array(gwf_inner.modelgrid.xcellcenters) - 50.0
    slp = (h_left - h_right) / (50.0 - 650.0)
    head_exact = slp * x + h_left
    head_exact_inner = slp * x_inner + h_left
    err = head - head_exact
    err_inner = head_inner - head_exact_inner
    vmin = min(np.nanmin(err), np.nanmin(err_inner))
    vmax = min(np.nanmax(err), np.nanmax(err_inner))
    cb = pmv.plot_array(err, cmap="jet", vmin=vmin, vmax=vmax)
    cb = pmv_inner.plot_array(err_inner, cmap="jet", vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(cb, shrink=0.25)
    cbar.ax.set_xlabel(r"Error, ($m$)")
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("y position (m)")
    fs.heading(ax, letter="B", heading="Error")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-head{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


def plot_results(idx, sim, silent=True):
    if config.plotModel:
        if idx == 0:
            plot_grid(idx, sim)
        plot_head(idx, sim)
    return


# Function that wraps all of the steps for the FHB model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(idx, silent=True):
    key = list(parameters.keys())[idx]
    params = parameters[key].copy()
    sim = build_model(key, **params)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(idx, sim, silent=silent)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(0, silent=False)


def test_02():
    simulation(1, silent=False)


def test_03():
    simulation(2, silent=False)


def test_04():
    simulation(3, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### USG-ex1 GWF-GWF Exchange Simulation
    #
    # Simulated heads without XT3D.

    simulation(0)

    # Simulated heads with XT3D enabled globally, but not at the exchange

    simulation(1)

    # Simulated heads with XT3D enabled globally

    simulation(2)

    # Simulated heads with XT3D enabled _only_ at the model interface.

    simulation(3)
