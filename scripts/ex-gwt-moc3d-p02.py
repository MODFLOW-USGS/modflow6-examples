# ## Three-Dimensional Steady Flow with Transport
#
# MOC3D Problem 2
#
#


# ### Three-Dimensional Steady Flow with Transport Problem Setup

# Imports

import os
import sys

import flopy
import matplotlib.pyplot as plt
import numpy as np

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# Import common functionality

import analytical
import config
from figspecs import USGSFigure

mf6exe = config.mf6_exe
exe_name_mf = config.mf2005_exe
exe_name_mt = config.mt3dms_exe

# Set figure properties specific to this problem

figure_size = (6, 4)

# Base simulation and model name and workspace

ws = config.base_ws
example_name = "ex-gwt-moc3d-p02"

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


def build_mf6gwf(sim_folder):
    print(f"Building mf6gwf model...{sim_folder}")
    name = "flow"
    sim_ws = os.path.join(ws, sim_folder, "mf6gwf")
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=sim_ws, exe_name=config.mf6_exe)
    tdis_ds = ((total_time, 1, 1.0),)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
    flopy.mf6.ModflowIms(sim, print_option="summary", inner_maximum=300)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
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
            rec = [(k, i, ncol - 1), 0.0]
            chdspd.append(rec)
            rec = [(k, i, 0), specific_discharge * delc * delv]
            welspd.append(rec)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=welspd)
    head_filerecord = f"{name}.hds"
    budget_filerecord = f"{name}.bud"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


# MODFLOW 6 flopy GWF simulation object (sim) is returned


def build_mf6gwt(sim_folder):
    print(f"Building mf6gwt model...{sim_folder}")
    name = "trans"
    sim_ws = os.path.join(ws, sim_folder, "mf6gwt")
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=sim_ws, exe_name=config.mf6_exe)
    tdis_ds = ((total_time, 400, 1.0),)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
    flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab")
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwtdis(
        gwt,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        alh=alpha_l,
        ath1=alpha_th,
        ath2=alpha_tv,
    )
    pd = [
        ("GWFHEAD", f"../mf6gwf/flow.hds", None),
        ("GWFBUDGET", "../mf6gwf/flow.bud", None),
    ]
    flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd)
    sourcerecarray = [[]]
    srcspd = [[source_location0, solute_mass_flux]]
    flopy.mf6.ModflowGwtsrc(gwt, stress_period_data=srcspd)
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    obs_data = {
        f"{name}.obs.csv": [
            ("SOURCELOC", "CONCENTRATION", source_location0),
        ],
    }
    obs_package = flopy.mf6.ModflowUtlobs(
        gwt, digits=10, print_input=True, continuous=obs_data
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
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


@config.timeit
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


def plot_results(sims):
    if config.plotModel:
        print("Plotting model results...")
        sim_mf6gwf, sim_mf6gwt = sims
        fs = USGSFigure(figure_type="map", verbose=False)

        conc = sim_mf6gwt.trans.output.concentration().get_data()

        fig, axs = plt.subplots(1, 1, figsize=figure_size, dpi=300, tight_layout=True)

        gwt = sim_mf6gwt.trans
        pmv = flopy.plot.PlotMapView(model=gwt, ax=axs)
        levels = [1, 3, 10, 30, 100, 300]
        cs1 = plot_analytical(axs, levels)
        cs2 = pmv.contour_array(conc, colors="blue", linestyles="--", levels=levels)
        axs.set_xlabel("x position (m)")
        axs.set_ylabel("y position (m)")
        axs.set_aspect(4.0)

        labels = ["Analytical", "MODFLOW 6"]
        lines = [cs1.collections[0], cs2.collections[0]]
        axs.legend(lines, labels, loc="upper left")

        # save figure
        if config.plotSave:
            sim_ws = sim_mf6gwt.simulation_data.mfpath.get_sim_path()
            sim_folder = os.path.split(sim_ws)[0]
            sim_folder = os.path.basename(sim_folder)
            fname = f"{sim_folder}-map{config.figure_ext}"
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
        plot_results(sims)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Model

    # Model run

    scenario(0)
