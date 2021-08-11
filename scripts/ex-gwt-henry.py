# ## Henry Problem
#
# Classic saltwater intrusion
#
#


# ### Henry Problem Setup

# Imports

import os
import sys
import matplotlib.pyplot as plt
import flopy
import numpy as np

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# Import common functionality

import config
from figspecs import USGSFigure

mf6exe = os.path.abspath(config.mf6_exe)
exe_name_mf = config.mf2005_exe
exe_name_mt = config.mt3dms_exe

# Set figure properties specific to this problem

figure_size = (6, 4)

# Base simulation and model name and workspace

ws = config.base_ws

# Scenario parameters - make sure there is at least one blank line before next item

parameters = {
    "ex-gwt-henry-a": {"inflow": 5.7024,},
    "ex-gwt-henry-b": {"inflow": 2.851,},
}

# Scenario parameter units - make sure there is at least one blank line before next item
# add parameter_units to add units to the scenario parameter table

parameter_units = {
    "inflow": "$m^3/d$",
}

# Model units

length_units = "cm"
time_units = "seconds"

# Table of model parameters

nper = 1  # Number of periods
nstp = 500  # Number of time steps
perlen = 0.5  # Simulation time length ($d$)
nlay = 40  # Number of layers
nrow = 1  # Number of rows
ncol = 80  # Number of columns
system_length = 2.0  # Length of system ($m$)
delr = 0.025  # Column width ($m$)
delc = 1.0  # Row width ($m$)
delv = 0.025  # Layer thickness
top = 1.0  # Top of the model ($m$)
hydraulic_conductivity = 864.0  # Hydraulic conductivity ($m d^{-1}$)
initial_concentration = 35.0  # Initial concentration (unitless)
porosity = 0.35  # porosity (unitless)
diffusion_coefficient = 0.57024  # diffusion coefficient ($m^2/d$)

botm = [top - k * delv for k in range(1, nlay + 1)]

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-6, 0.97


# ### Functions to build, write, run, and plot models
#
# MODFLOW 6 flopy GWF simulation object (sim) is returned
#


def build_model(sim_folder, inflow):
    print("Building model...{}".format(sim_folder))
    name = "flow"
    sim_ws = os.path.join(ws, sim_folder)
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=sim_ws, exe_name=config.mf6_exe
    )
    tdis_ds = ((perlen, nstp, 1.0),)
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )
    sim.register_ims_package(ims, [gwf.name])
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
        icelltype=0,
        k=hydraulic_conductivity,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=initial_concentration)
    pd = [(0, 0.7, 0.0, "trans", "concentration")]
    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)
    ghbcond = hydraulic_conductivity * delv * delc / (0.5 * delr)
    ghbspd = [[(k, 0, ncol - 1), top, ghbcond, 35.0] for k in range(nlay)]
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
        pname="GHB-1",
        auxiliary="CONCENTRATION",
    )

    welspd = [[(k, 0, 0), inflow / nlay, 0.0] for k in range(nlay)]
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=welspd,
        pname="WEL-1",
        auxiliary="CONCENTRATION",
    )
    head_filerecord = "{}.hds".format(name)
    budget_filerecord = "{}.bud".format(name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwt.name),
    )
    sim.register_ims_package(imsgwt, [gwt.name])
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
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtic(gwt, strt=initial_concentration)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=diffusion_coefficient)
    sourcerecarray = [
        ("GHB-1", "AUX", "CONCENTRATION"),
        ("WEL-1", "AUX", "CONCENTRATION"),
    ]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord="{}.cbc".format(gwt.name),
        concentration_filerecord="{}.ucn".format(gwt.name),
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )
    return sim


# Function to write model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)
    return


# Function to run the model
# True is returned if the model runs successfully


@config.timeit
def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success = False
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)
    return success


# Function to plot the model results


def plot_conc(sim, idx):
    fs = USGSFigure(figure_type="map", verbose=False)
    sim_name = list(parameters.keys())[idx]
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model("flow")
    gwt = sim.get_model("trans")

    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    # get MODFLOW 6 concentration
    conc = gwt.output.concentration().get_data()

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 0})
    pxs.plot_array(conc, cmap="jet")
    levels = [35 * f for f in [0.01, 0.1, 0.5, 0.9, 0.99]]
    cs = pxs.contour_array(
        conc, levels=levels, colors="w", linewidths=1.0, linestyles="-"
    )
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("z position (m)")
    plt.clabel(cs, fmt="%4.2f", fontsize=5)

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-conc{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


def plot_results(sim, idx):
    if config.plotModel:
        plot_conc(sim, idx)
    return


# Function that wraps all of the steps for each scenario
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def scenario(idx, silent=True):
    key = list(parameters.keys())[idx]
    parameter_dict = parameters[key]
    sim = build_model(key, **parameter_dict)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(sim, idx)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0, silent=False)


def test_02():
    scenario(1, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Henry Problem

    # Scenario 1 - Classic henry problem

    scenario(0)

    # Scenario 2 - Modified Henry problem with half the inflow rate

    scenario(1)
