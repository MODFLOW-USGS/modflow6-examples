# ## Stallman Problem
#
# Periodic heat boundary condition at surface
# Transient heat transfer problem in vertical
#

# ### Stallman Problem Setup

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
from analytical import Stallman

mf6exe = os.path.abspath(config.mf6_exe)

# Set figure properties specific to this problem

figure_size = (6, 8)

# Base simulation and model name and workspace

ws = config.base_ws
example_name = "ex-gwt-stallman"

# Model units

length_units = "m"
time_units = "seconds"

# Table of model parameters

nper = 600  # Number of periods
nstp = 6  # Number of time steps
perlen = 525600  # Simulation time length ($s$)
nlay = 120  # Number of layers
nrow = 1  # Number of rows
ncol = 1  # Number of columns
system_length = 60.0  # Length of system ($m$)
delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
delv_str = "ranges from 0.1 to 1"  # Layer thickness
top = 60.0  # Top of the model ($m$)
hydraulic_conductivity = 1.0e-4  # Hydraulic conductivity ($m s^{-1}$)
porosity = 0.35  # Porosity (unitless)
alphal = 0.0  # Longitudinal dispersivity ($m$)
alphat = 0.0  # Transverse dispersivity ($m$)
diffc = 1.02882E-06  # Diffusion coefficient ($m s^{-1}$)
T_az = 10  # Ambient temperature ($^o C$)
dT = 5  # Temperature variation ($^o C$)
bulk_dens = 2630  # Bulk density ($kg/m^3$)
kd = 0.000191663  # Distribution coefficient (unitless)

# Stress period input
per_data = []
for k in range(nper):
    per_data.append((perlen, nstp, 1.0))
per_mf6 = per_data

# Geometry input
tp = top
botm = []
for i in range(nlay):
    if i==0:botm.append(59.9)
    elif i==119:botm.append(0.0)
    else: botm.append(60-i*0.5)

# Head input
chd_data = {}
for k in range(nper):
    chd_data[k] = [[(0, 0, 0), 60.000000],[(119, 0, 0), 59.701801]]
chd_mf6 = chd_data

# Initial temperature input
strt_conc = T_az* np.ones((nlay, 1, 1), dtype=np.float32)

# Boundary temperature input
cnc_data = {}
for k in range(nper):
    cnc_temp = T_az+dT*np.sin(2*np.pi*k*perlen/365/86400)
    cnc_data[k] = [[(0, 0, 0), cnc_temp]]
cnc_mf6 = cnc_data

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-8, 1e-8, 0.97

# ### Functions to build, write, run, and plot models
#
# MODFLOW 6 flopy GWF simulation object (sim) is returned
#

def build_model(sim_folder):
    print("Building model...{}".format(sim_folder))
    name = "flow"
    sim_ws = os.path.join(ws, sim_folder)
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=sim_ws,
        exe_name=config.mf6_exe,
    )
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=per_mf6, time_units=time_units
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
        linear_acceleration="CG",
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

    flopy.mf6.ModflowGwfic(gwf, strt=top)

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_mf6,
    )

    head_filerecord = "{}.hds".format(name)
    budget_filerecord = "{}.bud".format(name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
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
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        sorption="linear",
        bulk_density=bulk_dens*(1-porosity),
        distcoef=kd,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=strt_conc)
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtdsp(
        gwt, xt3d_off=True, alh=alphal, ath1=alphat, diffc=diffc
    )
    flopy.mf6.ModflowGwtssm(gwt, sources=[[]])
    flopy.mf6.ModflowGwtcnc(
        gwt,
        stress_period_data=cnc_mf6,
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord="{}.cbc".format(gwt.name),
        concentration_filerecord="{}.ucn".format(gwt.name),
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "LAST")],
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
    print("Running model...")
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
    sim_name = example_name
    sim_ws = os.path.join(ws, sim_name)
    gwf = sim.get_model("flow")
    gwt = sim.get_model("trans")

    # create MODFLOW 6 head object
    cobj = gwt.output.concentration()
    times = cobj.get_times()
    times = np.array(times)

    time_in_pub = 284349600.0
    idx_conc = (np.abs(times - time_in_pub)).argmin()
    time_this_plot = times[idx_conc]
    conc = cobj.get_data(totim=time_this_plot)

    zconc = np.zeros(nlay)
    zbotm = np.zeros(nlay)
    for i in range(len(zconc)):
        zconc[i] = conc[i][0][0]
        zbotm[i] = -(60-botm[i])

    # Analytical solution - Stallman analysis
    tau = 365*86400
    t = 283824000.0
    c_w = 4174
    rho_w = 1000
    c_r = 800
    rho_r = bulk_dens
    c_rho = c_r*rho_r*(1-porosity) + c_w*rho_w*porosity
    darcy_flux = 5.00E-07
    ko = 1.503
    zanal = Stallman(T_az,dT,tau,t,c_rho,darcy_flux,ko,c_w,rho_w)

    # make conc figure
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)

    # configure plot and save
    ax.plot(zconc, zbotm, "k--", linewidth=0.5)
    ax.plot(zanal[:,1], zanal[:,0], "bo", mfc="none")
    ax.set_xlim(T_az-dT, T_az+dT)
    ax.set_ylim(-top, 0)
    ax.set_ylabel("Depth (m)")
    ax.set_xlabel("Temperature (deg C)")

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-conc{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


def plot_results(sim, idx):
    print("Plotting results...")
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
    sim = build_model(example_name)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(sim, idx)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Salt Lake Problem

    # Plot showing MODFLOW 6 results

    scenario(0)
