#
#
# # Original TWRI MODFLOW example
#
# This problem is described in McDonald and Harbaugh (1988) and duplicated in
# Harbaugh and McDonald (1996). This problem is also is distributed with
# MODFLOW-2005 (Harbaugh, 2005) and MODFLOW 6 (Langevin and others, 2017).
#
#

# Imports
import os
import sys
import matplotlib.pyplot as plt
import flopy

# Append to system path to include the common subdirectory
sys.path.append(os.path.join("..", "common"))

# import common functionality
import config
from figspecs import USGSFigure

# Set figure properties specific to the
figure_size = (3, 3)

# Base simulation and model name and workspace
ws = config.base_ws

# Scenario parameters - make sure there is at least one blank line before next item
parameters = {
    "ex-gwf-twri01": {"recharge": 3e-8, "evapotranspiration": 0.0},
    "ex-gwf-twri02": {"recharge": 2e-8, "squids": 15, "dogs": 0},
}

# Scenario parameter units - make sure there is at least one blank line before next item
# add parameter_units to add units to the scenario parameter table
parameter_units = {"recharge": "$ft/s$", "evapotranspiration": "$ft/s$",
                   "squids": "unitless", "dogs": "unitless"}

# ## Model units
length_units = "feet"
time_units = "seconds"

# Table
nper = 1  # Number of periods
nlay = 3  # Number of layers
ncol = 15  # Number of columns
nrow = 15  # Number of rows
delr = 5000.0  # Column width ($ft$)
delc = 5000.0  # Row width ($ft$)
top = 200.0  # Top of the model ($ft$)
botm_str = "-200.0, -300.0, -450.0"  # Layer bottom elevations ($ft$)
strt = 0.0  # Starting head ($ft$)
icelltype_str = "1, 0, 0"  # Cell conversion type
k11_str = "1.0e-3, 1e-4, 2.0e-4"  # Horizontal hydraulic conductivity ($ft/s$)
k33 = 2.0e-8  # Vertical hydraulic conductivity ($ft/s$)

# ## Static temporal data used by TDIS file
tdis_ds = ((8.640e04, 1, 1.000e00),)

# parse parameter strings into tuples
botm = tuple([float(value) for value in botm_str.split(",")])
k11 = tuple([float(value) for value in k11_str.split(",")])
icelltype = tuple([int(value) for value in icelltype_str.split(",")])

# ## Boundary Condition - CHD
# Constant head cells are specified on the west edge of the model
# in model layers 1 and 2 `(k, i, j)` = $(1 \rightarrow 2, 1 \rightarrow 15, 1)$
chd_spd = []
for k in range(2):
    chd_spd += [[(k, i, 0), 0.0] for i in range(nrow)]
chd_spd = {0: chd_spd}

# ## Boundary Condition - WEL
# Well boundary conditions
wel_spd = {
    0: [
        [(2, 4, 10), -5.0],
        [(1, 3, 5), -5.0],
        [(1, 5, 11), -5.0],
        [(0, 8, 7), -5.0],
        [(0, 8, 9), -5.0],
        [(0, 8, 11), -5.0],
        [(0, 8, 13), -5.0],
        [(0, 10, 7), -5.0],
        [(0, 10, 9), -5.0],
        [(0, 10, 11), -5.0],
        [(0, 10, 13), -5.0],
        [(0, 12, 7), -5.0],
        [(0, 12, 9), -5.0],
        [(0, 12, 11), -5.0],
        [(0, 12, 13), -5.0],
    ]
}

# ## Boundary Condition - DRN
# Drain boundary conditions
drn_spd = {
    0: [
        [(0, 7, 1), 0.0, 1.0e0],
        [(0, 7, 2), 0.0, 1.0e0],
        [(0, 7, 3), 10.0, 1.0e0],
        [(0, 7, 4), 20.0, 1.0e0],
        [(0, 7, 5), 30.0, 1.0e0],
        [(0, 7, 6), 50.0, 1.0e0],
        [(0, 7, 7), 70.0, 1.0e0],
        [(0, 7, 8), 90.0, 1.0e0],
        [(0, 7, 9), 100.0, 1.0e0],
    ]
}


# ## Function to build model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model
# recharge is the only variable
#
def build_model(sim_name, recharge=0.0, evapotranspiration=0.0, dogs=0, squids=0.0):
    if config.buildModel:
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
        )
        flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
        flopy.mf6.ModflowIms(sim)
        gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True)
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
            gwf, cvoptions="perched", perched=True, icelltype=icelltype, k=k11, k33=k33,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=strt)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=drn_spd)
        flopy.mf6.ModflowGwfwel(gwf, stress_period_data=wel_spd)
        flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
        head_filerecord = "{}.hds".format(sim_name)
        budget_filerecord = "{}.cbc".format(sim_name)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )
        return sim
    return None


# ## Function to write model files
def write_model(sim):
    if config.writeModel:
        sim.write_simulation(silent=True)


# ## Function to run the model
# True is returned if the model runs successfully
def run_model(sim):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=True)
        if not success:
            print(buff)
    return success


# ## Function to plot the model results
def plot_results(sim, idx):
    if config.plotModel:
        fs = USGSFigure(figure_type="map", verbose=False)
        sim_name = sim.name
        sim_ws = os.path.join(ws, sim_name)
        gwf = sim.get_model(sim_name)
        fpth = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, fpth)
        hobj = flopy.utils.HeadFile(fpth)

        # Create figure for scenario 1
        fig, axs = plt.subplots(1, 1, figsize=figure_size, dpi=300, tight_layout=True)
        fmp = flopy.plot.PlotMapView(model=gwf, ax=axs)
        fmp.plot_grid(lw=0.5)
        fmp.plot_array(hobj.get_data())
        title = "Recharge rate = " + "{} {}/{}".format(
            parameters[sim_name]["recharge"], length_units, time_units
        )
        letter = chr(ord("@") + idx + 1)
        fs.heading(letter=letter, heading=title)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}{}".format(sim_name, config.figure_ext)
            )
            fig.savefig(fpth)


# ## Function that wraps all of the steps for each scenario
#
# 1. build_model
# 2. write_model
# 3. run_model
# 4. plot_results
#
def scenario(idx):
    key = list(parameters.keys())[idx]
    parameter_dict = parameters[key]
    sim = build_model(key, **parameter_dict)

    write_model(sim)

    success = run_model(sim)

    if success:
        plot_results(sim, idx)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0)


def test_02():
    scenario(1)


# nosetest end

if __name__ == "__main__":
    # ## Scenario 1
    # This scenario uses the recharge rate defined in the original problem
    scenario(0)

    # ## Scenario 2
    # This uses a recharge rate less that the rate defined in the original problem
    scenario(1)
