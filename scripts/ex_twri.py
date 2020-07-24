#
#
# Describe the model...
#
#
#

# Imports
import os
import sys
import matplotlib.pyplot as plt
import flopy
import spnspecs as spn

# Append to system path
sys.path.append(os.path.join("..", "common"))

# import from common directory
import config

# Set figure properties
figure_size = (3, 3)
figure_ext = ".png"

# Base name and workspace
basename = "ex_twri"
ws = os.path.join("..", "examples")

# Scenario parameters
parameters = {"{}01".format(basename): {"recharge": 3e-8},
              "{}02".format(basename): {"recharge": 2e-8}}

# Model units
length_units = "feet"
time_units = "seconds"

# Table
nlay = 3  # Number of layers
ncol = 15  # Number of columns
nrow = 15  # Number of rows
delr = 5000.0  # Column width, in feet
delc = 5000.0  # Row width, in feet
top = 200.0  # Top of the model, in feet
bottom = "-200, -300, -450"  # Layer bottom elevations, in feet

# Table of variables
nper = 1  # Number of periods

botm = (-200.0, -300.0, -450.0)

# Static temporal data
tdis_ds = ((8.640e04, 1, 1.000e00),)

# Total simulation time
total_sim_time = 0.0
for n in range(nper):
    total_sim_time += tdis_ds[n][0]

# Initial Conditions
strt = 0.0  # table "Starting head (feet)"

# Hydraulic properties
icelltype = (1, 0, 0)
k11 = (1.0e-3, 1e-4, 2.0e-4)
k33 = 2.0e-8

# Boundary Condition - CHD
chd_spd = []
for k in range(2):
    chd_spd += [[(k, i, 0), 0.0] for i in range(nrow)]
chd_spd = {0: chd_spd}

# Boundary Condition - WEL
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
# Some text describing the _boundary_ $x = y$
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


# Routine to build model - recharge is the only variable
def build_model(sim_name, recharge=0.):
    sim_ws = os.path.join(ws, sim_name)
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name, sim_ws=sim_ws, exe_name=config.mf6_exe
    )
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
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
        gwf,
        cvoptions="perched",
        perched=True,
        icelltype=icelltype,
        k=k11,
        k33=k33,
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


def write_model(sim):
    if config.writeModel:
        sim.write_simulation(silent=True)


def run_model(sim):
    if config.runModel:
        success, buff = sim.run_simulation(silent=True)
        if not success:
            print(buff)
        return success


def plot_results(sim):
    if config.plotModel:
        spn.set_map_specifications()
        sim_name = sim.name
        sim_ws = os.path.join(ws, sim_name)
        gwf = sim.get_model(sim_name)
        fpth = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, fpth)
        hobj = flopy.utils.HeadFile(fpth)

        # Create figure for scenario 1
        fig, axs = plt.subplots(1, 1, figsize=figure_size,
                                dpi=300, tight_layout=True)
        fmp = flopy.plot.PlotMapView(model=gwf, ax=axs)
        fmp.plot_grid(lw=0.5)
        fmp.plot_array(hobj.get_data())
        title = "Recharge rate = " + \
                "{} {}/{}".format(parameters[sim_name]["recharge"],
                                  length_units, time_units)
        axs.set_title(title)

        # save figure
        if config.plotSave:
            fpth = os.path.join("..", "figures",
                                "{}{}".format(sim_name, figure_ext))
            fig.savefig(fpth)


def scenario(idx):
    key = list(parameters.keys())[idx]
    parameter_dict = parameters[key]
    sim = build_model(key, **parameter_dict)

    write_model(sim)

    success = run_model(sim)

    if success:
        plot_results(sim)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0)


def test_02():
    scenario(1)


# nosetest end

if __name__ == "__main__":
    # Scenario 1
    scenario(0)

    # Scenario 2
    scenario(1)
