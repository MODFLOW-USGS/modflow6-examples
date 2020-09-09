# ## Original TWRI MODFLOW example
#
# This problem is described in McDonald and Harbaugh (1988) and duplicated in
# Harbaugh and McDonald (1996). This problem is also is distributed with
# MODFLOW-2005 (Harbaugh, 2005) and MODFLOW 6 (Langevin and others, 2017).
#

# ### TWRI Problem Setup
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

figure_size = (6, 6)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-twri01"

# Scenario parameter units - make sure there is at least one blank line before next item
# add parameter_units to add units to the scenario parameter table

parameter_units = {"recharge": "$ft/s$"}

# Model units

length_units = "feet"
time_units = "seconds"

# Table TWRI Model Parameters

nper = 1  # Number of periods
nlay = 5  # Number of layers
ncol = 15  # Number of columns
nrow = 15  # Number of rows
delr = 5000.0  # Column width ($ft$)
delc = 5000.0  # Row width ($ft$)
top = 200.0  # Top of the model ($ft$)
botm_str = "-150.0, -200.0, -300.0, -350.0, -450.0"  # Layer bottom elevations ($ft$)
strt = 0.0  # Starting head ($ft$)
icelltype_str = "1, 0, 0, 0, 0"  # Cell conversion type
k11_str = "1.0e-3, 1.0e-8, 1.0e-4, 5.0e-7, 2.0e-4"  # Horizontal hydraulic conductivity ($ft/s$)
k33_str = "1.0e-3, 1.0e-8, 1.0e-4, 5.0e-7, 2.0e-4"  # Vertical hydraulic conductivity ($ft/s$)
recharge = 3e-8  # Recharge rate ($ft/s$)

# Static temporal data used by TDIS file

perlen = 8.640e04
nstp = 1
tsmult = 1.
tdis_ds = ((perlen, nstp, tsmult),)

# parse parameter strings into tuples

botm = [float(value) for value in botm_str.split(",")]
k11 = [float(value) for value in k11_str.split(",")]
k33 = [float(value) for value in k33_str.split(",")]
icelltype = [int(value) for value in icelltype_str.split(",")]

# ### Create TWRI Model Boundary Conditions
#
# Constant head cells are specified on the west edge of the model
# in model layers 1 and 2 `(k, i, j)` = $(1 \rightarrow 2, 1 \rightarrow 15, 1)$
#

chd_spd = []
for k in (0, 1, 2, 3):
    chd_spd += [[k, i, 0, 0.0] for i in range(nrow)]
chd_spd = {0: chd_spd}

# Constant head cells for MODFLOW-2005 model

chd_spd0 = []
for k in (0, 1):
    chd_spd0 += [[k, i, 0, 0.0, 0.0] for i in range(nrow)]
chd_spd0 = {0: chd_spd0}

# Well boundary conditions

wel_spd = {
    0: [
        [4, 4, 10, -5.0],
        [2, 3, 5, -5.0],
        [2, 5, 11, -5.0],
        [0, 8, 7, -5.0],
        [0, 8, 9, -5.0],
        [0, 8, 11, -5.0],
        [0, 8, 13, -5.0],
        [0, 10, 7, -5.0],
        [0, 10, 9, -5.0],
        [0, 10, 11, -5.0],
        [0, 10, 13, -5.0],
        [0, 12, 7, -5.0],
        [0, 12, 9, -5.0],
        [0, 12, 11, -5.0],
        [0, 12, 13, -5.0],
    ]
}

# Well boundary conditions for MODFLOW-2005 model

wel_spd0 = []
layer_map = {0: 0, 2: 1, 4: 2}
for (k, i, j, q) in wel_spd[0]:
    kk = layer_map[k]
    wel_spd0.append([kk, i, j, q])
wel_spd0 = {0: wel_spd0}


# Drain boundary conditions

drn_spd = {
    0: [
        [0, 7, 1, 0.0, 1.0e0],
        [0, 7, 2, 0.0, 1.0e0],
        [0, 7, 3, 10.0, 1.0e0],
        [0, 7, 4, 20.0, 1.0e0],
        [0, 7, 5, 30.0, 1.0e0],
        [0, 7, 6, 50.0, 1.0e0],
        [0, 7, 7, 70.0, 1.0e0],
        [0, 7, 8, 90.0, 1.0e0],
        [0, 7, 9, 100.0, 1.0e0],
    ]
}


# ### Functions to build, write, run, and plot the MODFLOW 6 TWRI model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model

def build_model():
    if config.buildModel:
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
            save_specific_discharge=True,
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

# MODFLOW-2005 model object (mf) is returned if building the model

def build_mf5model():
    if config.buildModel:
        sim_ws = os.path.join(ws, sim_name, "mf2005")
        mf = flopy.modflow.Modflow(modelname=sim_name, model_ws=sim_ws,
                                   exe_name=config.mf2005_exe)
        flopy.modflow.ModflowDis(mf, nlay=3, nrow=nrow, ncol=ncol,
                                 delr=delr, delc=delc, laycbd=[1, 1, 0],
                                 top=top, botm=botm,
                                 nper=1, perlen=perlen, nstp=nstp, tsmult=tsmult)
        flopy.modflow.ModflowBas(mf, strt=strt)
        flopy.modflow.ModflowLpf(mf,
                                 laytyp=[1, 0, 0],
                                 hk=[k11[0], k11[2], k11[4]],
                                 vka=[k11[0], k11[2], k11[4]],
                                 vkcb=[k11[1], k11[3], 0],
                                 ss=0, sy=0.,
                                 )
        flopy.modflow.ModflowChd(mf, stress_period_data=chd_spd0)
        flopy.modflow.ModflowDrn(mf, stress_period_data=drn_spd)
        flopy.modflow.ModflowWel(mf, stress_period_data=wel_spd0)
        flopy.modflow.ModflowRch(mf, rech=recharge)
        flopy.modflow.ModflowPcg(mf)
        oc = flopy.modflow.ModflowOc(mf,
                                     stress_period_data={(0,0):
                                                             ['save head',
                                                              'save budget']})
        oc.reset_budgetunit()
        return mf
    return None


# Function to write MODFLOW 6 TWRI model files


def write_model(sim, mf, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)
        mf.write_input()


# Function to run the TWRI model.
# True is returned if the model runs successfully
#


def run_model(sim, mf, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)
        else:
            success, buff = mf.run_model(silent=silent)
            if not success:
                print(buff)

    return success


# Function to plot the TWRI model results.
#


def plot_results(sim, mf):
    if config.plotModel:
        fs = USGSFigure(figure_type="map", verbose=False)
        sim_ws = os.path.join(ws, sim_name)
        gwf = sim.get_model(sim_name)

        # create MODFLOW 6 head object
        file_name = gwf.oc.head_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        hobj = flopy.utils.HeadFile(fpth)

        # create MODFLOW-2005 head object
        fpth = os.path.join(sim_ws, "mf2005", file_name)
        hobj0 = flopy.utils.HeadFile(fpth)

        # create MODFLOW 6 cell-by-cell budget object
        file_name = gwf.oc.budget_filerecord.get_data()[0][0]
        fpth = os.path.join(sim_ws, file_name)
        cobj = flopy.utils.CellBudgetFile(fpth, precision="double")

        # create MODFLOW-2005 cell-by-cell budget object
        fpth = os.path.join(sim_ws, "mf2005", file_name)
        cobj0 = flopy.utils.CellBudgetFile(fpth, precision="single")

        # extract heads
        head = hobj.get_data()
        head0 = hobj0.get_data()
        vmin, vmax = -25, 100

        # extract specific discharge
        spdis = cobj.get_data(text="DATA-SPDIS", kstpkper=(0, 0))
        frf = cobj0.get_data(text="FLOW RIGHT FACE", kstpkper=(0, 0))[0]
        fff = cobj0.get_data(text="FLOW FRONT FACE", kstpkper=(0, 0))[0]
        flf = cobj0.get_data(text="FLOW LOWER FACE", kstpkper=(0, 0))[0]

        # modflow 6 layers to extract
        layers_mf6 = [0, 2, 4]
        titles = ["Unconfined aquifer",
                  "Middle aquifer",
                  "Lower aquifer"]

        # Create figure for simulation
        extents = (0, ncol * delc, 0, nrow * delr)
        fig, axes = plt.subplots(
            3, 3, figsize=figure_size, dpi=300, constrained_layout=True, sharey=True,
        )
        for ax in axes.flatten():
            ax.set_aspect("equal")
            ax.set_xlim(extents[:2])
            ax.set_ylim(extents[2:])

        for idx, ax in enumerate(axes.flatten()[:3]):
            k = layers_mf6[idx]
            fmp = flopy.plot.PlotMapView(
                model=gwf, ax=ax, layer=k, extent=extents
            )
            ax.get_xaxis().set_ticks([])
            fmp.plot_grid(lw=0.5)
            plot_obj = fmp.plot_array(head, vmin=vmin, vmax=vmax)
            fmp.plot_bc("DRN", color="green")
            fmp.plot_bc("WEL", color="0.5")
            cv = fmp.contour_array(
                head,
                levels=[-25, 0, 25, 75, 100],
                linewidths=0.5,
                colors="black",
            )
            plt.clabel(cv, fmt="%1.0f", fontsize=9)
            fmp.plot_specific_discharge(spdis, normalize=True, color="0.75")
            title = titles[idx]
            letter = chr(ord("@") + idx + 1)
            fs.heading(letter=letter, heading=title, ax=ax)

        for idx, ax in enumerate(axes.flatten()[3:6]):
            fmp = flopy.plot.PlotMapView(
                model=mf, ax=ax, layer=idx, extent=extents
            )
            fmp.plot_grid(lw=0.5)
            plot_obj = fmp.plot_array(head0, vmin=vmin, vmax=vmax)
            fmp.plot_bc("DRN", color="green")
            fmp.plot_bc("WEL", color="0.5")
            cv = fmp.contour_array(
                head0,
                levels=[-25, 0, 25, 75, 100],
                linewidths=0.5,
                colors="black",
            )
            plt.clabel(cv, fmt="%1.0f", fontsize=9)
            fmp.plot_discharge(frf, fff, flf, head=head0, normalize=True, color="0.75")
            title = titles[idx]
            letter = chr(ord("@") + idx + 4)
            fs.heading(letter=letter, heading=title, ax=ax)

        # create legend
        ax = plt.subplot(313)
        ax.set_xlim(extents[:2])
        ax.set_ylim(extents[2:])

        # ax = axes.flatten()[-2]
        ax.axis("off")
        ax.plot(
            -10000,
            -10000,
            lw=0,
            marker="s",
            ms=10,
            mfc="green",
            mec="green",
            label="Drain",
        )
        ax.plot(
            -10000,
            -10000,
            lw=0,
            marker="s",
            ms=10,
            mfc="0.5",
            mec="0.5",
            label="Well",
        )
        ax.plot(
            -10000,
            -10000,
            lw=0,
            marker=u"$\u2192$",
            ms=10,
            mfc="0.75",
            mec="0.75",
            label="Normalized specific discharge",
        )
        ax.plot(
            -10000, -10000, lw=0.5, color="black", label=r"Head contour, $ft$"
        )
        fs.graph_legend(ax, loc="upper center")

        # plot colorbar
        cax = plt.axes([0.325, 0.125, 0.35, 0.025])
        cbar = plt.colorbar(
            plot_obj, shrink=0.8, orientation="horizontal", cax=cax
        )
        cbar.ax.tick_params(size=0)
        cbar.ax.set_xlabel(r"Head, $ft$", fontsize=9)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..", "figures", "{}{}".format(sim_name, config.figure_ext)
            )
            fig.savefig(fpth)


# Function that wraps all of the steps for the TWRI model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(silent=True):
    sim = build_model()
    mf = build_mf5model()

    write_model(sim, mf, silent=silent)

    success = run_model(sim, mf, silent=silent)

    if success:
        plot_results(sim, mf)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    simulation(silent=False)


# nosetest end

if __name__ == "__main__":
    # ### TWRI Simulation
    #
    # This simulation uses the recharge rate defined in the original problem.
    # Simulated heads in model layers 1, 2, and 3 are shown in the figure below.
    # The location of drain (green) and (gray) well boundary conditions, normalized
    # specific discharge, and head contours (25 ft contour intervals) are also shown.
    #

    simulation()
