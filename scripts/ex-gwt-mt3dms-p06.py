# ## Concentration at an Injection/Extraction Well, Comparison of MODFLOW 6 transport with MT3DMS
#
# The purpose of this script is to (1) recreate the example problems that were first
# described in the 1999 MT3DMS report, and (2) compare MF6-GWT solutions to the
# established MT3DMS solutions.
#
# Ten example problems appear in the 1999 MT3DMS manual, starting on page 130.  This
# notebook demonstrates example 10 from the list below:
#
#   1. One-Dimensional Transport in a Uniform Flow Field
#   2. One-Dimensional Transport with Nonlinear or Nonequilibrium Sorption
#   3. Two-Dimensional Transport in a Uniform Flow Field
#   4. Two-Dimensional Transport in a Diagonal Flow Field
#   5. Two-Dimensional Transport in a Radial Flow Field
#   6. _Concentration at an Injection/Extraction Well_
#   7. Three-Dimensional Transport in a Uniform Flow Field
#   8. Two-Dimensional, Vertical Transport in a Heterogeneous Aquifer
#   9. Two-Dimensional Application Example
#   10. Three-Dimensional Field Case Study


# ### MODFLOW 6 GWT MT3DMS Example 6 Problem Setup

# Append to system path to include the common subdirectory

import os
import sys

sys.path.append(os.path.join("..", "common"))

# Imports

import config
import flopy
import matplotlib.pyplot as plt
import numpy as np
from modflow_devtools.figspec import USGSFigure

mf6exe = "mf6"
exe_name_mf = "mf2005"
exe_name_mt = "mt3dusgs"

# Set figure properties specific to this problem

figure_size = (6, 4.5)

# Base simulation and model name and workspace

ws = config.base_ws
example_name = "ex-gwt-mt3dms-p06"

# Model units

length_units = "feet"
time_units = "days"

# Table

nlay = 1  # Number of layers
nrow = 31  # Number of rows
ncol = 31  # Number of columns
delr = 900.0  # Column width ($ft$)
delc = 900.0  # Row width ($ft$)
delz = 20.0  # Layer thickness ($ft$)
top = 0.0  # Top of the model ($ft$)
prsity = 0.35  # Porosity
dum1 = 2.5  # Length of the injection period ($years$)
dum2 = 7.5  # Length of the extraction period ($years$)
k11 = 432.0  # Horizontal hydraulic conductivity ($ft/d$)
qwell = 1.0  # Volumetric injection rate ($ft^3/d$)
cwell = 100.0  # Relative concentration of injected water ($\%$)
al = 100.0  # Longitudinal dispersivity ($ft$)
trpt = 1.0  # Ratio of transverse to longitudinal dispersitivity

# Additional model input

perlen = [912.5, 2737.5]
nper = len(perlen)
nstp = [365, 1095]
tsmult = [1.0, 1.0]
k11 = (
    0.005 * 86400
)  # established above, but explicitly writing out its origin here
sconc = 0.0
c0 = 0.0
dt0 = 56.25
dmcoef = 0
ath1 = al * trpt
botm = [top - delz]  # Model geometry
k33 = k11  # Vertical hydraulic conductivity ($m/d$)
icelltype = 0
mixelm = -1
strt = np.zeros((nlay, nrow, ncol), dtype=float)

# Active model domain

ibound_mf2k5 = np.ones((nlay, nrow, ncol), dtype=int) * -1
ibound_mf2k5[:, 1 : nrow - 1, 1 : ncol - 1] = 1
idomain = np.ones((nlay, nrow, ncol), dtype=int)
icbund = 1

# Boundary conditions
# MF2K5 pumping info:

qwell = 86400.0
welspd = {
    0: [[0, 15, 15, qwell]],  # Well pumping info for MF2K5
    1: [[0, 15, 15, -qwell]],
}
cwell = 100.0
spd = {
    0: [0, 15, 15, cwell, 2],  # Well pupming info for MT3DMS
    1: [0, 15, 15, 0.0, 2],
}

# MF6 pumping information
#          (k,  i,  j),  flow,  conc
spd_mf6 = {0: [[(0, 15, 15), qwell, cwell]], 1: [[(0, 15, 15), -qwell, 0.0]]}

# MF6 constant head boundaries:

chdspd = []
# Loop through the left & right sides.
for i in np.arange(nrow):
    chdspd.append([(0, i, 0), strt[0, i, 0]])
    chdspd.append([(0, i, ncol - 1), strt[0, i, ncol - 1]])
# Loop through the top & bottom while omitting the corner cells
for j in np.arange(1, ncol - 1):
    chdspd.append([(0, 0, j), strt[0, 0, j]])
    chdspd.append([(0, nrow - 1, j), strt[0, nrow - 1, j]])

chdspd = {0: chdspd}

# Solver settings

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0
percel = 1.0  # HMOC parameters
itrack = 3
wd = 0.5
dceps = 1.0e-5
nplane = 1
npl = 0
nph = 16
npmin = 2
npmax = 32
dchmoc = 1.0e-3
nlsink = nplane
npsink = nph

# Static temporal data used by TDIS file

tdis_rc = []
tdis_rc.append((perlen, nstp, 1.0))

# ### Functions to build, write, and run models and plot MT3DMS Example 10 Problem results
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(sim_name, mixelm=0, silent=False):
    if config.buildModel:
        mt3d_ws = os.path.join(ws, sim_name, "mt3d")
        modelname_mf = "p06-mf"

        # Instantiate the MODFLOW model
        mf = flopy.modflow.Modflow(
            modelname=modelname_mf, model_ws=mt3d_ws, exe_name=exe_name_mf
        )

        # Instantiate discretization package
        # units: itmuni=4 (days), lenuni=2 (m)
        flopy.modflow.ModflowDis(
            mf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
            nper=nper,
            nstp=nstp,
            perlen=perlen,
            itmuni=4,
            lenuni=1,
        )

        # Instantiate basic package
        flopy.modflow.ModflowBas(mf, ibound=ibound_mf2k5, strt=strt)

        # Instantiate layer property flow package
        flopy.modflow.ModflowLpf(mf, hk=k11, laytyp=icelltype)

        # Instantiate well package
        flopy.modflow.ModflowWel(mf, stress_period_data=welspd)

        # Instantiate solver package
        flopy.modflow.ModflowSip(mf)

        # Instantiate link mass transport package (for writing linker file)
        flopy.modflow.ModflowLmt(mf)

        # Transport
        modelname_mt = "p06-mt"
        mt = flopy.mt3d.Mt3dms(
            modelname=modelname_mt,
            model_ws=mt3d_ws,
            exe_name=exe_name_mt,
            modflowmodel=mf,
        )

        # Instantiate basic transport package
        flopy.mt3d.Mt3dBtn(
            mt,
            icbund=icbund,
            prsity=prsity,
            sconc=sconc,
            nper=nper,
            perlen=perlen,
            dt0=dt0,
            obs=[(0, 15, 15)],
        )

        # Instatiate the advection package
        flopy.mt3d.Mt3dAdv(
            mt,
            mixelm=mixelm,
            dceps=dceps,
            nplane=nplane,
            npl=npl,
            nph=nph,
            npmin=npmin,
            npmax=npmax,
            nlsink=nlsink,
            npsink=npsink,
            percel=percel,
            itrack=itrack,
            wd=wd,
        )

        # Instantiate the dispersion package
        flopy.mt3d.Mt3dDsp(mt, al=al, trpt=trpt, dmcoef=dmcoef)

        # Instantiate the source/sink mixing package
        flopy.mt3d.Mt3dSsm(mt, stress_period_data=spd)

        # Instantiate the GCG solver in MT3DMS
        flopy.mt3d.Mt3dGcg(mt)

        # MODFLOW 6
        name = "p06-mf6"
        gwfname = "gwf-" + name
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation(
            sim_name=sim_name, sim_ws=sim_ws, exe_name=mf6exe
        )

        # Instantiating MODFLOW 6 time discretization
        tdis_rc = []
        for i in range(nper):
            tdis_rc.append((perlen[i], nstp[i], tsmult[i]))
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_rc, time_units=time_units
        )

        # Instantiating MODFLOW 6 groundwater flow model
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=gwfname,
            save_flows=True,
            model_nam_file=f"{gwfname}.nam",
        )

        # Instantiating MODFLOW 6 solver for flow model
        imsgwf = flopy.mf6.ModflowIms(
            sim,
            print_option="SUMMARY",
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
            filename=f"{gwfname}.ims",
        )
        sim.register_ims_package(imsgwf, [gwf.name])

        # Instantiating MODFLOW 6 discretization package
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
            idomain=idomain,
            filename=f"{gwfname}.dis",
        )

        # Instantiating MODFLOW 6 node-property flow package
        flopy.mf6.ModflowGwfnpf(
            gwf,
            save_flows=False,
            icelltype=icelltype,
            k=k11,
            k33=k33,
            save_specific_discharge=True,
            filename=f"{gwfname}.npf",
        )

        # Instantiating MODFLOW 6 storage package (steady flow conditions, so no actual storage, using to print values in .lst file)
        flopy.mf6.ModflowGwfsto(gwf, ss=0, sy=0, filename=f"{gwfname}.sto")

        # Instantiating MODFLOW 6 initial conditions package for flow model
        flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

        # Instantiating MODFLOW 6 constant head package
        flopy.mf6.ModflowGwfchd(
            gwf,
            maxbound=len(chdspd),
            stress_period_data=chdspd,
            save_flows=False,
            pname="CHD-1",
            filename=f"{gwfname}.chd",
        )

        # Instantiate the wel package
        flopy.mf6.ModflowGwfwel(
            gwf,
            print_input=True,
            print_flows=True,
            stress_period_data=spd_mf6,
            save_flows=False,
            auxiliary="CONCENTRATION",
            pname="WEL-1",
            filename=f"{gwfname}.wel",
        )

        # Instantiating MODFLOW 6 output control package for flow model
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=f"{gwfname}.hds",
            budget_filerecord=f"{gwfname}.bud",
            headprintrecord=[
                ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
            ],
            saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
            printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        )

        # Instantiating MODFLOW 6 groundwater transport package
        gwtname = "gwt_" + name
        gwt = flopy.mf6.MFModel(
            sim,
            model_type="gwt6",
            modelname=gwtname,
            model_nam_file=f"{gwtname}.nam",
        )
        gwt.name_file.save_flows = True

        # create iterative model solution and register the gwt model with it
        imsgwt = flopy.mf6.ModflowIms(
            sim,
            print_option="SUMMARY",
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
            filename=f"{gwtname}.ims",
        )
        sim.register_ims_package(imsgwt, [gwt.name])

        # Instantiating MODFLOW 6 transport discretization package
        flopy.mf6.ModflowGwtdis(
            gwt,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
            idomain=idomain,
            filename=f"{gwtname}.dis",
        )

        # Instantiating MODFLOW 6 transport initial concentrations
        flopy.mf6.ModflowGwtic(gwt, strt=sconc, filename=f"{gwtname}.ic")

        # Instantiating MODFLOW 6 transport advection package
        if mixelm >= 0:
            scheme = "UPSTREAM"
        elif mixelm == -1:
            scheme = "TVD"
        else:
            raise Exception()
        flopy.mf6.ModflowGwtadv(gwt, scheme=scheme, filename=f"{gwtname}.adv")

        # Instantiating MODFLOW 6 transport dispersion package
        if al != 0:
            flopy.mf6.ModflowGwtdsp(
                gwt,
                xt3d_off=True,
                alh=al,
                ath1=ath1,
                filename=f"{gwtname}.dsp",
            )

        # Instantiating MODFLOW 6 transport mass storage package (formerly "reaction" package in MT3DMS)
        flopy.mf6.ModflowGwtmst(
            gwt,
            porosity=prsity,
            first_order_decay=False,
            decay=None,
            decay_sorbed=None,
            sorption=None,
            bulk_density=None,
            distcoef=None,
            filename=f"{gwtname}.mst",
        )

        # Instantiating MODFLOW 6 transport source-sink mixing package
        sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
        flopy.mf6.ModflowGwtssm(
            gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
        )

        # Instantiating MODFLOW 6 transport output control package
        flopy.mf6.ModflowGwtoc(
            gwt,
            budget_filerecord=f"{gwtname}.cbc",
            concentration_filerecord=f"{gwtname}.ucn",
            concentrationprintrecord=[
                ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
            ],
            saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
            printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        )

        # Instantiate observation package (for transport)
        obslist = [["bckgrnd_cn", "concentration", (0, 15, 15)]]
        obsdict = {f"{gwtname}.obs.csv": obslist}
        obs = flopy.mf6.ModflowUtlobs(
            gwt, print_input=False, continuous=obsdict
        )

        # Instantiating MODFLOW 6 flow-transport exchange mechanism
        flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=gwfname,
            exgmnameb=gwtname,
            filename=f"{name}.gwfgwt",
        )
        return mf, mt, sim
    return None


# Function to write model files


def write_model(mf2k5, mt3d, sim, silent=True):
    if config.writeModel:
        mf2k5.write_input()
        mt3d.write_input()
        sim.write_simulation(silent=silent)


# Function to run the model. True is returned if the model runs successfully.


@config.timeit
def run_model(mf2k5, mt3d, sim, silent=True):
    success = True
    if config.runModel:
        success, buff = mf2k5.run_model(silent=silent)
        success, buff = mt3d.run_model(silent=silent)
        success, buff = sim.run_simulation(silent=silent)
        if not success:
            print(buff)
    return success


# Function to plot the model results


def plot_results(mt3d, mf6, idx, ax=None):
    if config.plotModel:
        mt3d_out_path = mt3d.model_ws
        mf6_out_path = mf6.simulation_data.mfpath.get_sim_path()
        mf6.simulation_data.mfpath.get_sim_path()

        # Get the MT3DMS observation output file
        fname = os.path.join(mt3d_out_path, "MT3D001.OBS")
        if os.path.isfile(fname):
            cvt = mt3d.load_obs(fname)
        else:
            cvt = None

        # Get the MODFLOW 6 concentration observation output file
        fname = os.path.join(
            mf6_out_path, list(mf6.model_names)[1] + ".obs.csv"
        )
        mf6cobs = flopy.utils.Mf6Obs(fname).data

        # Create figure for scenario
        fs = USGSFigure(figure_type="graph", verbose=False)
        sim_name = mf6.name
        plt.rcParams["lines.dashed_pattern"] = [5.0, 5.0]
        if ax is None:
            fig = plt.figure(figsize=figure_size, dpi=300, tight_layout=True)
            ax = fig.add_subplot(1, 1, 1)

        x = cvt["time"] / 365.0
        y = cvt["(1, 16, 16)"]
        # Pare down the list length to clean plot
        x_pare = x[::20]
        y_pare = y[::20]
        ax.plot(x_pare, y_pare, label="Upstream FD", marker="^")

        # Add MF6 output
        x_mf6 = mf6cobs["totim"] / 365.0
        y_mf6 = mf6cobs["BCKGRND_CN"]
        x_mf6_pare = x_mf6[::20]
        y_mf6_pare = y_mf6[::20]
        ax.plot(
            x_mf6_pare,
            y_mf6_pare,
            label="MODFLOW 6",
            marker="x",
            linestyle=":",
        )

        plt.xlim(0, 10)
        plt.ylim(0, 100.0)
        plt.xlabel("Time, in years")
        plt.ylabel("Normalized Concentration, in percent")
        plt.legend()
        title = "Calculated Concentration at an Injection/Pumping Well"

        letter = chr(ord("@") + idx + 1)
        fs.heading(letter=letter, heading=title)

        # save figure
        if config.plotSave:
            fpth = os.path.join(
                "..",
                "figures",
                f"{sim_name}{config.figure_ext}",
            )
            fig.savefig(fpth)


# ### Function that wraps all of the steps for each MT3DMS Example 10 Problem scenario
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.


def scenario(idx, silent=True):
    mf2k5, mt3d, sim = build_model(example_name, mixelm=mixelm)

    write_model(mf2k5, mt3d, sim, silent=silent)

    success = run_model(mf2k5, mt3d, sim, silent=silent)

    if success:
        plot_results(mt3d, sim, idx)


# nosetest - exclude block from this nosetest to the next nosetest


def test_01():
    scenario(0, silent=False)


# nosetest end

if __name__ == "__main__":
    # ### Two-Dimensional Transport in a Diagonal Flow Field
    #
    # Compares the standard finite difference solutions between MT3D MF 6
    scenario(0)
