# # Comparison of MF6 transport with MT3DMS Example Problems
# The purpose of this script is to (1) recreate the example problems that were first described in the 1999 MT3DMS report, and (2) compare MF6-GWT solutions to the established MT3DMS solutions.
#
# Ten example problems appear in the 1999 MT3DMS manual, starting on page 130.  This notebook demonstrates example 10 from the list below:
#   1. One-Dimensional Transport in a Uniform Flow Field
#   2. One-Dimensional Transport with Nonlinear or Nonequilibrium Sorption
#   3. **Two-Dimensional Transport in a Uniform Flow Field**
#   4. Two-Dimensional Transport in a Diagonal Flow Field
#   5. Two-Dimensional Transport in a Radial Flow Field
#   6. Concentration at an Injection/Extraction Well
#   7. Three-Dimensional Transport in a Uniform Flow Field
#   8. Two-Dimensional, Vertical Transport in a Heterogeneous Aquifer
#   9. Two-Dimensional Application Example
#   10. Three-Dimensional Field Case Study

# Append to system path to include the common subdirectory
import os
import sys
sys.path.append(os.path.join("..", "common"))

# #### Imports
import matplotlib.pyplot as plt
import flopy
import numpy as np
import config
from figspecs import USGSFigure

mf6exe = os.path.abspath(config.mf6_exe)
assert os.path.isfile(mf6exe)
print(mf6exe)
exe_name_mf = config.mf2005_exe
exe_name_mt = config.mt3dms_exe

# Set figure properties specific to this problem
figure_size = (6, 4.5)

# Base simulation and model name and workspace
ws = config.base_ws
example_name = "ex-gwt-mt3dms-p03"

# #### Model units
length_units = "meters"
time_units = "days"

# Table
nlay = 1      # Number of layers
nrow = 31     # Number of rows
ncol = 46     # Number of columns
delr = 10.0   # Column width ($m$)
delc = 10.0   # Row width ($m$)
delz = 10.0   # Layer thickness ($m$)
top =  0.0    # Top of the model ($m$)
prsity = 0.3  # Porosity
perlen = 365  # Simulation time ($days$)
k11 = 1.0     # Horizontal hydraulic conductivity ($m/d$)
qwell = 1.    # Volumetric injection rate ($m^3/d$)
cwell = 1000. # Concentration of injected water ($mg/L$)
al = 10.      # Longitudinal dispersivity ($m$)
trpt = 0.3    # Ratio of transverse to longitudinal dispersivity

# #### Additional model input
perlen = [1, 365.]
nper = len(perlen)
nstp = [2, 730]
tsmult = [1., 1.]
sconc = 0.
dt0 = 0.3
ath1 = al * trpt
dmcoef = 0.
xt3d = [False]

botm = [top - delz] # Model geometry

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
icelltype = 0

# #### Initial conditions
Lx = (ncol - 1) * delr
v = 1. / 3.
prsity = 0.3
q = v * prsity
h1 = q * Lx
strt = np.zeros((nlay, nrow, ncol), dtype=np.float)
strt[0, :, 0] = h1

ibound_mf2k5 = np.ones((nlay, nrow, ncol), dtype=np.int)
ibound_mf2k5[0, :,  0] = -1
ibound_mf2k5[0, :, -1] = -1
idomain = np.ones((nlay, nrow, ncol), dtype=np.int)
icbund = 1
c0 = 0.
cncspd = [[(0, 0, 0), c0]]  
welspd = {0: [[0, 15, 15, qwell]]}  # Well pumping info for MF2K5
spd    = {0: [0, 15, 15, cwell, 2]} # Well pupming info for MT3DMS
#              (k,  i,  j),  flow, conc
spd_mf6 = {0:[[(0, 15, 15), qwell,  cwell]]} # MF6 pumping information


# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.
ttsmult = 1.0
percel = 1.0    # HMOC parameters in case they are invoked
itrack = 3      # HMOC
wd = 0.5        # HMOC
dceps = 1.e-5   # HMOC
nplane = 1      # HMOC
npl = 0         # HMOC
nph = 16        # HMOC
npmin = 4       # HMOC
npmax = 32      # HMOC   
dchmoc=1.e-3    # HMOC   
nlsink = nplane # HMOC
npsink = nph    # HMOC

# ## Static temporal data used by TDIS file
tdis_rc = []
tdis_rc.append((perlen, nstp, 1.0))            


# ## Function to build models
# MODFLOW 6 flopy simulation object (sim) is returned if building the model
def build_model(sim_name, mixelm=0, silent=False):
    if config.buildModel:
        
        mt3d_ws = os.path.join(ws, sim_name, 'mt3d')
        modelname_mf = 'p03_mf'
        
        # Instantiate the MODFLOW model
        mf = flopy.modflow.Modflow(modelname=modelname_mf,
                                   model_ws=mt3d_ws,
                                   exe_name=exe_name_mf
        )
        
        # Instantiate discretization package
        # units: itmuni=4 (days), lenuni=2 (m)
        flopy.modflow.ModflowDis(mf,
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
                                 lenuni=2
        )
        
        # Instantiate basic package
        flopy.modflow.ModflowBas(mf,
                                 ibound=ibound_mf2k5,
                                 strt=strt
        )
        
        # Instantiate layer property flow package
        flopy.modflow.ModflowLpf(mf,
                                 hk=k11,
                                 laytyp=icelltype
        )
        
        # Instantiate well package
        flopy.modflow.ModflowWel(mf, stress_period_data=welspd)
        
        # Instantiate solver package
        flopy.modflow.ModflowPcg(mf)
        
        # Instantiate link mass transport package (for writing linker file)
        flopy.modflow.ModflowLmt(mf)
        
        # Transport
        modelname_mt = 'p03_mt'
        mt = flopy.mt3d.Mt3dms(modelname=modelname_mt,
                               model_ws=mt3d_ws,
                               exe_name=exe_name_mt,
                               modflowmodel=mf
        )
        
        # Instantiate basic transport package
        flopy.mt3d.Mt3dBtn(mt,
                           icbund=icbund,
                           prsity=prsity,
                           sconc=sconc,
                           nstp=nstp,
                           perlen=perlen
        )
        
        # Instatiate the advection package
        flopy.mt3d.Mt3dAdv(mt,
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
                           wd=wd
        )
        
        # Instantiate the dispersion package
        flopy.mt3d.Mt3dDsp(mt,
                           al=al,
                           trpt=trpt,
                           dmcoef=dmcoef
        )
        
        # Instantiate the source/sink mixing package
        flopy.mt3d.Mt3dSsm(mt, stress_period_data=spd)
        
        # Instantiate the GCG solver in MT3DMS
        flopy.mt3d.Mt3dGcg(mt,
                           mxiter=10
        )
        
        # MODFLOW 6
        name = 'p01_mf6'
        gwfname = 'gwf_' + name
        sim_ws = os.path.join(ws, sim_name)
        sim = flopy.mf6.MFSimulation(sim_name=sim_name, 
                                     sim_ws=sim_ws, 
                                     exe_name=mf6exe
        )
        
        # Instantiating MODFLOW 6 time discretization
        tdis_rc = []
        for i in range(nper):
            tdis_rc.append((perlen[i], nstp[i], tsmult[i]))
        flopy.mf6.ModflowTdis(sim, 
                              nper=nper, 
                              perioddata=tdis_rc, 
                              time_units=time_units
        )
        
        # Instantiating MODFLOW 6 groundwater flow model
        gwf = flopy.mf6.ModflowGwf(sim,
                                   modelname=gwfname,
                                   save_flows=True,
                                   model_nam_file='{}.nam'.format(gwfname)
        )
        
        # Instantiating MODFLOW 6 solver for flow model
        imsgwf = flopy.mf6.ModflowIms(sim,
                                      print_option='SUMMARY',
                                      outer_dvclose=hclose,
                                      outer_maximum=nouter,
                                      under_relaxation='NONE',
                                      inner_maximum=ninner,
                                      inner_dvclose=hclose, rcloserecord=rclose,
                                      linear_acceleration='CG',
                                      scaling_method='NONE',
                                      reordering_method='NONE',
                                      relaxation_factor=relax,
                                      filename='{}.ims'.format(gwfname)
        )
        sim.register_ims_package(imsgwf, [gwf.name])
        
        # Instantiating MODFLOW 6 discretization package
        flopy.mf6.ModflowGwfdis(gwf,
                                length_units=length_units,
                                nlay=nlay,
                                nrow=nrow,
                                ncol=ncol,
                                delr=delr,
                                delc=delc,
                                top=top,
                                botm=botm,
                                idomain=np.ones((nlay, nrow, ncol),
                                                dtype=np.int),
                                filename='{}.dis'.format(gwfname)
        )
        
        # Instantiating MODFLOW 6 node-property flow package
        flopy.mf6.ModflowGwfnpf(gwf,
                                save_flows=False,
                                icelltype=icelltype, 
                                k=k11, 
                                k33=k33,
                                save_specific_discharge = True,
                                filename='{}.npf'.format(gwfname)
        )
        
        # Instantiating MODFLOW 6 initial conditions package for flow model
        flopy.mf6.ModflowGwfic(gwf, 
                               strt=strt,
                               filename='{}.ic'.format(gwfname)
        )
        
        # Instantiate MODFLOW 6 storage package 
        sto = flopy.mf6.ModflowGwfsto(gwf, 
                                      ss=0, 
                                      sy=0,
                                      filename='{}.sto'.format(gwfname))
        
        # Instantiating MODFLOW 6 constant head package
        rowList = np.arange(0,nrow).tolist()
        chdspd = []
        # Loop through the left & right sides.
        for itm in rowList:
            # first, do left side of model
            chdspd.append([(0,itm,0), h1])
            # finally, do right side of model
            chdspd.append([(0, itm, ncol - 1), 0.])
        
        chdspd = {0: chdspd} 
        flopy.mf6.ModflowGwfchd(gwf,
                                maxbound=len(chdspd), 
                                stress_period_data=chdspd,
                                save_flows=False,
                                pname='CHD-1',
                                filename='{}.chd'.format(gwfname)
        )
        
        # Instantiate the wel package
        flopy.mf6.ModflowGwfwel(gwf,
                                print_input=True,
                                print_flows=True,
                                stress_period_data=spd_mf6,
                                save_flows=False,
                                auxiliary='CONCENTRATION',
                                pname='WEL-1',
                                filename='{}.wel'.format(gwfname))
        
        # Instantiating MODFLOW 6 output control package for flow model
        flopy.mf6.ModflowGwfoc(gwf,
                               head_filerecord="{}.hds".format(gwfname),
                               budget_filerecord="{}.bud".format(gwfname),
                               headprintrecord=[
                                                ('COLUMNS', 10, 'WIDTH', 15,
                                                 'DIGITS', 6, 'GENERAL')],
                               saverecord=[('HEAD', 'LAST'),
                                           ('BUDGET', 'LAST')],
                               printrecord=[('HEAD', 'LAST'),
                                            ('BUDGET', 'LAST')]
        )
        
        # Instantiating MODFLOW 6 groundwater transport package
        gwtname = 'gwt_' + name
        gwt = flopy.mf6.MFModel(sim,
                                model_type='gwt6',
                                modelname=gwtname,
                                model_nam_file='{}.nam'.format(gwtname))
        gwt.name_file.save_flows = True
        
        # create iterative model solution and register the gwt model with it  
        imsgwt = flopy.mf6.ModflowIms(sim,
                             print_option='SUMMARY',
                             outer_dvclose=hclose,
                             outer_maximum=nouter,
                             under_relaxation='NONE',
                             inner_maximum=ninner,
                             inner_dvclose=hclose, rcloserecord=rclose,
                             linear_acceleration='BICGSTAB',
                             scaling_method='NONE',
                             reordering_method='NONE',
                             relaxation_factor=relax,
                             filename='{}.ims'.format(gwtname))
        sim.register_ims_package(imsgwt, [gwt.name])
        
        # Instantiating MODFLOW 6 transport discretization package
        flopy.mf6.ModflowGwtdis(gwt, 
                                nlay=nlay,
                                nrow=nrow,
                                ncol=ncol,
                                delr=delr, 
                                delc=delc,
                                top=top,
                                botm=botm,
                                idomain=1,
                                filename='{}.dis'.format(gwtname)
        )
        
        # Instantiating MODFLOW 6 transport initial concentrations
        flopy.mf6.ModflowGwtic(gwt, 
                               strt=sconc,
                               filename='{}.ic'.format(gwtname)
        )
        
        # Instantiating MODFLOW 6 transport advection package
        if mixelm == 0:
            scheme = 'UPSTREAM'
        elif mixelm == -1:
            scheme = 'TVD'
        else:
            raise Exception()
        flopy.mf6.ModflowGwtadv(gwt, 
                                scheme=scheme,
                                filename='{}.adv'.format(gwtname)
        )
        
        # Instantiating MODFLOW 6 transport dispersion package
        if al != 0:
            flopy.mf6.ModflowGwtdsp(gwt, 
                                    xt3d=xt3d,
                                    alh=al,
                                    ath1=ath1,
                                    filename='{}.dsp'.format(gwtname)
        )

        # Instantiating MODFLOW 6 transport mass storage package (formerly "reaction" package in MT3DMS)
        flopy.mf6.ModflowGwtmst(gwt, 
                                porosity=prsity,
                                first_order_decay=False,
                                decay=None, 
                                decay_sorbed=None,
                                sorbtion=False, 
                                bulk_density=None, 
                                distcoef=None,
                                filename='{}.mst'.format(gwtname)
        )

        # Instantiating MODFLOW 6 transport constant concentration package
        flopy.mf6.ModflowGwtcnc(gwt, 
                                maxbound=len(cncspd),
                                stress_period_data=cncspd,
                                save_flows=False,
                                pname='CNC-1',
                                filename='{}.cnc'.format(gwtname)
        )

        # Instantiating MODFLOW 6 transport source-sink mixing package
        sourcerecarray = [('WEL-1', 'AUX', 'CONCENTRATION')]
        flopy.mf6.ModflowGwtssm(gwt, 
                                sources=sourcerecarray,
                                filename='{}.ssm'.format(gwtname)
        )

        # Instantiating MODFLOW 6 transport output control package
        flopy.mf6.ModflowGwtoc(gwt,
                               budget_filerecord='{}.cbc'.format(gwtname),
                               concentration_filerecord='{}.ucn'.format(
                                   gwtname),
                               concentrationprintrecord=[
                                   ('COLUMNS', 10, 'WIDTH', 15,
                                    'DIGITS', 6, 'GENERAL')],
                               saverecord=[('CONCENTRATION', 'LAST'),
                                           ('BUDGET', 'LAST')],
                               printrecord=[('CONCENTRATION', 'LAST'),
                                            ('BUDGET', 'LAST')]
        )

        # Instantiating MODFLOW 6 flow-transport exchange mechanism
        flopy.mf6.ModflowGwfgwt(sim, 
                                exgtype='GWF6-GWT6',
                                exgmnamea=gwfname, 
                                exgmnameb=gwtname,
                                filename='{}.gwfgwt'.format(name)
        )
        return mf, mt, sim
    return None

# ## Function to write model files
def write_model(mf2k5, mt3d, sim, silent=False):
    if config.writeModel:
        mf2k5.write_input()
        mt3d.write_input()
        sim.write_simulation(silent=True)

# ## Function to run the model
# _True_ is returned if the model runs successfully
def run_model(mf2k5, mt3d, sim, silent=False):
    success = True
    if config.runModel:
        success, buff = mf2k5.run_model(silent=silent)
        success, buff = mt3d.run_model(silent=silent)
        success, buff = sim.run_simulation(silent=False)
        if not success:
            print(buff)
    return success

# ## Function to plot the model results
def plot_results(mt3d, mf6, idx, ax=None):
    if config.plotModel:
        mt3d_out_path = mt3d.model_ws
        mf6_out_path = mf6.simulation_data.mfpath.get_sim_path()
        mf6.simulation_data.mfpath.get_sim_path()

        # Get the MT3DMS concentration output
        fname_mt3d = os.path.join(mt3d_out_path, 'MT3D001.UCN')
        ucnobj_mt3d = flopy.utils.UcnFile(fname_mt3d)
        times_mt3d = ucnobj_mt3d.get_times()
        conc_mt3d = ucnobj_mt3d.get_alldata()

        # Get the MF6 concentration output
        fname_mf6 = os.path.join(mf6_out_path, list(mf6.model_names)[1] + '.ucn')
        ucnobj_mf6 = flopy.utils.HeadFile(fname_mf6, precision='double',
                                          text='CONCENTRATION')
        
        times_mf6 = ucnobj_mf6.get_times() 
        conc_mf6 = ucnobj_mf6.get_alldata()
        
        # Create figure for scenario 
        fs = USGSFigure(figure_type="graph", verbose=False)
        sim_name = mf6.name
        plt.rcParams['lines.dashed_pattern'] = [5.0, 5.0]
        if ax is None:
            fig = plt.figure(figsize=figure_size, dpi=300, tight_layout=True)
            ax = fig.add_subplot(1, 1, 1,  aspect='equal')
        
        mm = flopy.plot.PlotMapView(model=mt3d)
        mm.plot_grid(color='.5', alpha=0.2)
        cs1 = mm.contour_array(conc_mt3d[1], levels=[0.1, 1., 10., 50.], colors='k')
        plt.clabel(cs1, inline=1, fontsize=10)
        cs2 = mm.contour_array(conc_mf6[1], levels=[0.1, 1., 10., 50.], colors='r', linestyles='--')
        plt.clabel(cs2, inline=1, fontsize=10)
        labels = ['MT3DMS', 'MODFLOW6']
        lines = [cs1.collections[0], cs2.collections[0]]
        
        plt.xlabel('DISTANCE ALONG X-AXIS, IN METERS')
        plt.ylabel('DISTANCE ALONG Y-AXIS, IN METERS')
        title = 'Plume at Time = 365 ' + '{}'.format(time_units)
        
        ax.legend(lines, labels, loc='upper left')
    
        #ax.plot(np.linspace(0, l, ncol), conc_mt3d[0,0,0,:], color='k', label='MT3DMS')
        #ax.plot(np.linspace(0, l, ncol), conc_mf6[0,0,0,:], '^', color='b', label='MF6')
        #ax.set_ylim(0, 1.2)
        #ax.set_xlim(0, 1000)
        #ax.set_xlabel('Distance, in m')
        #ax.set_ylabel('Concentration')
        #title = 'Concentration Profile at Time = 1,000 ' + '{}'.format(
        #                                                        time_units)
        #ax.legend()
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
    mf2k5, mt3d, sim = build_model(example_name)
    write_model(mf2k5, mt3d, sim, silent=False)
    success = run_model(mf2k5, mt3d, sim, silent=False)

    if success:
        plot_results(mt3d, sim, idx)


# nosetest - exclude block from this nosetest to the next nosetest
def test_01():
    scenario(0)

# nosetest end

if __name__ == "__main__":
    # ## Scenario 1
    # Two-dimensional transport in a uniform flow field
    scenario(0)
