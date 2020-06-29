# A script for generating problem 1 from the MT3DMS documentation (Zheng & 
# Wang, 1999; pg 130).  

import os
import sys
import numpy as np
import argparse
import distutils   # used in argparse script
import matplotlib.pyplot as plt
from flopy.utils.util_array import read1d

# Import the analytical solutions for verification
import analytical

try:
    import pymake
except:
    msg = 'Error. Pymake package is not available.\n'
    msg += 'Try installing using the following command:\n'
    msg += ' pip install https://github.com/modflowpy/pymake/zipball/master'
    raise Exception(msg)

try:
    import flopy
except:
    msg = 'Error. FloPy package is not available.\n'
    msg += 'Try installing using the following command:\n'
    msg += ' pip install flopy'
    raise Exception(msg)

mf6path = '../../modflow6.git'
assert os.path.isdir(mf6path)

exe_ext = '.exe'

mf6exe = os.path.join(mf6path, 'bin', 'mf6' + exe_ext)
mf6exe = os.path.abspath(mf6exe)
assert os.path.isfile(mf6exe)
print(mf6exe)

exe_name_mf = os.path.join(mf6path,'autotest','temp','mfexes','mf2005'+exe_ext)
exe_name_mt = os.path.join(mf6path,'autotest','temp','mfexes','mt3dms'+exe_ext)
exe_name_mf = os.path.abspath(exe_name_mf)
exe_name_mt = os.path.abspath(exe_name_mt)

assert os.path.isfile(exe_name_mt)
assert os.path.isfile(exe_name_mf)
print(exe_name_mt)
print(exe_name_mf)

testdir = '../examples/ex_mt3dms_p01'

# -------------------------------------------------------
# Layout variables used by both MT3DMS and MF6-transport
# -------------------------------------------------------

# Discretization related variables
# --------------------------------
nlay = 1
nrow = 1
ncol = 101
delr = 10.
delc = 1
delv = 1.
# Model geometry
top  = 0.
botm = [top - delv]
# Time discretization
perlen = 2000.
nstp = 100.
dt0 = perlen / nstp

# Aquifer properties
# ------------------
hk = 1.
laytyp = 1

# Initial conditions
# ------------------
# Starting Heads:
prsity = 0.25
Lx = (ncol - 1) * delr
v = 0.24
q = v * prsity
h1 = q * Lx
strt = np.zeros((nlay, nrow, ncol), dtype=np.float)
strt[0, 0, 0] = h1
l = 1000.            # Needed for plots

# Active model domain
ibound = np.ones((nlay, nrow, ncol), dtype=np.int)
ibound[0, 0, 0] = -1
ibound[0, 0, -1] = -1

# Transport related variables
# ---------------------------
mixelm = 0                              # TVD
# Additional reactive transport related terms defined when calling the 
# model creation function
rhob = 0.25
sp2 = 0.                                # read, but not used in this problem
    
# Starting concentrations:
sconc = np.zeros((nlay, nrow, ncol), dtype=np.float)

# Molecular diffusion coefficient
dmcoef = 0.

# Solver settings (and related)
# -----------------------------
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.
ttsmult = 1.0
# HMOC parameters in case they are invoked
dceps = 1.e-5
nplane = 1
npl = 0
nph = 4
npmin = 0
npmax = 8
nlsink = nplane
npsink = nph

def gen_p01mt3dms(mt3d_ws, al, retardation, zero_order_decay, mixelm, silent=False):
    
    model_ws = mt3d_ws
    modelname_mf = 'p01_mf'
    
    # Instantiate the MODFLOW model
    mf = flopy.modflow.Modflow(modelname=modelname_mf, model_ws=model_ws, 
                               exe_name=exe_name_mf)
    
    # Instantiate discretization package
    # units: itmuni=4 (days), lenuni=2 (m)
    dis = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol,
                                   delr=delr, delc=delc, top=top, nstp=nstp,
                                   botm=botm, perlen=perlen, itmuni=4, lenuni=2)
    
    # Instantiate basic package
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
    
    # Instantiate layer property flow package
    lpf = flopy.modflow.ModflowLpf(mf, hk=hk, laytyp=laytyp)
    
    # Instantiate solver package
    pcg = flopy.modflow.ModflowPcg(mf)
    
    # Instantiate link mass transport package (for writing linker file for MT3D)
    lmt = flopy.modflow.ModflowLmt(mf)
    
    #-------------
    #
    #  Transport
    #
    #-------------
    
    # Instantiate the MT3DMS model
    modelname_mt = 'p01_mt'
    mt = flopy.mt3d.Mt3dms(modelname=modelname_mt, model_ws=model_ws,
                           exe_name=exe_name_mt, modflowmodel=mf)
    
    c0 = 1.
    icbund = np.ones((nlay, nrow, ncol), dtype=np.int)
    icbund[0, 0, 0] = -1
    sconc = np.zeros((nlay, nrow, ncol), dtype=np.float)
    sconc[0, 0, 0] = c0
    btn = flopy.mt3d.Mt3dBtn(mt, laycon=laytyp, icbund=icbund,
                             prsity=prsity, sconc=sconc, dt0=dt0, ifmtcn=1)
    
    # Instatiate the advection package
    adv = flopy.mt3d.Mt3dAdv(mt, mixelm=mixelm, dceps=dceps, nplane=nplane,
                             npl=npl, nph=nph, npmin=npmin, npmax=npmax,
                             nlsink=nlsink, npsink=npsink, percel=0.5)
    
    # Instantiate the dispersion package
    dsp = flopy.mt3d.Mt3dDsp(mt, al=al)
    
    # Set reactive variables and instantiate chemical reaction package
    if retardation==1.:
        isothm=0.
        ireact=0.
        rc1 = 0.
    else:
        isothm=1
        ireact=1
    
    if ireact==1:
        rc1 = zero_order_decay
    
    kd = (retardation - 1.) * prsity / rhob
    rct = flopy.mt3d.Mt3dRct(mt, isothm=isothm, ireact=ireact, igetsc=0, 
                             rhob=rhob, sp1=kd, rc1=rc1, rc2=rc1)
    
    # Instantiate the source/sink mixing package
    ssm = flopy.mt3d.Mt3dSsm(mt)
    
    # Instantiate the GCG solver in MT3DMS
    gcg = flopy.mt3d.Mt3dGcg(mt, mxiter=10)
    
    return (mf, mt, modelname_mt)

def gen_p01mf6(model_ws, al, retardation, lambda1, mixelm, xt3d=[False], 
           silent=False):
    
    # Name of associated model input files
    name = 'p01_mf6'
    
    # Instantiate the MODFLOW 6 simulation
    ws = model_ws
    exe_name = os.path.abspath(mf6exe)
    sim = flopy.mf6.MFSimulation(sim_name=name, version='mf6',
                                 exe_name=mf6exe,
                                 sim_ws=ws)
    
    from flopy.mf6.mfbase import VerbosityLevel
    sim.simulation_data.verbosity_level = VerbosityLevel.quiet
    sim.name_file.memory_print_option = 'all'
    
    # create tdis package
    tdis_rc = []
    tdis_rc.append((perlen, nstp, 1.0))
    tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS',
                                 nper=len(tdis_rc), perioddata=tdis_rc)
    
        # create gwf model
    gwfname = 'gwf_' + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True,
                               model_nam_file='{}.nam'.format(gwfname))
    
    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(sim, print_option='SUMMARY',
                                  outer_dvclose=hclose,
                                  outer_maximum=nouter,
                                  under_relaxation='NONE',
                                  inner_maximum=ninner,
                                  inner_dvclose=hclose, rcloserecord=rclose,
                                  linear_acceleration='CG',
                                  scaling_method='NONE',
                                  reordering_method='NONE',
                                  relaxation_factor=relax,
                                  filename='{}.ims'.format(gwfname))
    sim.register_ims_package(imsgwf, [gwf.name])
    
    # Instantiate discretization package
    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol,
                                  delr=delr, delc=delc,
                                  top=top, botm=botm,
                                  idomain=np.ones((nlay, nrow, ncol),
                                                  dtype=np.int),
                                  filename='{}.dis'.format(gwfname))
    
    # Instantiate initial conditions package
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt,
                                filename='{}.ic'.format(gwfname))
    
    # Instantiate node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False,
                                  icelltype=laytyp,
                                  k=hk,
                                  k33=hk, save_specific_discharge=True)
    
    # Instantiate constant head files
    chdspd = [[(0, 0, 0), h1], [(0, 0, ncol - 1), 0.0]]
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf,
                                                   maxbound=len(chdspd),
                                                   stress_period_data=chdspd,
                                                   save_flows=False,
                                                   pname='CHD-1')
    
    # output control
    oc = flopy.mf6.ModflowGwfoc(gwf,
                                budget_filerecord='{}.bud'.format(gwfname),
                                head_filerecord='{}.hds'.format(gwfname),
                                headprintrecord=[
                                    ('COLUMNS', 10, 'WIDTH', 15,
                                     'DIGITS', 6, 'GENERAL')],
                                saverecord=[('HEAD', 'LAST'),
                                            ('BUDGET', 'LAST')],
                                printrecord=[('HEAD', 'LAST'),
                                             ('BUDGET', 'LAST')])
    
    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False,
                                  icelltype=laytyp,
                                  k=hk,
                                  k33=hk, save_specific_discharge=True)
    
    # chd files
    chdspd = [[(0, 0, 0), h1], [(0, 0, ncol - 1), 0.0]]
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf,
                                                   maxbound=len(chdspd),
                                                   stress_period_data=chdspd,
                                                   save_flows=False,
                                                   pname='CHD-1')
    
    # output control
    oc = flopy.mf6.ModflowGwfoc(gwf,
                                budget_filerecord='{}.bud'.format(gwfname),
                                head_filerecord='{}.hds'.format(gwfname),
                                headprintrecord=[
                                    ('COLUMNS', 10, 'WIDTH', 15,
                                     'DIGITS', 6, 'GENERAL')],
                                saverecord=[('HEAD', 'LAST'),
                                            ('BUDGET', 'LAST')],
                                printrecord=[('HEAD', 'LAST'),
                                             ('BUDGET', 'LAST')])
    
    # create gwt model
    gwtname = 'gwt_' + name
    gwt = flopy.mf6.MFModel(sim, model_type='gwt6', modelname=gwtname,
                            model_nam_file='{}.nam'.format(gwtname))
    gwt.name_file.save_flows = True
    
    # create iterative model solution and register the gwt model with it
    imsgwt = flopy.mf6.ModflowIms(sim, print_option='SUMMARY',
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
    
    # Instantiate discretization package for transport
    dis = flopy.mf6.ModflowGwtdis(gwt, nlay=nlay, nrow=nrow, ncol=ncol,
                                  delr=delr, delc=delc,
                                  top=top, botm=botm,
                                  idomain=1,
                                  filename='{}.dis'.format(gwtname))
    
    # Instantiate initial conditions for transport problem
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.,
                                filename='{}.ic'.format(gwtname))
    
    # advection
    if mixelm == 0:
        scheme = 'UPSTREAM'
    elif mixelm == -1:
        scheme = 'TVD'
    else:
        raise Exception()
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme=scheme,
                                  filename='{}.adv'.format(gwtname))
    
    # dispersion
    if al != 0:
        dsp = flopy.mf6.ModflowGwtdsp(gwt, alh=al, ath1=al)
    
    # Set reactive variables and instantiate mass storage and transfer package
    if retardation != 1.:
        sorbtion = True
        kd = (retardation - 1.) * prsity / rhob  # prsity & rhob defined in 
                                                 # global variable section
    else:
        sorbtion = False
        kd = 1.
    
    if lambda1 != None:
        lmbda = lambda1
        fod = True
    else:
        lmbda = None
        fod = False
    
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=prsity,
                                  first_order_decay=fod,
                                  decay=lmbda, 
                                  sorbtion=sorbtion, 
                                  bulk_density=rhob, 
                                  distcoef=kd)
    
    # Instantiate constant concentration
    c0 = 1.
    cncspd = [[(0, 0, 0), c0]]
    cnc = flopy.mf6.ModflowGwtcnc(gwt, maxbound=len(cncspd),
                                  stress_period_data=cncspd,
                                  save_flows=False,
                                  pname='CNC-1')
    
    # Instantiate source/sink mixing package
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=[[]],
                                  filename='{}.ssm'.format(gwtname))
    
    # Instantiate output control for transport model
    oc = flopy.mf6.ModflowGwtoc(gwt,
                                budget_filerecord='{}.cbc'.format(gwtname),
                                concentration_filerecord='{}.ucn'.format(
                                    gwtname),
                                concentrationprintrecord=[
                                    ('COLUMNS', 10, 'WIDTH', 15,
                                     'DIGITS', 6, 'GENERAL')],
                                saverecord=[('CONCENTRATION', 'LAST'),
                                            ('BUDGET', 'LAST')],
                                printrecord=[('CONCENTRATION', 'LAST'),
                                             ('BUDGET', 'LAST')])
    
    # Setup GWF-GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(sim, exgtype='GWF6-GWT6',
                                     exgmnamea=gwfname, exgmnameb=gwtname,
                                     filename='{}.gwfgwt'.format(name))
    
    return (sim, gwtname)

def clean_model_output(mt3d_model, mf6_model, gwtname):
    # Get the paths to the model output
    mt3d_out_path = mt3d_model.model_ws
    mf6_out_path = mf6_model.simulation_data.mfpath.get_sim_path() 
    
    # Clear previous model output (to ensure reading of latest model output)
    fname_mt3d = os.path.join(os.path.abspath(mt3d_out_path), 'MT3D001.UCN')
    if os.path.isfile(fname_mt3d):
        os.remove(fname_mt3d)
    
    # Clear previous MF6 output (to ensure reading of latest model output)
    fname_mf6 = os.path.join(os.path.abspath(mf6_out_path), gwtname + '.ucn')
    if os.path.isfile(fname_mf6):
        os.remove(fname_mf6)

def get_model_output(mt3d, mf6, gwtname):
    # Get the paths to the model output
    mt3d_out_path = mt3d.model_ws
    mf6_out_path = mf6.simulation_data.mfpath.get_sim_path()
    
    # Get the MT3DMS concentration output
    fname_mt3d = os.path.join(mt3d_out_path, 'MT3D001.UCN')
    ucnobj_mt3d = flopy.utils.UcnFile(fname_mt3d)
    times_mt3d = ucnobj_mt3d.get_times()
    conc_mt3d = ucnobj_mt3d.get_alldata()
    
    # Get the MF6 concentration output
    fname_mf6 = os.path.join(mf6_out_path, gwtname + '.ucn')
    ucnobj_mf6 = flopy.utils.HeadFile(fname_mf6, precision='double',
                                      text='CONCENTRATION')
    times_mf6 = ucnobj_mf6.get_times()
    conc_mf6 = ucnobj_mf6.get_alldata()
    
    return conc_mt3d, conc_mf6 

def make_figure(mt3d, mf6, gwtname, ax=None, figname=None):
    
    ncol = mt3d.ncol
    conc_mt3d, conc_mf6 = get_model_output(mt3d, mf6, gwtname)
    
    if ax is None:
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(1, 1, 1)
    
    ax.plot(np.linspace(0, l, ncol), conc_mt3d[0,0,0,:], color='k', label='MT3DMS')
    ax.plot(np.linspace(0, l, ncol), conc_mf6[0,0,0,:], '^', color='b', label='MF6')
    ax.set_ylim(0, 1.2)
    ax.set_xlim(0, 1000)
    ax.set_xlabel('DISTANCE, IN M')
    ax.set_ylabel('CONCENTRATION')
    plt.title('CONCENTRATION PROFILE @ TIME = 1,000 DAYS'.format(3))
    ax.legend()
    
    if figname is not None:
        fig.savefig(fname=figname)
    
    return ax

def test_p01():
    main()
    
def main(write_models=True, run_models=True, plot_results=True):
    
    # General options
    silent = False
    
    times = 2000
    
    # Define all models to be used regardless of which options the user selects
    # p01 specifics
    al = [0., 10., 10., 10.]           # Longitudinal dispersion
    retardation = [1.0, 1.0, 5.0, 5.0] # Original problem specifies retardation
    lambda1 = [0., 0., 0., 0.002]      # Decay term
    mixelm = [-1, -1, -1, 0]           # Use TVD since MF6 doesn't have a MOC
                                       # method at this point, which is best 
                                       # for advection only.
    xt3d=[False, True, False, False]   # For MF6
    testgroup = ['mt3dms_p01_a', 'mt3dms_p01_b', 'mt3dms_p01_c', 'mt3dms_p01_d']
    
    mt3d_test_models = []
    mf6_test_models = []
    
    for i in range(len(al)):
        mt3d_ws = os.path.join(testdir, testgroup[i], 'mt3d')
        mf6_ws = os.path.join(testdir, testgroup[i], 'mf6')
        
        mt3d_test_models.append(gen_p01mt3dms(mt3d_ws=mt3d_ws, al=al[i], 
                                              retardation=retardation[i], 
                                              zero_order_decay=lambda1[i], 
                                              mixelm=mixelm[i], silent=True))
        
        mf6_test_models.append(gen_p01mf6(mf6_ws, al=al[i], 
                                          retardation=retardation[i], 
                                          lambda1=lambda1[i], mixelm=mixelm[i],
                                          xt3d=xt3d[i],
                                          silent=True))
    
    if write_models:
        # Name all models something like "ex_..."	
        print("Writing model input")
        
        for mt3d_model, mf6_model in zip(mt3d_test_models, mf6_test_models):
            
            mf2k5 = mt3d_model[0]
            mt3d  = mt3d_model[1]
            
            # Write the MF2k5/MT3DMS input
            mf2k5.write_input()
            mt3d.write_input()
            
            # Next write the MF6 flow & transport input
            mf6 = mf6_model[0]
            mf6.write_simulation()
    
    if run_models:
        
        print("Running MF2k5/MT3DMS and MF6 transport models")
        
        for mt3d_model, mf6_model in zip(mt3d_test_models, mf6_test_models):
            
            # If it exists, clean-out the existing model output to ensure old
            # output isn't accidently read if the model fails to completely run
            clean_model_output(mt3d_model[1], mf6_model[0], mf6_model[1])
            
            # Run the MF2k5/MT3DMS models
            mt3d_model[0].run_model(silent=silent)  # 0th position is the modflow model
            mt3d_model[1].run_model(silent=silent)  # 1st position is the mt3d model
            
            # Now run the MF6 model
            success, buff = mf6_model[0].run_simulation(silent=silent, report=True)
            if not success:
                print(buff)
            assert(success), 'MF6 bombed'
    
    if plot_results:
        print("Plotting up MF2k5/MT3DMS comparison with MF6 transport")
        
        for i, (mt3d, mf6) in enumerate(zip(mt3d_test_models, mf6_test_models)):
            gwtname = mf6[1]
            figname = os.path.join('..','figures','p01mt3d-f' + str(i+1)+'.pdf')
            ax = make_figure(mt3d[1], mf6[0], gwtname, figname=figname)
        
        
if __name__ == "__main__":
    
    # Determine which mode(s) to run in
    write_models = True
    run_models = True
    plot_results = True
    
    for idx, arg in enumerate(sys.argv):
        if '--no_write' in arg.lower():
            write_models = False
        if '--no_run' in arg.lower():
            run_models = False
        if '--no_plot' in arg.lower():
            plot_results = False

    main(write_models, run_models, plot_results)