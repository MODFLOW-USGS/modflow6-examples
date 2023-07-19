# ## USG1DISU example
#
# This example, ex-gwf-radial, shows how the MODFLOW 6 DISU Package
# can be used to simulate an axisymmetric radial model.
#
# The example corresponds to the first example described in:
#    Bedekar, V., Scantlebury, L., and Panday, S. (2019).
#       Axisymmetric Modeling Using MODFLOW-USG.Groundwater, 57(5), 772-777.
#
# And the numerical result is compared against the analytical solution
# presented in Equation 17 of
#    Neuman, S. P. (1974). Effect of partial penetration on flow in
#    unconfined aquifers considering delayed gravity response.
#    Water resources research, 10(2), 303-312

# Imports

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import flopy

# Append to system path to include the common subdirectory

sys.path.append(os.path.join("..", "common"))

# import common functionality

import config
from figspecs import USGSFigure

from get_disu_radial_kwargs import get_disu_radial_kwargs
from neuman1974_soln import RadialUnconfinedDrawdown

# Utility function to return DISU node for given radial band and layer


def get_radial_node(rad, lay, nradial):
    """
    Given nradial dimension (bands per layer),
    returns the 0-based disu node number for given 0-based radial
    band and layer

    Parameters
    ----------
    rad : int
        radial band number (0 to nradial-1)
    lay : float or ndarray
        layer number (0 to nlay-1)
    nradial : int
        total number of radial bands

    Returns
    -------
    result : int
        0-based disu node number located at rad and lay

    """

    return nradial * lay + rad


def get_radius_lay_from_node(node, nradial):
    """
    Given nradial dimension (bands per layer),
    returns 0-based layer and radial band indices for given disu node number

    Parameters
    ----------
    node : int
        disu node number
    nradial : int
        total number of radial bands

    Returns
    -------
    result : int
        0-based disu node number located at rad and lay

    """
    #
    lay = node // nradial
    rad = node - (lay * nradial)
    return rad, lay


# Run Analytical Model - Very slow
# If True, solves the Neuman 1974 analytical model (very slow)
# else uses stored results from solving the Neuman 1974 analytical model

solve_analytical_solution = False


# Set default figure properties

figure_size = (6, 6)

# Base simulation and model name and workspace

ws = config.base_ws

# Simulation name

sim_name = "ex-gwf-rad-disu"

# Model units

length_units = "feet"
time_units = "days"

# Table Model Parameters

nper = 1  # Number of periods
_ = 24  # Number of time steps
_ = "10"  # Simulation total time ($day$)

nlay = 25  # Number of layers
nradial = 22  # Number of radial direction cells (radial bands)

initial_head = 50.0  # Initial water table elevation ($ft$)

surface_elevation = 50.0  # Top of the radial model ($ft$)
_ = 0.0  # Base of the radial model ($ft$)
layer_thickness = 2.0  # Thickness of each radial layer ($ft$)
_ = "0.25 to 2000"  # Outer radius of each radial band ($ft$)

k11 = 20.0  # Horizontal hydraulic conductivity ($ft/day$)
k33 = 20.0  # Vertical hydraulic conductivity ($ft/day$)

ss = 1.0e-5  # Specific storage ($1/day$)
sy = 0.1  # Specific yield (unitless)

_ = "0.0 to 10"  # Well screen elevation ($ft$)
_ = "1"  # Well radial band location (unitless)
_ = "-4000.0"  # Well pumping rate ($ft^3/day$)

_ = "40"  # Observation distance from well ($ft$)
_ = "1"  # ``Top'' observation elevation ($ft$)
_ = "25"  # ``Middle'' observation depth ($ft$)
_ = "49"  # ``Bottom'' observation depth ($ft$)


# Outer Radius for each radial band

radius_outer = [
    0.25,
    0.75,
    1.5,
    2.5,
    4.0,
    6.0,
    9.0,
    13.0,
    18.0,
    23.0,
    33.0,
    47.0,
    65.0,
    90.0,
    140.0,
    200.0,
    300.0,
    400.0,
    600.0,
    1000.0,
    1500.0,
    2000.0,
]  # Outer radius of each radial band ($ft$)

# Well boundary conditions
# Well must be located on central radial band (rad = 0)
# and have a contiguous screen interval and constant pumping rate.
# This example has the well screen interval from
# layer 20 to 24 (zero-based index)
wel_spd = {
    sp: [
        [(get_radial_node(0, lay, nradial),), -800.0] for lay in range(20, 25)
    ]
    for sp in range(nper)
}

# Static temporal data used by TDIS file
# Simulation has 1 ten-day stress period with 24 time steps.
# The multiplier for the length of successive time steps is 1.62

tdis_ds = ((10.0, 24, 1.62),)

# Setup observation location and times

obslist = [
    ["h_top", "head", (get_radial_node(11, 0, nradial),)],
    ["h_mid", "head", (get_radial_node(11, (nlay - 1) // 2, nradial),)],
    ["h_bot", "head", (get_radial_node(11, nlay - 1, nradial),)],
]
obsdict = {"{}.obs.head.csv".format(sim_name): obslist}

# Solver parameters

nouter = 500
ninner = 300
hclose = 1e-4
rclose = 1e-4


# ### Functions to build, write, run, and plot the MODFLOW 6 Axisymmetric Model
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_model(name):
    if config.buildModel:
        sim_ws = os.path.join(ws, name)
        sim = flopy.mf6.MFSimulation(
            sim_name=name, sim_ws=sim_ws, exe_name=config.mf6_exe
        )
        flopy.mf6.ModflowTdis(
            sim, nper=nper, perioddata=tdis_ds, time_units=time_units
        )
        flopy.mf6.ModflowIms(
            sim,
            print_option="summary",
            complexity="complex",
            outer_maximum=nouter,
            outer_dvclose=hclose,
            inner_maximum=ninner,
            inner_dvclose=hclose,
        )

        gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

        disukwargs = get_disu_radial_kwargs(
            nlay,
            nradial,
            radius_outer,
            surface_elevation,
            layer_thickness,
            get_vertex=True,
        )

        disu = flopy.mf6.ModflowGwfdisu(
            gwf, length_units=length_units, **disukwargs
        )

        npf = flopy.mf6.ModflowGwfnpf(
            gwf,
            k=k11,
            k33=k33,
            save_flows=True,
            save_specific_discharge=True,
        )

        flopy.mf6.ModflowGwfsto(
            gwf,
            iconvert=1,
            sy=sy,
            ss=ss,
            save_flows=True,
        )

        flopy.mf6.ModflowGwfic(gwf, strt=initial_head)

        flopy.mf6.ModflowGwfwel(
            gwf, stress_period_data=wel_spd, save_flows=True
        )

        flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=f"{name}.cbc",
            head_filerecord=f"{name}.hds",
            headprintrecord=[
                ("COLUMNS", nradial, "WIDTH", 15, "DIGITS", 6, "GENERAL")
            ],
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            filename=f"{name}.oc",
        )

        flopy.mf6.ModflowUtlobs(gwf, print_input=False, continuous=obsdict)
        return sim
    return None


# Function to write model files


def write_model(sim, silent=True):
    if config.writeModel:
        sim.write_simulation(silent=silent)


# Function to run the Axisymmetric model.
# True is returned if the model runs successfully.


@config.timeit
def run_model(sim, silent=True):
    success = True
    if config.runModel:
        success, buff = sim.run_simulation(silent=silent, report=True)
        if not success:
            print("\n".join(buff))

    return success


# Function to solve Axisymmetric model using analytical equation.


def solve_analytical(obs2ana, times=None, no_solve=False):
    # obs2ana = {obsdict[file][0] : analytical_name}
    disukwargs = get_disu_radial_kwargs(
        nlay, nradial, radius_outer, surface_elevation, layer_thickness
    )
    model_bottom = disukwargs["bot"][get_radial_node(0, nlay - 1, nradial)]
    sat_thick = initial_head - model_bottom

    key = next(iter(wel_spd))
    nodes = []
    rates = []
    for nod, rat in wel_spd[key]:
        nodes.append(nod[0])
        rates.append(rat)
    nodes.sort()
    well_top = disukwargs["top"][nodes[0]]
    well_bot = disukwargs["bot"][nodes[-1]]
    pump = abs(sum(rates))

    ana_model = RadialUnconfinedDrawdown(
        bottom_elevation=model_bottom,
        hydraulic_conductivity_radial=k11,
        hydraulic_conductivity_vertical=k33,
        specific_storage=ss,
        specific_yield=sy,
        well_screen_elevation_top=well_top,
        well_screen_elevation_bottom=well_bot,
        saturated_thickness=sat_thick,
    )

    build_times = times is None
    if build_times:
        totim = 0.0
        for pertim in tdis_ds:
            totim += pertim[0]
        times_sy_base = np.logspace(-3, max([np.log10(totim), 2]), 100)

    analytical = {}
    prop = {}
    if not no_solve:
        print("Solving Analytical Model (Very Slow)")
    for file in obsdict:
        for row in obsdict[file]:
            obs = row[0].upper()
            if obs not in obs2ana:
                continue
            ana = obs2ana[obs]
            nod = row[2][0]
            obs_top = disukwargs["top"][nod]
            obs_bot = disukwargs["bot"][nod]

            rad, lay = get_radius_lay_from_node(nod, nradial)

            if lay == 0:
                # Uppermost layer has obs elevation at top,
                # otherwise cell center
                obs_el = obs_top
            else:
                obs_el = 0.5 * (obs_top + obs_bot)

            if rad == 0:
                obs_rad = 0.0
            else:
                # radius_outer[r-1] + 0.5*(radius_outer[r] - radius_outer[r-1])
                obs_rad = 0.5 * (radius_outer[rad - 1] + radius_outer[rad])

            if build_times:
                times_sy = times_sy_base
                times = [
                    ty * sy * obs_rad * obs_rad / (k11 * sat_thick)
                    for ty in times_sy_base
                ]
            else:
                times_sy = [
                    ty * k11 * sat_thick / (sy * obs_rad * obs_rad)
                    for ty in times
                ]

            times_ss = [ty * k11 / (ss * obs_rad * obs_rad) for ty in times]

            if not no_solve:
                print(f"Solving {ana}")
                analytical[ana] = ana_model.drawdown_times(
                    pump,
                    times,
                    obs_rad,
                    obs_el,
                    sumrtol=1.0e-6,
                    u_n_rtol=1.0e-5,
                )
            prop[ana] = [
                times,
                times_sy,
                times_ss,
                pump,
                obs_rad,
                sat_thick,
                model_bottom,
            ]

    return analytical, prop


# Function to plot the Axisymmetric model results.


def plot_ts(sim, verbose=False, solve_analytical_solution=False):
    pi = 3.141592653589793
    gwf = sim.get_model(sim_name)
    obs_csv_name = gwf.obs.output.obs_names[0]
    obs_csv_file = gwf.obs.output.obs(f=obs_csv_name)

    tsdata = obs_csv_file.data
    fmt = {
        "H_TOP": "og",
        "H_MID": "or",
        "H_BOT": "ob",
        "a_top": "-g",
        "a_mid": "-r",
        "a_bot": "-b",
    }
    obsnames = {
        "H_TOP": "MF6 (Top)",
        "H_MID": "MF6 (Middle)",
        "H_BOT": "MF6 (Bottom)",
        "a_top": "Analytical (Top)",
        "a_mid": "Analytical (Middle)",
        "a_bot": "Analytical (Bottom)",
    }

    obs2ana = {"H_TOP": "a_top", "H_MID": "a_mid", "H_BOT": "a_bot"}

    if solve_analytical_solution:
        analytical, ana_prop = solve_analytical(obs2ana)
        analytical_time = []
    else:
        analytical, ana_prop = solve_analytical(obs2ana, no_solve=True)
        analytical_time = [
            0.00016,
            0.000179732,
            0.000201897,
            0.000226796,
            0.000254765,
            0.000286184,
            0.000321477,
            0.000361123,
            0.000405658,
            0.000455686,
            0.000511883,
            0.00057501,
            0.000645923,
            0.000725581,
            0.000815062,
            0.000915579,
            0.001028492,
            0.001155329,
            0.001297809,
            0.00145786,
            0.00163765,
            0.001839611,
            0.002066479,
            0.002321326,
            0.002607601,
            0.002929181,
            0.00329042,
            0.003696208,
            0.004152039,
            0.004664085,
            0.005239279,
            0.005885408,
            0.00661122,
            0.007426542,
            0.008342413,
            0.009371233,
            0.010526932,
            0.011825155,
            0.013283481,
            0.014921654,
            0.016761852,
            0.018828991,
            0.021151058,
            0.023759492,
            0.026689609,
            0.029981079,
            0.033678466,
            0.037831831,
            0.042497405,
            0.047738356,
            0.053625642,
            0.060238973,
            0.067667886,
            0.076012963,
            0.085387188,
            0.09591748,
            0.107746411,
            0.121034132,
            0.13596055,
            0.152727753,
            0.171562756,
            0.192720566,
            0.216487644,
            0.243185773,
            0.273176424,
            0.306865642,
            0.34470955,
            0.387220522,
            0.434974119,
            0.488616881,
            0.548875086,
            0.616564575,
            0.692601805,
            0.778016253,
            0.873964355,
            0.981745164,
            1.102817937,
            1.238821892,
            1.391598404,
            1.563215932,
            1.755998025,
            1.972554783,
            2.215818194,
            2.48908183,
            2.79604544,
            3.14086504,
            3.528209184,
            3.96332217,
            4.452095044,
            5.00114536,
            5.617906775,
            6.310729695,
            7.088994332,
            7.963237703,
            8.945296292,
            10.04846631,
            11.2876837,
            12.67972637,
            14.24344137,
            16,
        ]

        analytical["a_top"] = [
            9.14e-06,
            1.48e-05,
            2.32e-05,
            3.51e-05,
            5.16e-05,
            7.39e-05,
            0.000103386,
            0.000141564,
            0.000190115,
            0.000250872,
            0.000325813,
            0.000417056,
            0.00052686,
            0.000657611,
            0.000811829,
            0.000992179,
            0.001201482,
            0.001442755,
            0.001719253,
            0.002034538,
            0.002392557,
            0.002797739,
            0.003255095,
            0.003770324,
            0.004349925,
            0.005001299,
            0.005732852,
            0.006554096,
            0.007475752,
            0.008509865,
            0.009669933,
            0.010971051,
            0.01243007,
            0.014065786,
            0.015899134,
            0.017953402,
            0.020254465,
            0.022831026,
            0.025714866,
            0.028941104,
            0.032548456,
            0.03657948,
            0.041080814,
            0.046103371,
            0.051702489,
            0.057938,
            0.064874212,
            0.072579742,
            0.081127182,
            0.090592554,
            0.101054493,
            0.112593132,
            0.125288641,
            0.139219391,
            0.154459742,
            0.171077492,
            0.189131034,
            0.208666373,
            0.229714163,
            0.252287007,
            0.276377285,
            0.301955794,
            0.328971456,
            0.357352252,
            0.387007461,
            0.417831084,
            0.449706211,
            0.482509951,
            0.516118451,
            0.550411552,
            0.585276682,
            0.620611726,
            0.656326766,
            0.692344745,
            0.72860122,
            0.765043446,
            0.801629049,
            0.838324511,
            0.875103648,
            0.911946208,
            0.94883664,
            0.985763066,
            1.022716447,
            1.059689924,
            1.0966783,
            1.133677645,
            1.170684993,
            1.207698108,
            1.24471531,
            1.281735341,
            1.318757264,
            1.355780385,
            1.392804192,
            1.429828316,
            1.466852493,
            1.503876535,
            1.540900317,
            1.577923757,
            1.614946805,
            1.651969438,
        ]

        analytical["a_mid"] = [
            0.042573248,
            0.053215048,
            0.065047371,
            0.077900907,
            0.091562277,
            0.10578645,
            0.120308736,
            0.134861006,
            0.149179898,
            0.163017203,
            0.176149026,
            0.188382628,
            0.199562503,
            0.209575346,
            0.218353885,
            0.22587733,
            0.232171831,
            0.237306589,
            0.241387585,
            0.244548545,
            0.246940521,
            0.24872049,
            0.250040497,
            0.25103848,
            0.251831919,
            0.252514802,
            0.253157938,
            0.253811816,
            0.254511215,
            0.255280058,
            0.256135833,
            0.257093057,
            0.258165369,
            0.259366956,
            0.260713168,
            0.262220999,
            0.263909262,
            0.265798794,
            0.267912645,
            0.270276262,
            0.272917675,
            0.27586768,
            0.279160018,
            0.28283153,
            0.286922294,
            0.291475727,
            0.296538626,
            0.302161152,
            0.308396727,
            0.31530182,
            0.322935596,
            0.331359479,
            0.340636352,
            0.350829851,
            0.362003212,
            0.374217985,
            0.387532593,
            0.402000693,
            0.417669407,
            0.434577538,
            0.452753827,
            0.472215504,
            0.492966953,
            0.514999008,
            0.538288679,
            0.562799606,
            0.588483025,
            0.615279515,
            0.643120962,
            0.671933085,
            0.701637983,
            0.732156526,
            0.763410706,
            0.795325415,
            0.827829825,
            0.860858334,
            0.894351053,
            0.928253957,
            0.962518747,
            0.997102451,
            1.031967356,
            1.067080013,
            1.102411044,
            1.137934722,
            1.17362841,
            1.209472238,
            1.245448747,
            1.281542585,
            1.317740238,
            1.354029805,
            1.390400789,
            1.426843932,
            1.463351049,
            1.499914938,
            1.53652918,
            1.573188127,
            1.609886766,
            1.64662066,
            1.683385875,
            1.720178921,
        ]

        analytical["a_bot"] = [
            0.086684154,
            0.103649206,
            0.12188172,
            0.141163805,
            0.16124535,
            0.181846061,
            0.202667182,
            0.223392774,
            0.243699578,
            0.263275092,
            0.281824849,
            0.299088528,
            0.31485153,
            0.328955429,
            0.341304221,
            0.351868481,
            0.360684432,
            0.367848506,
            0.373508853,
            0.377853369,
            0.381094083,
            0.383450902,
            0.385136877,
            0.386345028,
            0.387238663,
            0.387947536,
            0.388568938,
            0.389170336,
            0.389796684,
            0.390476869,
            0.39123089,
            0.392073081,
            0.393015962,
            0.39407278,
            0.395256673,
            0.396583156,
            0.398068337,
            0.399731151,
            0.401591476,
            0.403672105,
            0.405997893,
            0.408596189,
            0.411497001,
            0.414733164,
            0.418340479,
            0.422357833,
            0.426826999,
            0.431793796,
            0.437306073,
            0.443415589,
            0.450176656,
            0.457646089,
            0.46588273,
            0.47494648,
            0.484898814,
            0.495799471,
            0.507707617,
            0.520679174,
            0.534765662,
            0.550012634,
            0.566458073,
            0.584130854,
            0.60304937,
            0.623219948,
            0.644638018,
            0.667285046,
            0.691130674,
            0.716132499,
            0.742239678,
            0.76939134,
            0.797520861,
            0.826557513,
            0.85642935,
            0.887062398,
            0.918386519,
            0.950334816,
            0.982842103,
            1.015851071,
            1.049306941,
            1.083160734,
            1.117368294,
            1.151890127,
            1.186690934,
            1.221739246,
            1.257007172,
            1.292471077,
            1.32810677,
            1.363895877,
            1.399822391,
            1.435868577,
            1.472023401,
            1.508273605,
            1.544607161,
            1.581017341,
            1.617494414,
            1.654030939,
            1.690620311,
            1.72725511,
            1.763933202,
            1.800648435,
        ]

    fs = USGSFigure(figure_type="graph", verbose=verbose)

    obs_fig = "obs-head"
    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_subplot()
    ax.set_xlabel("time (d)")
    ax.set_ylabel("head (ft)")
    for name in tsdata.dtype.names[1:]:
        ax.plot(
            tsdata["totim"],
            tsdata[name],
            fmt[name],
            label=obsnames[name],
            markerfacecolor="none",
        )
        # , markersize=3

    for name in analytical:
        n = len(analytical[name])
        if solve_analytical_solution:
            ana_times = ana_prop[name][0]
        else:
            ana_times = analytical_time

        ax.plot(
            ana_times[:n],
            [50.0 - h for h in analytical[name]],
            fmt[name],
            label=obsnames[name],
        )

    fs.graph_legend(ax)

    fig.tight_layout()

    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-{}{}".format(sim_name, obs_fig, config.figure_ext),
        )
        fig.savefig(fpth)

    obs_fig = "obs-dimensionless"
    fig = plt.figure(figsize=(5, 3))
    fig.tight_layout()
    ax = fig.add_subplot()
    ax.set_xlim(0.001, 100.0)
    ax.set_ylim(0.001, 100.0)
    ax.grid(visible=True, which="major", axis="both")
    ax.set_ylabel("Dimensionless Drawdown, $s_d$")
    ax.set_xlabel("Dimensionless Time, $t_y$")
    for name in tsdata.dtype.names[1:]:
        q = ana_prop[obs2ana[name]][3]
        r = ana_prop[obs2ana[name]][4]
        b = ana_prop[obs2ana[name]][5]
        ax.loglog(
            [k11 * b * ts / (sy * r * r) for ts in tsdata["totim"]],
            [4 * pi * k11 * b * (initial_head - h) / q for h in tsdata[name]],
            fmt[name],
            label=obsnames[name],
            markerfacecolor="none",
        )

    for name in analytical:
        q = ana_prop[name][3]
        b = ana_prop[name][5]  # [pump, radius, sat_thick, model_bottom]
        if solve_analytical_solution:
            ana_times = ana_prop[name][0]
        else:
            ana_times = analytical_time

        n = len(analytical[name])
        time_sy = [k11 * b * ts / (sy * r * r) for ts in ana_times[:n]]
        ana = [4 * pi * k11 * b * s / q for s in analytical[name]]
        ax.plot(time_sy, ana, fmt[name], label=obsnames[name])

    fs.graph_legend(ax)

    fig.tight_layout()

    if config.plotSave:
        fpth = os.path.join(
            "..",
            "figures",
            "{}-{}{}".format(sim_name, obs_fig, config.figure_ext),
        )
        fig.savefig(fpth)


# Function to plot the model radial bands.


def plot_grid(verbose=False):
    fs = USGSFigure(figure_type="map", verbose=verbose)

    # Print all radial bands
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6.4, 3.1))
    # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.5))
    ax = axs[0]

    max_rad = radius_outer[-1]
    max_rad = max_rad + (max_rad * 0.1)
    ax.set_xlim(-max_rad, max_rad)
    ax.set_ylim(-max_rad, max_rad)
    ax.set_aspect("equal", adjustable="box")

    circle_center = (0.0, 0.0)
    for r in radius_outer:
        circle = Circle(circle_center, r, color="black", fill=False, lw=0.3)
        ax.add_artist(circle)

    ax.set_xlabel("x-position (ft)")
    ax.set_ylabel("y-position (ft)")
    ax.annotate(
        "A",
        (-0.11, 1.02),
        xycoords="axes fraction",
        fontweight="black",
        fontsize="xx-large",
    )

    # Print first 5 radial bands
    nband = 5
    ax = axs[1]

    radius_subset = radius_outer[:nband]
    max_rad = radius_subset[-1]
    max_rad = max_rad + (max_rad * 0.3)

    ax.set_xlim(-max_rad, max_rad)
    ax.set_ylim(-max_rad, max_rad)
    ax.set_aspect("equal", adjustable="box")

    circle_center = (0.0, 0.0)

    r = radius_subset[0]
    circle = Circle(circle_center, r, color="red", label="Well")
    ax.add_artist(circle)
    for r in radius_subset:
        circle = Circle(circle_center, r, color="black", lw=1, fill=False)
        ax.add_artist(circle)

    ax.set_xlabel("x-position (ft)")
    ax.set_ylabel("y-position (ft)")

    ax.annotate(
        "B",
        (-0.06, 1.02),
        xycoords="axes fraction",
        fontweight="black",
        fontsize="xx-large",
    )

    fs.graph_legend(ax)

    fig.tight_layout()

    # save figure
    if config.plotSave:
        fpth = os.path.join(
            "..", "figures", "{}-grid{}".format(sim_name, config.figure_ext)
        )
        fig.savefig(fpth)
    return


# Function to plot the model results.


def plot_results(silent=True):
    if not config.plotModel:
        return

    if silent:
        verbosity_level = 0
    else:
        verbosity_level = 1

    sim_ws = os.path.join(ws, sim_name)
    sim = flopy.mf6.MFSimulation.load(
        sim_name=sim_name, sim_ws=sim_ws, verbosity_level=verbosity_level
    )

    verbose = not silent

    if config.plotModel:
        plot_grid(verbose)
        plot_ts(
            sim, verbose, solve_analytical_solution=solve_analytical_solution
        )
    return


# Function that wraps all of the steps for the Axisymmetric model
#
# 1. build_model,
# 2. write_model,
# 3. run_model, and
# 4. plot_results.
#


def simulation(silent=True):
    # key = list(parameters.keys())[idx]
    # params = parameters[key].copy()

    sim = build_model(sim_name)

    write_model(sim, silent=silent)

    success = run_model(sim, silent=silent)
    assert success, "could not run...{}".format(sim_name)


# nosetest - exclude block from this nosetest to the next nosetest
def test_and_plot():
    simulation(silent=False)
    plot_results(silent=False)
    return


# nosetest end


if __name__ == "__main__":
    # ### Axisymmetric Example

    # MF6 Axisymmetric Model
    simulation()

    # Solve analytical and plot results with MF6 results
    plot_results()
