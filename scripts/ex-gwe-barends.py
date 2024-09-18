# +
# An analytical solution for this example was originally published in
#
#  Barends, F.B.J., 2010. Complete Solution for Transient Heat Transport in
#    Porous Media, Following Lauwerier's Concept. Society of Petroleum '
#    Engineers, Annual Technical Conference and Exhibition, Florence, Italy,
#    19–22 September 2010.
#    https://doi.org/10.2118/134670-MS
#
# Below is a diagram of the model cell numbering, which is helpful for this
# DISU setup since the groundwater reservoir is flowing in a 1D, left-to-right
# manner with no flow in the overburden.  Additionally, heat "bleeds" into
# the overburden from the flowing reservoir.  Heat bleeds in a 1D upward
# manner and does not conduct laterally once in the overburden.  Thus, the
# connections in the DISU grid are such that cells in the overburden are not
# connected to their lateral neighbors, but only the cells above and below.
#
#           +-------+-------+-------+     +-------+-------+-------+          \
# Cell IDs  |   0   |   1   |   2   | ... |  997  |  998  |  999  | Layer 1   |
# (0-based) +-------+-------+-------+     +-------+-------+-------+           |
#           | 1,000 | 1,001 | 1,002 | ... | 1,997 | 1,998 | 1,999 | Layer 2   |
#           +-------+-------+-------+     +-------+-------+-------+           | -- "Overburden"
#           | 2,000 | 2,001 | 2,002 | ... | 2,997 | 2,998 | 2,999 | Layer 3   |
#           +-------+-------+-------+     +-------+-------+-------+           .
#               .       .       .             .       .       .               .
#               .       .       .             .       .       .               .
#               .       .       .             .       .       .              /
#           +-------+-------+-------+     +-------+-------+-------+          \
#           |       |       |       |     |       |       |       | Layer 1   |
#           |       |       |       |     |       |       |       |           |
#           |100,000|100,001|100,002| ... |100,997|100,998|100,999| Layer 2   | -- "Groundwater reservoir"
#           |       |       |       |     |       |       |       |           |
#           |       |       |       |     |       |       |       | Layer 3   |
#           +-------+-------+-------+     +-------+-------+-------+          /
#                       ----> gw flow direction ---->
# -

# +
import math
import os
import pathlib as pl

import flopy
import flopy.discretization as fgrid
import git
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from modflow_devtools.misc import get_env, timed
from scipy.integrate import quad

# Example name and base workspace
sim_name = "ex-gwe-barends"
try:
    root = pl.Path(git.Repo(".", search_parent_directories=True).working_dir)
except:
    root = None

workspace = root / "examples" if root else pl.Path.cwd()
figs_path = root / "figures" if root else pl.Path.cwd()

# Settings from environment variables
write = get_env("WRITE", True)
run = get_env("RUN", True)
plot = get_env("PLOT", True)
plot_show = get_env("PLOT_SHOW", True)
plot_save = get_env("PLOT_SAVE", True)
# -

# ### Define parameters
#
# Define model units, parameters and other settings.

# +
# Model units
length_units = "meters"
time_units = "seconds"

nrow = 1
ncol = 250
L = 1000  # Length of simulation in X direction ($m$)
delc = 1.0  # Width along the column ($m$)
nlay_res = 1  # Number of layers in the groundwater reservoir ($-$)
nlay_overburden = 100  # Number of layer representing the overburden ($-$)
nlay = nlay_overburden + nlay_res  # Total number of layers ($-$)
top_res = 100.0  # Elevation of the top of the groundwater reservoir ($m$)
H = 100.0  # Reservoir height ($m$)
thk = 1.0  # Unit thickness "into the page" ($m$)
k11 = 9.80670e-7  # Hydraulic conductivity of the groundwater reservoir($m/d$)
k33 = 1e-50  # Vertical hydraulic conductivity of the entire domain ($m/s$)
perlen = 86400.0 * 365  # Period length ($seconds$)
nstp = 1  # Number of time steps per stress period ($-$)
tsmult = 1  # Time step multiplier ($-$)
nper = 300  # Number of simulated stress periods ($-$)
prsity = 0.25  # Porosity ($-$)
scheme = "TVD"  # Advection solution scheme ($-$)
ktw = 0.6  # Thermal conductivity of the fluid ($\dfrac{W}{m \cdot ^{\circ}C}$)
kts = 0.06667  # Thermal conductivity of the aquifer matrix (which is also the dry overburden material) ($\dfrac{W}{m \cdot ^{\circ}C}$)
rhow = 1000.0  # Density of water ($\frac{kg}{m^3}$) # 3282.296651
cpw = (
    5000.0  # Mass-based heat capacity of the fluid (($\dfrac{J}{kg \cdot $^{\circ}C}$))
)
rhos = 2000.0  # Density of the solid material ($\dfrac{kg}{m^3}$)
cps = 500.0  # Mass-based heat capacity of the solid material ($\dfrac{J}{kg \cdot $^{\circ}C}$)
alpha_l = 0.0  # Longitudinal dispersivity ($m$)
ath1 = 1.0e9  # Transverse mechanical dispersivity ($m$)
ath2 = 1.0e9  # Transverse mechanical dispersivity ($m$)
T0 = 80  # Initial temperature of the active domain ($^{\circ}C$)
T1 = 30  # Injected temperature ($^{\circ}C}$)
q = 1.2649e-8  # Darcy velocity ($m/s$)


# delr is a calculated value based on the number of desired columns
assert (
    L / ncol % 1 == 0
), "reconsider specification of NCOL such that length of the simulation divided by NCOL results in an even number value"
delr = L / ncol  # Width along the row ($m$)

# Some values for the analytical solution
Q_well = q * H * thk  # Injection flow rate ($m^3/day$)
# Calculate gw reservoir heat capacity (heat capacity by volume)
rhoCp_bulk = prsity * rhow * cpw + (1 - prsity) * rhos * cps
# Calculate "thermal conduction velocity"
v = q * (rhow * cpw) / rhoCp_bulk
# Calculate bulk thermal conductivity
kt_bulk = prsity * ktw + (1 - prsity) * kts
D = kt_bulk / rhoCp_bulk + alpha_l * v

# The overburden thermal diffusion coefficient is the same as for the gw reservoir
# by virtue of how we are setting everything up (same porosity, both are saturated,
# physical properties are the same in both zones)
D_prime = D
# Heat capacity ratio
# (if different rho and Cp values are specified among the overburden and gw
#  reservoir, then this value will need to be calculated based on those values.
h_prime = 1

# Calculate grid characteristics for the DISU instantiation
top_overburden = top_res * 2
icelltype = 0
perlen = [perlen] * nper  # Convert to a list

# The flow simulation just needs 1 steady-state stress period
# Transport simulation is transient
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp, tsmult))

# GWE output control


# Solver parameters
nouter = 50
ninner = 100
hclose = 1e-7
rclose = 1e-6
# -


# +
# A few steps to get data ready for setting up the model grid
def dis_mult(thk, nlay, dis_mult):
    lay_elv = thk * ((dis_mult - 1) / (dis_mult**nlay - 1))
    return lay_elv


# For the vertical discretization, use a geometric multiplier to refine discretization at
# the reservoir-overburden interface
dismult = 1.0673002615
Delta_zbot1 = dis_mult(top_res, nlay_overburden, dismult)
botm_seq = [0.01]  # Starting thickness for the geometric expansion of layer thicknesses

for k in np.arange(1, nlay_overburden):
    botm_seq.append(botm_seq[-1] * dismult)

# Reverse the order for the overburden since botm is top -> down
botm_seq_overburden = list(reversed(botm_seq))

# Define botm for overburden
botm_overburden = [top_overburden - botm_seq_overburden[0]]
for thk_interval in botm_seq_overburden[1:]:
    next_elv = botm_overburden[-1] - thk_interval
    botm_overburden.append(next_elv)

botm_overburden = [round(itm, 2) for itm in botm_overburden]

# Continue moving down the profile and determine discretization for the sat zone
if nlay_res == 1:
    botm_res = [0]
else:
    botm_res = [top_res - botm_seq[0]]
    for thk_interval in botm_seq[1:]:
        next_elv = botm_res[-1] - thk_interval
        botm_res.append(next_elv)

    botm_res = [abs(round(itm, 2)) for itm in botm_res]


# Merge the two botm lists for coming up with a single unified botm object
# to pass to the discretization constructor
botm = botm_overburden + botm_res

# Establish a 2D array of the nodal locations
z_tops = [top_overburden] + botm
z_diffs = abs(np.diff(z_tops))
z_mids = z_tops[:-1] - (z_diffs / 2)

# Based on the way the analytical solution works, the interface altitude
# equals 0.  Moreover, the temperature at the interface equals the
# temperature in the 1D groundwater reservoir.  With that in view,
# set the last element to the interface altitude and subtract off the
# elevation of the top of the gw reservoir, which is z-axis point of
# reference
z_mids[-1] = top_res
z_mids = z_mids - 100

x_mids = [j * delr + (delr / 2) for j in np.arange(ncol)]
X, Y = np.meshgrid(x_mids, z_mids)

# -


# +
# Some functions for calculating the Barends analytical solution
def barends_eqn4(sigma, x, z, t):
    exponent = -((sigma - x * v / (4 * D * sigma)) ** 2)
    term1 = math.exp(exponent)
    term2 = x**2 * h_prime * math.sqrt(D_prime) / (8 * D * H * sigma**2)
    term3 = z / (2 * math.sqrt(D_prime))
    term4 = t - x**2 / (4 * D * sigma**2)

    # put it all together
    eqn4_val = term1 * math.erfc((term2 + term3) * (term4) ** (-0.5))

    return eqn4_val


def calc_analytical_sln(times):
    # Remember X, Y are output from MeshGrid
    myarr = np.zeros((len(times), X.shape[0], X.shape[1]))
    integral_rslt = np.zeros_like(myarr)
    prefix_rslt = np.zeros_like(myarr)
    for t_idx, t in enumerate(times):
        for k, z in enumerate(z_mids):  # rows (which are effectively z)
            for j, x in enumerate(x_mids):  # columns (which are effectively x)
                lower_lim = x / (2 * math.sqrt(D * t))
                integral = quad(barends_eqn4, lower_lim, np.inf, args=(x, z, t))
                integral_rslt[t_idx, k, j] = integral[0]

                # multiply the prefix by the solution to the integral
                prefix = (2 * (T1 - T0)) / (math.sqrt(np.pi))
                prefix_rslt[t_idx, k, j] = prefix
                result = prefix * integral[0]
                # store result for plotting
                myarr[t_idx, k, j] = result + T0

                # myarr[t_idx, k, j] = result + T0

    return myarr


# -

# +
# Function needed to establish the DISU grid with specific requirements
# for analogy with the Barends analytical solution


def flatten(xss):
    return [x for xs in xss for x in xs]


def set_connectiondata(n, lay, top, left, right, bottom):
    # Instatiate empty lists
    jas = [n]
    ihc = [lay]
    cl12 = [n]
    hwva = [n]
    angldeg = [360.0]

    # Calculate half-cell thickness up front for vertical connections
    if lay == 0:
        cl12_val = (top_overburden - botm[lay]) / 2
    else:
        cl12_val = (botm[lay - 1] - botm[lay]) / 2

    if top:
        jas.append(n - ncol)
        ihc.append(0)  # ihc = 0 for vertical connection
        cl12.append(cl12_val)  # half the cell thickness or a vertical connection
        hwva.append(delr * delc)  # for vertical connection, area is delr * delc
        angldeg.append(0.0)  # placeholder only for vertical connections

    if left:
        jas.append(n - 1)  # left
        ihc.append(delc)  # ihc = 1 for horizontal connection
        cl12.append(delr / 2)  # half the cell width along a horizontal connection
        hwva.append(
            delc
        )  # for horizontal connection, value of hwva is width along a row
        angldeg.append(
            180.0
        )  # for horizontal connection, value is 180.0 along negative x-axis

    if right:
        jas.append(n + 1)  # right
        ihc.append(delc)  # ihc = 1 for horizontal connection
        cl12.append(delr / 2)  # half the cell width along a horizontal connection
        hwva.append(
            delc
        )  # for horizontal connection, value of hwva is width along a row
        angldeg.append(
            0.0
        )  # for horizontal connection, value is 0.0 along positive x-axis

    if bottom:
        jas.append(n + ncol)  # below
        ihc.append(0)  # ihc = 0 for vertical connection
        cl12.append(cl12_val)  # half the cell thickness or a vertical connection
        hwva.append(delr * delc)  # for vertical connection, value of hwva is area
        angldeg.append(0.0)  # placeholder only for vertical connections

    return jas, ihc, cl12, hwva, angldeg


def get_conndat(n, lay, col):
    top = False
    left = False
    right = False
    bottom = False

    # iac is the number of connections (plus 1) for each cell
    iac = 1  # start with the "plus 1"

    if lay == 0:
        # For upper-most layer, only connection cell underneath it
        iac += 1
        # Bottom
        bottom = True

    elif lay <= 99:
        # For layers 2-100, connections fixed at 2 (above and below)
        iac += 2
        # Above
        top = True
        # Bottom
        bottom = True

    elif lay > 99:
        # Layer 100 is the upper most layer in the gw reservoir and where
        # horizontal connections start
        iac += 2
        if col == 0:
            # In addition to vertical connections, include 1 horizontal connecion
            # (horizontal connection will be to the right only in this case)
            iac += 1
            # Above
            top = True
            # Right
            right = True
            # Bottom
            bottom = True

        elif col == ncol - 1:
            # Horizontal connection will be to the left only
            iac += 1
            # Above
            top = True
            # Left
            left = True
            # Below
            bottom = True

        else:
            # If interior column, there will be two horizontal connections
            iac += 2
            # Above
            top = True
            # Left
            left = True
            # Right
            right = True
            # Below
            bottom = True

    jas_vals, ihc_vals, cl12_vals, hwva_vals, angldeg_vals = set_connectiondata(
        n, lay, top, left, right, bottom
    )

    # If bottom most layer, need to .pop() the last values out of the respective lists
    # This should be done because the bottom connections will always be represented by the final
    # values in the list (at this point in the development anyway)
    if lay == nlay - 1:
        iac -= 1
        jas_vals.pop(-1)
        ihc_vals.pop(-1)
        cl12_vals.pop(-1)
        hwva_vals.pop(-1)
        angldeg_vals.pop(-1)

    return iac, jas_vals, ihc_vals, cl12_vals, hwva_vals, angldeg_vals


def filter_nodes(xvc, yvc, iv, xv, yv):
    # Check if iv is empty
    if iv:
        vert_id = max(iv)

        # Add nodes represented by xvc and yvc if not already contained in xv, yv
        for new_pair in zip(xvc, yvc):
            found = False
            for old_pair in zip(xv, yv):
                if old_pair[0] == new_pair[0] and old_pair[1] == new_pair[1]:
                    found = True
                    break

            if not found:
                # if not already present, add the vertex
                iv.append(vert_id + 1)
                vert_id += 1
                xv.append(new_pair[0])
                yv.append(new_pair[1])

    else:
        iv.extend(list(range(0, 4)))
        xv.extend(xvc)
        yv.extend(yvc)

    return iv, xv, yv


def populate_linked_IDs(xvc, yvc, xv, yv):
    # Ensure 4 items contained in the passed list of nodes to get positional IDs for
    assert len(xvc) == 4, "Number of processed vertices is off"

    vertices = []
    for pair in zip(xvc, yvc):  # xvc should have 4 items
        for i, existing_pair in enumerate(zip(xv, yv)):
            if existing_pair[0] == pair[0] and existing_pair[1] == pair[1]:
                vertices.append(i)
                break

    assert len(vertices) == 4, "Number of vertices should be 4"
    return vertices


def buildout_vertex_locations():
    n = -1
    cell_vert_lkup = {}
    # Only need 1 layer's worth of vertices since they can be repeated
    iv = []
    xv = []
    yv = []

    ly = 0
    for j in np.arange(ncol):
        # There are two X locations (leverages the fact that delr = 1.0)
        vert_left = float(j * delr)
        vert_right = vert_left + delr
        # There are two Y locations (only 1 row with unit width)
        vert_front = 0.0
        vert_back = 1.0

        # First define each vertex for the cell in the current column
        if ly == 0:
            # left, back
            xvc = [vert_left]
            yvc = [vert_back]
            # right, back
            xvc.append(vert_right)
            yvc.append(vert_back)
            # right, front
            xvc.append(vert_right)
            yvc.append(vert_front)
            # left, front
            xvc.append(vert_left)
            yvc.append(vert_front)

        # Second, only keep vertices that don't already appear in the respective lists
        iv, xv, yv = filter_nodes(xvc, yvc, iv, xv, yv)

        # Store dictionary entry linking cell ID with its vertices
        vertices = populate_linked_IDs(xvc, yvc, xv, yv)
        n += 1
        cell_vert_lkup[n] = vertices

    # Now loop for the remaining layers
    for j in np.arange(ncol):
        # Only need to find the top layer's vertices once
        verts_to_use = cell_vert_lkup[j]
        for ly in np.arange(1, nlay):
            n = (ly * ncol) + j
            cell_vert_lkup[n] = verts_to_use

    return iv, xv, yv, cell_vert_lkup


def append_cell2d(n, xv_lst, yv_lst, cell_vert_lkup):
    # Get the vertex IDs for the current cell

    # Start with calculating the x location of the cell
    col_id = n % ncol
    cell_x = delr / 2 + col_id * delr
    # The y location of the cell is fixed at delc/2
    cell_y = delc / 2

    # Get associated vertices
    vertices = cell_vert_lkup[n]

    # Every cell will have 4 vertices
    new_cell2d = [[n, cell_x, cell_y, 4], vertices]

    # Flatten
    new_cell2d = flatten(new_cell2d)

    return new_cell2d


# +
# Setup model grid
iac_lst = []
ja_lst = []
ihc_lst = []
cl12_lst = []
hwva_lst = []
angldeg_lst = []

# VERTICES block input
iv_lst = []
xv_lst = []
yv_lst = []
cell_2d_lst = []
top_lst = []
bot_lst = []

# For top and bottoms, can manually build-out the lists
top_lst = [top_overburden] * ncol
for lay in np.arange(0, nlay - 1):
    top_lst.extend([botm[lay]] * ncol)

bot_lst = []
for lay in np.arange(nlay):
    bot_lst.extend([botm[lay]] * ncol)

topnp = np.array(top_lst)
botnp = np.array(bot_lst)

# Build out VERTICES block (only needs to happen once)
iv, xv_lst, yv_lst, cell_vert_lkup = buildout_vertex_locations()

vertices = []
for i in np.arange(len(iv)):
    vertices.append([iv[i], xv_lst[i], yv_lst[i]])

# Cycle through each layer and column.
for lay in np.arange(nlay):
    for col in np.arange(ncol):
        n = lay * ncol + col  # n will be zero based

        # Values for CONNECTIONDATA block
        iac, ja_cell, ihc_cell, cl12_cell, hwva_cell, angldeg_cell = get_conndat(
            n, lay, col
        )

        # accumulate connection information in lists
        iac_lst.append(iac)
        ja_lst.append(ja_cell)
        ihc_lst.append(ihc_cell)
        cl12_lst.append(cl12_cell)
        hwva_lst.append(hwva_cell)
        angldeg_lst.append(angldeg_cell)

        # Add Cell2D information
        # An example of what cell2d information should look like.
        # [997, 0.09539255, 0.022023124999999998, 4, 997, 1177, 1176, 996],
        # [998, 0.09456585, 0.025338825, 4, 998, 1178, 1177, 997],
        # [999, 0.093623925, 0.028623675, 4, 999, 1179, 1178, 998],
        cell_2d = append_cell2d(n, xv_lst, yv_lst, cell_vert_lkup)
        cell_2d_lst.append(cell_2d)


iacnp = np.array(iac_lst)
janp = np.array(flatten(ja_lst))
ihcnp = np.array(flatten(ihc_lst))
cl12np = np.array(flatten(cl12_lst))
hwvanp = np.array(flatten(hwva_lst))
angldegnp = np.array(flatten(angldeg_lst))
# -


# +
# Establish flow boundary conditions. To get a uniform flow field, need to
# thickness-weight the WEL-inserted/extracted flows since the layer thicknesses
# are not uniform
wel_spd_left = []
wel_spd_right = []
res_thickness = top_res - botm_res[-1]
for ly in np.arange(0, nlay_res):
    # calculate the thickness-weighting factor
    if ly < 1:
        thick_interval = top_res - botm_res[ly]
    else:
        thick_interval = botm_res[ly - 1] - botm_res[ly]

    Qcell = Q_well * thick_interval / res_thickness
    id_left = (ly + nlay_overburden) * ncol
    id_right = id_left + (ncol - 1)
    wel_spd_left.append(
        [id_left, Qcell, T1]
    )  # 30.0 is inflow temperature (auxiliary var)
    wel_spd_right.append([id_right, -Qcell])

wel_spd_left = {0: wel_spd_left}
wel_spd_right = {0: wel_spd_right}

# -


# +
@timed
def build_mf6_flow_model():
    gwf_name = sim_name
    sim_ws = os.path.join(workspace, sim_name, "mf6gwf")

    # Instantiate a MODFLOW 6 simulation
    sim = flopy.mf6.MFSimulation(sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6")

    # Instantiate time discretization package
    flopy.mf6.ModflowTdis(
        sim,
        nper=1,  # just one steady state stress period for gwf
        perioddata=[(1.0, 1, 1.0)],
        time_units=time_units,
    )

    # Instantiate Iterative model solution package
    flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="bicgstab",
        outer_maximum=nouter,
        outer_dvclose=hclose,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=f"{rclose} strict",
    )

    # Instantiate a groundwater flow model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)

    # Instantiate an unstructured discretization package
    flopy.mf6.ModflowGwfdisu(
        gwf,
        length_units=length_units,
        nogrb=True,
        nodes=len(topnp),
        nja=len(janp),
        nvert=len(iv),
        top=topnp,
        bot=botnp,
        area=delr * delc,
        idomain=1,
        iac=iacnp,
        ja=janp,
        ihc=ihcnp,
        cl12=cl12np,
        hwva=hwvanp,
        angldegx=angldegnp,
        vertices=vertices,
        cell2d=cell_2d_lst,
        pname="DISU-Barends",
        filename="{}.disu".format(gwf_name),
    )

    # Instantiate node-property flow (NPF) package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_saturation=True,
        save_specific_discharge=True,
        xt3doptions=False,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        pname="NPF-1",
        filename="{}.npf".format(gwf_name),
    )

    # Instantiate initial conditions package for the GWF model
    flopy.mf6.ModflowGwfic(
        gwf, strt=top_overburden, pname="IC-1", filename="{}.ic".format(gwf_name)
    )

    # Instantiating MODFLOW 6 storage package
    # (steady flow conditions, so no actual storage,
    # using to print values in .lst file)
    flopy.mf6.ModflowGwfsto(
        gwf,
        ss=0,
        sy=0,
        steady_state={0: True},
        pname="STO",
        filename="{}.sto".format(gwf_name),
    )

    # Instantiate WEL package for adding water on the left and removing
    # it on the right
    flopy.mf6.ModflowGwfwel(
        gwf,
        auxiliary="TEMPERATURE",
        stress_period_data=wel_spd_left,
        pname="WEL-left",
        filename="{}.wel-left".format(gwf_name),
    )
    # Instantiate WEL package for extracting water on the right
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wel_spd_right,
        pname="WEL-right",
        filename="{}.wel-right".format(gwf_name),
    )

    # Instantiating MODFLOW 6 output control package (flow model)
    head_filerecord = "{}.hds".format(sim_name)
    budget_filerecord = "{}.cbc".format(sim_name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim


# -


# +
def build_mf6_heat_model():
    print("Building mf6gwe model...{}".format(sim_name))
    gwename = sim_name
    sim_ws = os.path.join(workspace, sim_name, "mf6gwe")

    sim = flopy.mf6.MFSimulation(sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6")

    # Instantiating MODFLOW 6 groundwater transport model
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file="{}.nam".format(gwename),
    )

    # Instantiate Iterative model solution package
    imsgwe = flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="bicgstab",
        complexity="SIMPLE",
        outer_maximum=nouter,
        outer_dvclose=hclose,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=f"{rclose} strict",
    )
    sim.register_ims_package(imsgwe, [gwe.name])

    # MF6 time discretization differs from corresponding flow simulation
    flopy.mf6.ModflowTdis(
        sim, nper=len(tdis_rc), perioddata=tdis_rc, time_units=time_units
    )

    # Instantiate an unstructured discretization package
    flopy.mf6.ModflowGwedisu(
        gwe,
        length_units=length_units,
        nogrb=True,
        nodes=len(topnp),
        nja=len(janp),
        nvert=len(iv),
        top=topnp,
        bot=botnp,
        area=delr * delc,
        idomain=1,
        iac=iacnp,
        ja=janp,
        ihc=ihcnp,
        cl12=cl12np,
        hwva=hwvanp,
        angldegx=angldegnp,
        vertices=vertices,
        cell2d=cell_2d_lst,
        pname="DISU-Barends",
        filename="{}.disu".format(gwename),
    )

    # Instantiating MODFLOW 6 heat transport initial temperature
    flopy.mf6.ModflowGweic(
        gwe, strt=T0, pname="IC-gwe", filename="{}.ic".format(gwename)
    )

    # Instantiating MODFLOW 6 heat transport advection package
    flopy.mf6.ModflowGweadv(
        gwe, scheme=scheme, pname="ADV", filename="{}.adv".format(gwename)
    )

    # Instantiating MODFLOW 6 heat transport dispersion package
    if ktw != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            alh=alpha_l,
            ath1=ath1,
            ath2=ath2,
            ktw=ktw,
            kts=kts,
            pname="CND",
            filename="{}.cnd".format(gwename),
        )

    # Instantiating MODFLOW 6 heat transport mass storage package (consider renaming to est)
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=prsity,
        heat_capacity_solid=cps,
        heat_capacity_water=cpw,
        density_solid=rhos,
        density_water=rhow,
        pname="EST",
        filename="{}.est".format(gwename),
    )

    flopy.mf6.ModflowGwectp(
        gwe,
        save_flows=True,
        maxbound=1,
        stress_period_data={0: [(ncol * nlay_overburden), T1]},
        pname="CTP",
        filename="{}.ctp".format(gwename),
    )

    # Instantiating MODFLOW 6 source/sink mixing package for dealing with
    # auxiliary temperature specified in WEL boundary package.
    sourcerecarray = [("WEL-left", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwessm(
        gwe, sources=sourcerecarray, pname="SSM", filename="{}.ssm".format(gwename)
    )

    # Instantiating MODFLOW 6 heat transport output control package
    # e.g., day 100 = 100 * 86400 = 8,640,000 = 8.64e6
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord="{}.cbc".format(gwename),
        temperature_filerecord="{}.ucn".format(gwename),
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord={
            99: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
            100: [],  # Toggle output off
            199: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
            200: [],  # Toggle output off
            299: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        },
        # 300: [],  # Toggle output off
        # 399: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        # 400: [],  # Toggle output off
        # 499: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        # 500: [],  # Toggle output off
        # 599: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        # 600: [],  # Toggle output off
        # 699: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        # 700: [],  # Toggle output off
        # 799: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        # 800: [],  # Toggle output off
        # 899: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        # 900: [],  # Toggle output off
        # 999: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],},
        printrecord={99: [("BUDGET", "LAST")]},  # ,
        #             999: [("BUDGET", "LAST")]}
    )

    # Instantiating MODFLOW 6 Flow-Model Interface package
    pd = [
        ("GWFHEAD", "../mf6gwf/" + sim_name + ".hds", None),
        ("GWFBUDGET", "../mf6gwf/" + sim_name + ".cbc", None),
    ]
    flopy.mf6.ModflowGwefmi(gwe, packagedata=pd)

    return sim


# -


# +
def write_mf6_models(sim_gwf, sim_gwe, silent=True):
    # Run the steady-state flow model
    if sim_gwf is not None:
        sim_gwf.write_simulation(silent=silent)

    # Second, run the heat transport model
    if sim_gwe is not None:
        sim_gwe.write_simulation(silent=silent)


def run_model(sim, silent=True):
    # Attempting to run model
    success, buff = sim.run_simulation(silent=silent, report=True)
    if not success:
        print(buff)
    return success


# -


# +
def gen_faux_grid():
    # Mimicking DISU grid with a DIS grid. DISU plotting is much
    # too slow for whatever reason.
    delc_local = delc * np.ones(nrow)
    delr_local = delr * np.ones(ncol)
    top_overburden_local = top_overburden * np.ones((nrow, ncol))
    botm_local = []
    for bot in botm:
        bot_arr = bot * np.ones((nrow, ncol))
        botm_local.append(bot_arr)

    botm_local = np.array(botm_local)

    sgr = fgrid.StructuredGrid(
        delc_local,
        delr_local,
        top=top_overburden_local,
        botm=botm_local,
        xoff=0.0,
        yoff=0.0,
        angrot=0.0,
    )

    return sgr


# -


# +
@timed
def plot_thermal_bleeding(sim_gwe):
    # Select times to compare (year 100 is 100 * 365 * 86400
    #                          year 200 is 200 * 365 * 86400, etc.)
    # The following is the full slate of times which could be compared
    # times = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    times = [300]
    times = [86400 * 365 * tm for tm in times]
    analytical_answers = calc_analytical_sln(times)

    gwename = sim_name
    gwe = sim_gwe.get_model(gwename)

    # Retrieve simulated temperature field
    temperature = gwe.output.temperature().get_alldata()

    # Get a DIS representation of the DISU grid
    sgr = gen_faux_grid()

    # Reshape results to fit with DIS results
    dis_temps = []
    temp_arr = temperature[-1].squeeze()  # Corresponds to first requested output time
    temp_arr = temp_arr.reshape((nlay, nrow, ncol))
    dis_temps.append(temp_arr)
    dis_temps = np.array(dis_temps)

    figsize = (6, 4)
    fig, ax = plt.subplots(figsize=figsize)
    mx = flopy.plot.PlotCrossSection(
        ax=ax, modelgrid=sgr, line={"Row": 0}, extent=(0, 100, 70, 125)
    )
    mx.plot_grid()
    # cb_alt = mx.contour_array(analytical_answers[-1], levels=[30, 40, 50, 60, 70, 80], linewidth=0.1)

    # save figure
    if plot_show:
        plt.show()
    if plot_save:
        fpth = os.path.join(figs_path / "{}-gridView{}".format(sim_name, ".png"))
        fig.savefig(fpth, dpi=600)

    # Next, plot model output compared to analytical solution
    fontsize = 8
    fig2, ax2 = plt.subplots(figsize=figsize)
    mx2 = flopy.plot.PlotCrossSection(
        ax=ax2, modelgrid=sgr, line={"Row": 0}, extent=(0, 500, 50, 175)
    )
    levels = [32, 40, 50, 60, 70, 75]
    manual_locations = [
        (25, 101),
        (50, 110),
        (85, 118),
        (120, 130),
        (150, 140),
        (175, 145),
    ]

    cb_alt = mx2.contour_array(
        analytical_answers[-1], levels=levels, linewidths=1.0, colors="k"
    )
    cb_sim = mx2.contour_array(
        dis_temps[-1, :, 0, :],
        levels=levels,
        linewidths=1.4,
        colors="r",
        linestyles="dashed",
    )
    ax2.tick_params(axis="both", which="major", labelsize=fontsize)
    labels = ax2.clabel(
        cb_sim,
        cb_sim.levels,
        inline=False,
        inline_spacing=0.0,
        fontsize=fontsize,
        fmt="%1d",
        colors="k",
        manual=manual_locations,
    )

    for label in labels:
        label.set_bbox(dict(facecolor="white", pad=2, ec="none"))

    ax2.text(350, 103, "Overburden", fontsize=fontsize)
    ax2.text(350, 87, "Groundwater\nReservoir", fontsize=fontsize)
    ax2.set_xlabel("X coordinate, m", fontsize=fontsize)
    ax2.set_ylabel("Elevation, m", fontsize=fontsize)
    custom_lines = [
        Line2D([0], [0], color="k", lw=1),
        Line2D([0], [0], color="r", lw=1.4, linestyle="dashed"),
    ]

    ax2.legend(
        custom_lines,
        ["Barends Solution", "MODFLOW 6 GWE"],
        loc="upper right",
        frameon=False,
        fontsize=fontsize,
    )

    cb_sim.set_dashes([(0, (5.0, 6.0))])

    ax2.axhline(y=100, color="grey", linestyle="-")
    # save figure
    if plot_show:
        plt.show()
    if plot_save:
        fpth = os.path.join(figs_path / "{}-300yrs{}".format(sim_name, ".png"))
        fig2.savefig(fpth, dpi=600)

    return


# -


# +
def scenario(idx, silent=False):
    sim_gwf = build_mf6_flow_model()
    sim_gwe = build_mf6_heat_model()

    if sim_gwf is not None and sim_gwe is not None:
        write_mf6_models(sim_gwf, sim_gwe, silent=silent)

    success = run_model(sim_gwf, silent=silent)
    if success:
        success = run_model(sim_gwe, silent=silent)

    if success:
        plot_thermal_bleeding(sim_gwe)


# -

# +
scenario(0, silent=False)
# -
