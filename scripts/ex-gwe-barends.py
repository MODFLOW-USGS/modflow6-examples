#
# An analytical solution for this example was originally published in
#
#  Barends, F.B.J., 2010. Complete Solution for Transient Heat Transport in
#    Porous Media, Following Lauwerier's Concept. Society of Petroleum '
#    Engineers, Annual Technical Conference and Exhibition, Florence, Italy,
#    19–22 September 2010.
#    https://doi.org/10.2118/134670-MS
#
# Below is a diagram of the model cell numbering, which is helpful for this
# DISU setup since the groundwater reservoir is flowing in a 1D left-to-right
# manner with no flow in the overburden.  Additionally, heat "bleeds" into
# the overburden from the flowing saturated zone.  Heat bleeds in a 1D upward
# manner and does not conduct laterally once in the overburden.  Thus, the
# connections in the DISU grid are such that cells in the overburden are not
# connected to their lateral neighbors, but only the cells above and below.
#
#           +-------+-------+-------+     +-------+-------+-------+
# Cell IDs  |   0   |   1   |   2   | ... |  997  |  998  |  999  | Layer 1
# (0-based) +-------+-------+-------+     +-------+-------+-------+
#           | 1,000 | 1,001 | 1,002 | ... | 1,997 | 1,998 | 1,999 | Layer 2
#           +-------+-------+-------+     +-------+-------+-------+
#           | 2,000 | 2,001 | 2,002 | ... | 2,997 | 2,998 | 2,999 | Layer 3
#           +-------+-------+-------+     +-------+-------+-------+
#               .       .       .             .       .       .
#               .       .       .             .       .       .
#               .       .       .             .       .       .
#           +-------+-------+-------+     +-------+-------+-------+        
#           |197,000|197,001|197,002| ... |197,997|197,998|197,999| Layer 1
#           +-------+-------+-------+     +-------+-------+-------+        
#           |198,000|198,001|198,002| ... |198,997|198,998|198,999| Layer 2
#           +-------+-------+-------+     +-------+-------+-------+        
#           |199,000|199,001|199,002| ... |199,997|199,998|199,999| Layer 3
#           +-------+-------+-------+     +-------+-------+-------+        
#                       ----> gw flow direction ---->
# +
import os
import pathlib as pl
from pprint import pformat

import flopy
import flopy.discretization as fgrid
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import quad
from flopy.plot.styles import styles
from modflow_devtools.misc import get_env, timed

import datetime

# Example name and base workspace
sim_name = "barends"
workspace = pl.Path("../examples")

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
time_units = "days"

nrow = 1
ncol = 1000  # Number of columns ($-$)
delr = 1.0  # Width along the row ($m$)
delc = 1.0  # Width along the column ($m$)
nlay_satzn = 1  # Number of layers in the saturated zone ($-$)
nlay_overburden = 100  # Number of layer representing the overburden ($-$)
nlay = nlay_overburden + nlay_satzn  # Total number of layers ($-$)

top_satzn = 100.0  # Elevation of saturated zone ($m$)
H = 100.0  # Reservoir height ($m$)
thk = 1.0  # Unit thickness ("into the page" $m$)

k11 = 1.0  # Hydraulic conductivity ($m/d$)
k33 = 1e-5  # Vertical hydraulic conductivity ($m/d$)
Q = 10  # Injection flow rate ($m^3/day$)
perlen = 1.0  # Period length ($days$)
nstp = 1  # Number of time steps per stress period ($-$)
tsmult = 1.0  # Time step multiplier ($-$)
nper = 1000  # Number of simulated stress periods
prsity = 0.2  # Porosity ($-$)
scheme = "TVD"  # Advection solution scheme ($-$)
ktw = 0.6  # Thermal conductivity of the fluid ($\dfrac{W}{m \cdot ^{\circ}C}$)
kts = 3.0  # Thermal conductivity of the aquifer matrix (which is also the dry overburden material) ($\dfrac{W}{m \cdot ^{\circ}C}$)
kth = 0.25  # Longitudinal thermal conductivity in the saturated zone ($\frac{W}{m \cdot ^{\circ}C}$)
kttransverse = 2.50  # Transverse thermal conductivity in the saturated zone ($\frac{W}{m \cdot ^{\circ}C}$)
rhow = 1000.0  # Density of water ($\frac{kg}{m^3}$) # 3282.296651
cpw = 4180.0  # Mass-based heat capcity of the fluid (($\dfrac{J}{kg \cdot $^{\circ}C}$))
rhos = 2650.0  # Density of the aquifer material ($\dfrac{kg}{m^3}$)
cps = 800.0  # Mass-based heat capacity of the overburden material ($\dfrac{J}{kg \cdot $^{\circ}C}$)
al = 0.0  # Longitudinal mechanical dispersivity ($\frac{m^2}{day}$)
alpha_l = 0.2  # Longitudinal dispersivity ($m$)
ath2 = 0.0  # Transverse mechanical dispersivity ($\frac{m^2}{day}$)
lhv = 2500.0  # Latent heat of vaporization ($$)
T0 = 80  # Initial temperature of the active domain ($^{\circ}C$)
T1 = 30  # Injected temperature ($^{\circ}C}$)

# Some values for the analytical solution
days2seconds = 86400
q = Q / (prsity * H * thk)  # Darcy velocity
q *= 1/86400  # Convert from units of days to seconds
# Calculate "soil" heat capacity (it seems Barends refers to the sat zone as the soil zone
rhoCp_bulk = prsity * rhow * cpw + (1-prsity) * rhos * cps
# Calculate thermal conduction velocity
v = q * (rhow * cpw) / rhoCp_bulk
# Calculate thermal longitudinal dispersion-diffusion coefficient
kt_bulk = prsity * ktw + (1 - prsity) * kts
# In Barends (2010), D = K_t_b/(rho_b * C_p_b) + alpha_l * v
# For the present purpose, D is calculated without the alpha_l * v term
# since, for the GWE model, alpha_l is going to be used to get the model
# to do the equivalent of the analytical solution
D = kt_bulk/rhoCp_bulk  # + alpha_l * v
# Calculate an alpha_l to be used by the model to get an equivalent
# thermal diffusive spread in the model
alpha_l = D / q
# Calculate overburden thermal diffusion coefficient
D_prime = kts / (rhos * cps)
# Calculate an alpha_th1 (and alpha_th2) value that is essentially equal
# to D_prime for representing anisotropic thermal diffusion (conduction)
# with the mechanical dispersion terms
alpha_th = D_prime / q
# Heat capacity ratio
h_prime = (rhos * cps) / rhoCp_bulk

# For the thermal conductivity of the solid material, need to set up lists
# equal to the number of cells such that Kts is what it needs to be in the
# overburden, but 0.0 in the aquifer
alh = [0.0 if i < (nlay_overburden * ncol) else alpha_l for i in np.arange(nlay * nrow * ncol)]
ath1 = [0.0 if i < (nlay_overburden * ncol) else alpha_th for i in np.arange(nlay * nrow * ncol)]
ath2 = ath1.copy()
Kts = [kts if i < (nlay_overburden * ncol) else 0.0 for i in np.arange(nlay * nrow * ncol)]
Cps = [cps if i < (nlay_overburden * ncol) else 0.0 for i in np.arange(nlay * nrow * ncol)]
Rhos = [rhos if i < (nlay_overburden * ncol) else 0.0 for i in np.arange(nlay * nrow * ncol)]


# Calculate grid characteristics for the DISU instantiation
top_overburden = top_satzn * 2
icelltype = [1 for i in np.arange((nlay-nlay_satzn)*ncol)] + [0 for i in np.arange(nlay_satzn * ncol)]
wetdry = [0 for i in np.arange((nlay-nlay_satzn)*ncol)] + [0.01 for i in np.arange(nlay_satzn * ncol)]
perlen = [perlen] * nper  # Convert to a list

# The flow simulation just needs 1 steady-state stress period

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp, tsmult))


# Solver parameters
nouter = 50
ninner = 100
hclose = 1e-7
rclose = 1e-6
# -

# 
def dis_mult(thk, nlay, dis_mult):
    lay_elv = thk * ((dis_mult - 1) / (dis_mult ** nlay - 1))
    return lay_elv


# For the vertical discretization, use a geometric multiplier to refine discretization at
# the reservoir-overburden interface
dismult = 1.0673002615
Delta_zbot1 = dis_mult(top_satzn, nlay_overburden, dismult)
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
if nlay_satzn == 1:
    botm_satzn = [0]
else:
    botm_satzn = [top_satzn - botm_seq[0]]
    for thk_interval in botm_seq[1:]:
        next_elv = botm_satzn[-1] - thk_interval
        botm_satzn.append(next_elv)

    botm_satzn = [abs(round(itm, 2)) for itm in botm_satzn]


# Merge the two botm lists for coming up with a single unified botm object
# to pass to the discretization constructor
botm = botm_overburden + botm_satzn

# Establish a 2D array of the nodal locations
z_tops = [top_overburden] + botm
z_diffs = abs(np.diff(z_tops))
z_mids = z_tops[:-1] - (z_diffs / 2)

# Based on the way the analytical solution works, the interface altitude
# equals 0.  Moreover, the temperature at the interface equals the
# temperature in the 1D reservoir.  With that in view, set the last element
# to the interface altitude and subtract off the elevation of the top of
# the saturated zone, which is y-axis point of reference
z_mids[-1] = top_satzn
z_mids = z_mids - 100

x_mids = [j * delr + 0.5 for j in np.arange(ncol)]
X, Y = np.meshgrid(x_mids, z_mids)


# +
# Some functions for calculating the Barends analytical solution
def barends_eqn4(sigma, x, z, t):
    exponent = -(sigma - x * v / ( 4 * D * sigma ) ) ** 2
    term1 = math.exp(exponent)
    term2 = x ** 2 * h_prime * math.sqrt(D_prime) / (8 * D * H * sigma ** 2)
    term3 = z / (2 * math.sqrt(D_prime))
    term4 = t - x ** 2 / (4 * D * sigma ** 2)

    # put it all together
    eqn4_val = term1 * math.erfc((term2 + term3) * (term4) ** (-0.5))

    return eqn4_val


def calc_analytical_sln(times):

    times = [tm * days2seconds for tm in times]
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

                #myarr[t_idx, k, j] = result + T0

    return myarr

# -

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
        hwva.append(delr*delc)  # for vertical connection, area is 1.0m x 1.0m
        angldeg.append(0.0)  # placeholder only for vertical connections
    
    if left:
        jas.append(n - 1)  # left
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delc / 2)  # half the cell width along a horizontal connection
        hwva.append(delr)  # for horizontal connection, value of hwva is width along a row
        angldeg.append(180.0)  # for horizontal connection, value is 180.0 along negative x-axis
    
    if right:
        jas.append(n + 1)  # right
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delc / 2)  # half the cell width along a horizontal connection
        hwva.append(delc)  # for horizontal connection, value of hwva is width along a row
        angldeg.append(0.0)  # for horizontal connection, value is 0.0 along positive x-axis
    
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
        # Layer 100 is the upper most layer in the saturated zone and where
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

    jas_vals, ihc_vals, cl12_vals, hwva_vals, angldeg_vals = set_connectiondata(n, lay, top, left, right, bottom)
    
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
        vert_left = float(j)
        vert_right = float(j) + delr
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
for lay in np.arange(0, nlay-1):
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
        iac, ja_cell, ihc_cell, cl12_cell, hwva_cell, angldeg_cell = get_conndat(n, lay, col)

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

# Establish flow boundary conditions. To get a uniform flow field, need to
# thickness-weight the WEL-inserted/extracted flows since the layer thicknesses
# are not uniform
wel_spd_left = []
wel_spd_right = []
satzn_thickness = (top_satzn - botm_satzn[-1])
for ly in np.arange(0, nlay_satzn):
    # calculate the thickness-weighting factor
    if ly < 1:
        thick_interval = top_satzn - botm_satzn[ly]
    else:
        thick_interval = botm_satzn[ly - 1] - botm_satzn[ly]

    Qcell = Q * thick_interval / satzn_thickness
    id_left = (ly + nlay_overburden) * ncol
    id_right = id_left + (ncol - 1)
    wel_spd_left.append([id_left, Qcell, T1])  # 30.0 is inflow temperature (auxiliary var)
    wel_spd_right.append([id_right, -Qcell])

wel_spd_left = {0: wel_spd_left}
wel_spd_right = {0: wel_spd_right}

# +
def build_mf6_flow_model():

    gwf_name = 'gwf-' + sim_name
    sim_ws = os.path.join(workspace, "ex-gwe-" + sim_name, "mf6gwf")

    # Instantiate a MODFLOW 6 simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name,
        sim_ws=sim_ws,
        exe_name="mf6"
    )

    # Instantiate time discretization package
    flopy.mf6.ModflowTdis(
        sim,
        nper=1,  # just one steady state stress period for gwf
        perioddata=[(1.0, 1, 1.0)],
        time_units=time_units
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
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name,
        save_flows=True
    )

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
        area=delr*delc,
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
        filename="{}.disu".format(gwf_name)
    )

    # Instantiate node-property flow (NPF) package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_saturation=True,
        save_specific_discharge=True,
        xt3doptions=False,
        rewet_record=[("WETFCT", 0.01, "IWETIT", 1, "IHDWET", 0)],
        icelltype=icelltype,
        k=k11,
        k33=k33,
        wetdry=wetdry,
        pname="NPF-1",
        filename="{}.npf".format(gwf_name)
    )

    # Instantiate initial conditions package for the GWF model
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=top_satzn,
        pname="IC-1",
        filename="{}.ic".format(gwf_name)
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
        filename="{}.sto".format(gwf_name)
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
        filename="{}.wel-right".format(gwf_name)
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


def build_mf6_heat_model():

    print("Building mf6gwe model...{}".format(sim_name))
    gwename = 'gwe-' + sim_name
    sim_ws = os.path.join(workspace, "ex-gwe-" + sim_name, "mf6gwe")

    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name, sim_ws=sim_ws, exe_name="mf6"
    )

    # Instantiating MODFLOW 6 groundwater transport model
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file="{}.nam".format(gwename)
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
        sim,
        nper=len(tdis_rc),
        perioddata=tdis_rc,
        time_units=time_units
    )

    # Instantiate an unstructured discretization package
    flopy.mf6.ModflowGwfdisu(
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
        filename="{}.disu".format(gwename)
    )

    # Instantiating MODFLOW 6 heat transport initial temperature
    flopy.mf6.ModflowGweic(
        gwe,
        strt=T0,
        pname="IC-gwe",
        filename="{}.ic".format(gwename)
    )

    # Instantiating MODFLOW 6 heat transport advection package
    flopy.mf6.ModflowGweadv(
        gwe,
        scheme=scheme,
        pname="ADV-1",
        filename="{}.adv".format(gwename)
    )

    # Instantiating MODFLOW 6 heat transport dispersion package
    if ktw != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            alh=alh,
            ath1=ath1,
            ath2=ath2,
            ktw=0.0,
            kts=Kts,
            pname="CND-1",
            filename="{}.cnd".format(gwename),
        )

    # Instantiating MODFLOW 6 heat transport mass storage package (consider renaming to est)
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=prsity,
        cps=Cps,
        heat_capacity_water=cpw,
        rhos=Rhos,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        pname="EST-1",
        filename="{}.est".format(gwename),
    )

    # Instantiating MODFLOW 6 source/sink mixing package for dealing with
    # auxiliary temperature specified in WEL boundary package.
    sourcerecarray = [("WEL-left", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwessm(
        gwe,
        sources=sourcerecarray,
        pname="SSM-1",
        filename="{}.ssm".format(gwename)
    )

    # Instantiating MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord="{}.cbc".format(gwename),
        temperature_filerecord="{}.ucn".format(gwename),
        temperatureprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord={99: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    100: [],  # Toggle output off
                    199: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    200: [],  # Toggle output off
                    299: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    300: [],  # Toggle output off
                    399: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    400: [],  # Toggle output off
                    499: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    500: [],  # Toggle output off
                    599: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    600: [],  # Toggle output off
                    699: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    700: [],  # Toggle output off
                    799: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    800: [],  # Toggle output off
                    899: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
                    900: [],  # Toggle output off
                    999: [("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],},
                    printrecord={99: [("BUDGET", "LAST")],
                                 999: [("BUDGET", "LAST")]}
    )

    # Instantiating MODFLOW 6 Flow-Model Interface package
    pd = [
        ("GWFHEAD", "../mf6gwf/" + sim_name + ".hds", None),
        ("GWFBUDGET", "../mf6gwf/" + sim_name + ".cbc", None),
    ]
    flopy.mf6.ModflowGwefmi(gwe, packagedata=pd)

    return sim

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
        angrot=0.0
    )

    return sgr

def plot_thermal_bleeding(sim_gwe):

    time_indexs = {100: 0, 200: 1, 300: 2, 400: 3, 500: 4, 600: 5, 700: 6, 800: 7, 900: 8, 1000: 9}
    # Select times to compare
    times = [200, 400]
    #now1 = datetime.datetime.now()
    analytical_answers = calc_analytical_sln(times)
    #now2 = datetime.datetime.now()
    #print("Time taken: " + str(now))

    gwename = 'gwe-' + sim_name
    gwe = sim_gwe.get_model(gwename)

    # Retrieve simulated temperature field
    temperature = gwe.output.temperature().get_alldata()

    # Get a DIS representation of the DISU grid
    sgr = gen_faux_grid()

    # Reshape results to fit with DIS results
    dis_temps = []
    for tm in times:
        tm_idx = time_indexs[tm]
        temp_arr = temperature[tm_idx].squeeze()
        temp_arr = temp_arr.reshape((nlay, nrow, ncol))
        dis_temps.append(temp_arr)

    dis_temps = np.array(dis_temps)


    # modelgrid = gwe.modelgrid
    #modelgrid._ncpl = np.array([modelgrid.nnodes, ])
    #csv = modelgrid.cross_section_vertices

    #nnodes = modelgrid.nnodes
    #ncpl = modelgrid.ncpl
    #assert nnodes == ncpl, "nnodes and ncpl should be equal, but are not - rats!"

    figsize = (9, 6)
    fig, ax = plt.subplots(figsize=figsize)
    mx = flopy.plot.PlotCrossSection(ax=ax, modelgrid=sgr, line={"Row": 0})
    mx.plot_grid()
    #cb = mx.plot_array(dis_temps[-1], cmap="jet", vmin=30, vmax=80)
    cb_alt = mx.contour_array(analytical_answers[-1], levels=[30, 40, 50, 60, 70, 80])
    #cb = mx.plot_array(analytical_answers[-1], cmap="jet", vmin=30, vmax=80)
    plt.show()

    # Example script to run
    #pd.Series(analytical_answers[1, :-1, 4]).plot(logy=True)
    #pd.Series(dis_temps[1, :-1, 0, 4]).plot(logy=True)
    return

def scenario(idx, silent=False):
    sim_gwf = build_mf6_flow_model()
    sim_gwe = build_mf6_heat_model()

    if sim_gwf is not None and sim_gwe is not None:
        write_mf6_models(sim_gwf, sim_gwe)

    success = run_model(sim_gwf)
    if success:
        success = run_model(sim_gwe)

    if success:
        plot_thermal_bleeding(sim_gwe)


# +
scenario(0)
# -