============================
MODFLOW 6 – Example problems
============================

TWRI Example
============

This example is a modified version of the original MODFLOW example
(TWRI) described in (McDonald and Harbaugh 1988) and duplicated in
(Harbaugh and McDonald 1996). This problem is also is distributed with
MODFLOW-2005(Harbaugh 2005) and MODFLOW 6(Langevin et al. 2017). The
problem has been modified from a quasi-3D problem, where confining beds
are not explicitly simulated, to an equivalent three-dimensional
problem.

Example description
-------------------

There are three simulated aquifers, which are separated from each other
by confining layers (fig. `1 <#fig:ex-gwf-twri-01>`__). The confining
beds are 50 :math:`ft` thick and are explicitly simulated as model
layers 2 and 4, respectively. Each layer is a square 75,000 :math:`ft`
on a side and is divided into a grid with 15 rows and 15 column, which
forms squares 5,000 :math:`ft` on a side. A single steady-stress period
with a total length of 86,400 seconds (1 day) is simulated.

.. figure:: _images/twri-system.png
   :name: fig:ex-gwf-twri-01

   Illustration of the system simulated in the TWRI example problem
   (from (McDonald and Harbaugh 1988)).

The transmissivity of the middle and lower aquifers
(fig. `1 <#fig:ex-gwf-twri-01>`__) was converted to a horizontal
hydraulic conductivity using the layer thickness
(table `1 <#tab:ex-gwf-twri-01>`__). The vertical hydraulic conductivity
in the aquifers was set equal to the horizontal hydraulic conductivity.
The vertical hydraulic conductivity of the confining units was
calculated from the vertical conductance of the confining beds defined
in the original problem and the confining unit thickness
(table `1 <#tab:ex-gwf-twri-01>`__); the horizontal hydraulic
conductivity of the confining bed was set to the vertical hydraulic
conductivity and results in vertical flow in the confining unit.

.. container::
   :name: tab:ex-gwf-twri-01

   .. table:: Model Parameters

      +----------------------------------+----------------------------------+
      | **Parameter**                    | **Value**                        |
      +==================================+==================================+
      | Number of periods                | 1                                |
      +----------------------------------+----------------------------------+
      | Number of layers                 | 5                                |
      +----------------------------------+----------------------------------+
      | Number of columns                | 15                               |
      +----------------------------------+----------------------------------+
      | Number of rows                   | 15                               |
      +----------------------------------+----------------------------------+
      | Column width (:math:`ft`)        | 5000.0                           |
      +----------------------------------+----------------------------------+
      | Row width (:math:`ft`)           | 5000.0                           |
      +----------------------------------+----------------------------------+
      | Top of the model (:math:`ft`)    | 200.0                            |
      +----------------------------------+----------------------------------+
      | Layer bottom elevations          | -150.0, -200.0, -300.0, -350.0,  |
      | (:math:`ft`)                     | -450.0                           |
      +----------------------------------+----------------------------------+
      | Starting head (:math:`ft`)       | 0.0                              |
      +----------------------------------+----------------------------------+
      | Cell conversion type             | 1, 0, 0, 0, 0                    |
      +----------------------------------+----------------------------------+
      | Horizontal hydraulic             | 1.0e-3, 1.0e-8, 1.0e-4, 5.0e-7,  |
      | conductivity (:math:`ft/s`)      | 2.0e-4                           |
      +----------------------------------+----------------------------------+
      | Vertical hydraulic conductivity  | 1.0e-3, 1.0e-8, 1.0e-4, 5.0e-7,  |
      | (:math:`ft/s`)                   | 2.0e-4                           |
      +----------------------------------+----------------------------------+
      | Recharge rate (:math:`ft/s`)     | 3e-8                             |
      +----------------------------------+----------------------------------+

An initial head of zero :math:`ft` was specified in all model layers.
Any initial head exceeding the bottom of model layer 1 (-150 :math:`ft`)
could be specified since the model is steady-state.

Flow into the system is from infiltration from precipitation and was
represented using the recharge (RCH) package. A constant recharge rate
of :math:`3 \times 10^{-7}` :math:`ft/s` was specified for every cell in
model layer 1. Flow out of the model is from buried drain tubes
represented by drain (DRN) package cells in model layer 1, discharging
wells represented by well (WEL) package cells in all three aquifers, and
a lake represented by constant head (CHD) packages cells in the
unconfined and middle aquifers (fig. `1 <#fig:ex-gwf-twri-01>`__).

Example Results
---------------

Simulated results in the unconfined, middle, and lower aquifers are
shown in figure `2 <#fig:ex-gwf-twri-02>`__. Simulated results for a
quasi-3D MODFLOW-2005 simulation are also shown in
figure `2 <#fig:ex-gwf-twri-02>`__. MODFLOW 6 and MODFLOW-2005 results
differ by less that 0.05 :math:`ft` in any aquifer unit.

.. figure:: _images/ex-gwf-twri01.png
   :name: fig:ex-gwf-twri-02

   Simulated water levels and normalized specific discharge vectors in
   the unconfined, upper, and lower aquifers. *A*. MODFLOW 6 unconfined
   aquifer results. *B*. MODFLOW 6 middle aquifer results. *C*.
   MODFLOW 6 lower aquifer results. *D*. MODFLOW-2005 unconfined aquifer
   results. *E*. MODFLOW-2005 middle aquifer results. *F*. MODFLOW-2005
   lower aquifer results.

.. _multi-model-simulation-with-local-grid-refinement-example-3-in:

Multi-model Simulation with Local Grid Refinement (Example 3 in (Mehl and Hill 2013)
====================================================================================

In this example, two models are arranged in a nested parent-child grid
setup, where the parent model uses a coarse grid that surrounds the
finely-gridded child model (Figure `3 <#fig:mvr_lgr_planview>`__).
Parent-child model setups limit grid refinement to localized regions of
interest within large model domains to keep model runtimes to a minimum.
Owing to the relative ease of setup for this type of model arrangement
in MODFLOW 6, particularly when developing a model (i.e., simulation)
with a support utility like flopy (Bakker et al. 2016), the use of
multi-model arrangements necessitates a generalized supporting package
like MVR to transfer water among features within a simulation. To
demonstrate, this example uses MVR to cascade streamflow from an
upstream SFR reach in the parent model to a downstream SFR reach in the
child model. Further downstream MVR again cascades streamflow from the
last child model SFR reach to the appropriate SFR reach in the parent
model.

.. figure:: _images/mvr_lgr_planview.png
   :name: fig:mvr_lgr_planview

   Plan view of a three-dimensional aquifer system used to test the
   local grid refinement method in MODFLOW 6. A 15\ :math:`\times`\ 15
   horizontal grid discretization is shown for the parent grid and the
   locally refined grid (15\ :math:`\times`\ 18) spacing is equivalent
   to a 45\ :math:`\times`\ 45 discretization over the whole domain. Two
   between-model MVR connections are invoked for the SFR packages
   (parent and child) where the river crosses between domains

.. _example-description-1:

Example description
-------------------

The parent model is 15 rows by 15 columns by 3 layers with the child
model occupying a 5 row by 6 column by 2 layer portion of the parent
grid as shown in figure `3 <#fig:mvr_lgr_planview>`__. The child model
applies a 3:1 refinement to each parent grid cell, including in the
vertical direction, resulting in a 15 row by 18 column by 6 layer child
model. In this way, each parent grid cell is replaced by 27
(:math:`3^3`) child model grid cells. Aquifer properties are uniform and
equal throughout both model domains. Interested readers are referred to
(Mehl and Hill 2013) for additional model details.

.. container::
   :name: tab:ex-gwf-lgr-01

   .. table:: Model Parameters

      ====================================================== =========
      **Parameter**                                          **Value**
      ====================================================== =========
      Number of layers in parent model                       3
      Number of rows in parent model                         15
      Number of columns in parent model                      15
      Parent model column width (:math:`m`)                  102.94
      Parent model row width (:math:`m`)                     68.63
      Number of child model columns per parent model columns 5
      Number of child model rows per parent model rows       5
      Child model column width (:math:`m`)                   20.588
      Child model row width (:math:`m`)                      13.725
      Horizontal hydraulic conductivity (:math:`m/d`)        1.0
      Vertical hydraulic conductivity (:math:`m/d`)          1.0
      ====================================================== =========

Figure `4 <#fig:ex-gwf-lgr>`__ demonstrates that both MVR connections
transfer streamflow to the appropriate stream reach in the receiving
downstream model. This is evidenced by the continuous and smooth
downward trend in streamflow at the first MVR connection as well as by
the uninterrupted upward trend in streamflow at the second MVR
connection.

.. figure:: _images/ex-gwf-lgr.png
   :name: fig:ex-gwf-lgr

   Simulated streamflow and groundwater surface-water exchange for each
   stream reach in figure `3 <#fig:mvr_lgr_planview>`__. The transfer of
   water between adjacent stream reaches in the parent versus child
   models occurs at locations identified by the vertical black lines on
   the plot. Bar widths are indicative of the relative stream reach
   lengths within the host grid cell (figure
   `3 <#fig:mvr_lgr_planview>`__).

One-Dimensional Steady Flow with Transport (MOC3D Problem 1)
============================================================

This problem corresponds to the first problem presented in the MOC3D
report (Konikow, Goode, and Hornberger 1996), which involves the
transport of a dissolved constituent in a steady, one-dimensional flow
field. An analytical solution for this problem is given by (Wexler
1992). This example is simulated with the GWT Model in MODFLOW 6, which
receives flow information from a separate simulation with the GWF Model
in MODFLOW 6. Results from the GWT Model are compared with the results
from the (Wexler 1992) analytical solution.

.. _example-description-2:

Example description
-------------------

The parameters used for this problem are listed in
table `3 <#tab:ex-gwt-moc3d-p01-01>`__. The model grid for this problem
consists of one layer, 120 rows, and 1 columns. The top for each cell is
assigned a value of 1.0 :math:`cm` and the bottom is assigned a value of
zero. DELR is set to 1.0 :math:`cm` and DELC is specified with a
constant value of 0.1 :math:`cm`. The simulation consists of one stress
period that is 120 :math:`s` in length, and the stress period is divided
into 240 equally sized time steps. By using a uniform porosity value of
0.1, a velocity value of 0.1 :math:`cm/s` results from the injection of
water at a rate of 0.001 :math:`cm^3/s` into the first cell. The last
cell is assigned a constant head with a value of zero, though this value
is not important as the cells are marked as being confined. The
concentration of the injected water is assigned a value of 1.0, and any
water that leaves through the constant-head cell leaves with the
simulated concentration of the water in that last cell. Advection is
solved using the TVD scheme to reduce numerical dispersion.

.. container::
   :name: tab:ex-gwt-moc3d-p01-01

   .. table:: Model Parameters

      ========================================== =========
      **Parameter**                              **Value**
      ========================================== =========
      Number of periods                          1
      Number of layers                           1
      Number of rows                             1
      Number of columns                          122
      Length of system (:math:`cm`)              12.0
      Column width (:math:`cm`)                  0.1
      Row width (:math:`cm`)                     0.1
      Top of the model (:math:`cm`)              1.0
      Layer bottom elevation (:math:`cm`)        0
      Specific discharge (:math:`cm s^{-1}`)     0.1
      Hydraulic conductivity (:math:`cm s^{-1}`) 0.01
      Porosity of mobile domain (unitless)       0.1
      Simulation time (:math:`s`)                120.0
      Source concentration (unitless)            1.0
      Initial concentration (unitless)           0.0
      ========================================== =========

Example Scenarios
-----------------

This example problem consists of several different scenarios, as listed
in table `4 <#tab:ex-gwt-moc3d-p01-scenario>`__. Two different levels of
dispersion were simulated, and these simulations are referred to as the
low dispersion case and the high dispersion case. The low dispersion
case has a dispersion coefficient of 0.01 :math:`cm^2/s`, which, for the
specified velocity, corresponds to a dispersivity value of 0.1
:math:`cm`. The high-dispersion case has a dispersion coefficient of 0.1
:math:`cm^2/s`, which corresponds to a dispersivity value of 1.0
:math:`cm`.

.. container::
   :name: tab:ex-gwt-moc3d-p01-scenario

   .. table:: Model Scenario Parameters

      +--------------+-------------------+-------------------+-----------+
      | **Scenario** | **Scenario Name** | **Parameter**     | **Value** |
      +==============+===================+===================+===========+
      | 1            | ex-gwt-moc3d-p01a | longitud          | 0.1       |
      |              |                   | inal dispersivity |           |
      |              |                   | (:math:`cm`)      |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | r                 | 1.0       |
      |              |                   | etardation factor |           |
      |              |                   | (unitless)        |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay rate        | 0.0       |
      |              |                   | (:math:`s^{-1}`)  |           |
      +--------------+-------------------+-------------------+-----------+
      | 2            | ex-gwt-moc3d-p01b | longitud          | 1.0       |
      |              |                   | inal dispersivity |           |
      |              |                   | (:math:`cm`)      |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | r                 | 1.0       |
      |              |                   | etardation factor |           |
      |              |                   | (unitless)        |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay rate        | 0.0       |
      |              |                   | (:math:`s^{-1}`)  |           |
      +--------------+-------------------+-------------------+-----------+
      | 3            | ex-gwt-moc3d-p01c | longitud          | 1.0       |
      |              |                   | inal dispersivity |           |
      |              |                   | (:math:`cm`)      |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | r                 | 2.0       |
      |              |                   | etardation factor |           |
      |              |                   | (unitless)        |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay rate        | 0.0       |
      |              |                   | (:math:`s^{-1}`)  |           |
      +--------------+-------------------+-------------------+-----------+
      | 4            | ex-gwt-moc3d-p01d | longitud          | 1.0       |
      |              |                   | inal dispersivity |           |
      |              |                   | (:math:`cm`)      |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | r                 | 1.0       |
      |              |                   | etardation factor |           |
      |              |                   | (unitless)        |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay rate        | 0.01      |
      |              |                   | (:math:`s^{-1}`)  |           |
      +--------------+-------------------+-------------------+-----------+

Scenario Results
~~~~~~~~~~~~~~~~

For the first scenario with a relatively small dispersivity value (0.1
:math:`cm`), plots of concentration versus time and concentration versus
distance are shown in figures `5 <#fig:ex-gwt-moc3d-p01a-ct>`__ and
 `6 <#fig:ex-gwt-moc3d-p01a-cd>`__, respectively.
Figure `5 <#fig:ex-gwt-moc3d-p01a-ct>`__ can be compared to figure 18 in
(Konikow, Goode, and Hornberger 1996). The three separate concentration
versus time curves represent the three different distances (0.05, 4.05,
and 11.05 cm). For this low-dispersion case, the MODFLOW 6 solution is
in relatively good agreement with the analytical solution, but some
slight differences are observed. These differences are due primarily to
limitations with the second-order TVD scheme implemented in MODFLOW 6,
which can suffer from numerical dispersion.

.. figure:: _images/ex-gwt-moc3d-p01a-ct.png
   :name: fig:ex-gwt-moc3d-p01a-ct

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the low dispersion case. Circles are for the GWT Model results; the
   lines represent the analytical solution by (Wexler 1992). Results are
   shown for three different distances (0.05, 4.05, and 11.05 :math:`cm`
   from the end of the first cell). Every fifth time step is shown for
   the MODFLOW 6 results.

.. figure:: _images/ex-gwt-moc3d-p01a-cd.png
   :name: fig:ex-gwt-moc3d-p01a-cd

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the low dispersion case. Circles are for the GWT Model results; the
   lines represent the analytical solution by (Wexler 1992). Results are
   shown for three different times (6, 60, and 120 :math:`s`). Every
   fifth cell is shown for the MODFLOW 6 results.

For the second scenario, which has a relatively large dispersivity value
(1.0 :math:`cm`), plots of concentration versus time and concentration
versus distance are shown in figures `7 <#fig:ex-gwt-moc3d-p01b-ct>`__
and  `8 <#fig:ex-gwt-moc3d-p01b-cd>`__, respectively. For the
high-dispersion case, the results from the MODFLOW 6 simulation are in
better agreement with the analytical solution than for the low
dispersion case.

.. figure:: _images/ex-gwt-moc3d-p01b-ct.png
   :name: fig:ex-gwt-moc3d-p01b-ct

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the high dispersion case. Circles are for the GWT Model results; the
   lines represent the analytical solution by (Wexler 1992). Results are
   shown for three different distances (0.05, 4.05, and 11.05 :math:`cm`
   from the end of the first cell). Every fifth time step is shown for
   the MODFLOW 6 results.

.. figure:: _images/ex-gwt-moc3d-p01b-cd.png
   :name: fig:ex-gwt-moc3d-p01b-cd

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the high dispersion case. Circles are for the GWT Model results; the
   lines represent the analytical solution by (Wexler 1992). Results are
   shown for three different times (6, 60, and 120 :math:`s`). Every
   fifth cell is shown for the MODFLOW 6 results.

For the remaining scenarios, the results from MODFLOW 6 are compared
with the (Wexler 1992) analytical solution for two variations of the
high-dispersion case. The effects sorption are included in scenario 3
and the effects of decay are included in scenario 4. Plots of
concentration versus time and concentration versus distance for a
simulation with a retardation factor of 2.0 are shown in
figures `9 <#fig:ex-gwt-moc3d-p01c-ct>`__ and
 `10 <#fig:ex-gwt-moc3d-p01c-cd>`__, respectively. Plots of
concentration versus time and concentration versus distance for a
simulation with a decay rate of 0.01 :math:`s^{-1}` are shown in
figures `11 <#fig:ex-gwt-moc3d-p01d-ct>`__ and
 `12 <#fig:ex-gwt-moc3d-p01d-cd>`__, respectively.

.. figure:: _images/ex-gwt-moc3d-p01c-ct.png
   :name: fig:ex-gwt-moc3d-p01c-ct

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the high dispersion case and a retardation factor of 2. Circles are
   for the GWT Model results; the lines represent the analytical
   solution by (Wexler 1992). Results are shown for three different
   distances (0.05, 4.05, and 11.05 :math:`cm` from the end of the first
   cell). Every fifth time step is shown for the MODFLOW 6 results.

.. figure:: _images/ex-gwt-moc3d-p01c-cd.png
   :name: fig:ex-gwt-moc3d-p01c-cd

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the high dispersion case and a retardation factor of 2. Circles are
   for the GWT Model results; the lines represent the analytical
   solution by (Wexler 1992). Results are shown for three different
   times (6, 60, and 120 :math:`s`). Every fifth cell is shown for the
   MODFLOW 6 results.

.. figure:: _images/ex-gwt-moc3d-p01d-ct.png
   :name: fig:ex-gwt-moc3d-p01d-ct

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the high dispersion case and a decay rate of 0.01 :math:`s^{-1}`.
   Circles are for the GWT Model results; the lines represent the
   analytical solution by (Wexler 1992). Results are shown for three
   different distances (0.05, 4.05, and 11.05 :math:`cm` from the end of
   the first cell). Every fifth time step is shown for the MODFLOW 6
   results.

.. figure:: _images/ex-gwt-moc3d-p01d-cd.png
   :name: fig:ex-gwt-moc3d-p01d-cd

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for one-dimensional flow with transport for
   the high dispersion case and a decay rate of 0.01 :math:`s^{-1}`.
   Circles are for the GWT Model results; the lines represent the
   analytical solution by (Wexler 1992). Results are shown for three
   different times (6, 60, and 120 :math:`s`). Every fifth cell is shown
   for the MODFLOW 6 results.

Three-Dimensional Steady Flow with Transport (MOC3D Problem 2)
==============================================================

This problem corresponds to the second problem presented in the MOC3D
report (Konikow, Goode, and Hornberger 1996), which involves the
transport of a dissolved constituent in a steady, three-dimensional flow
field. An analytical solution for this problem is given by (Wexler
1992). This example is simulated with the GWT Model in MODFLOW 6, which
receives flow information from a separate simulation with the GWF Model
in MODFLOW 6. Results from the GWT Model are compared with the results
from the (Wexler 1992) analytical solution.

.. _example-description-3:

Example description
-------------------

(Wexler 1992) presents an analytical solution for three dimensional
solute transport from a point source in a one-dimensional flow field. As
described by (Konikow, Goode, and Hornberger 1996), only one quadrant of
the three-dimensional domain is represented by the numerical model.
Thus, the solute mass flux specified for the model is one quarter of the
solute mass flux used in the analytical solution.

The parameters used for this problem are listed in
table `5 <#tab:ex-gwt-moc3d-p02-01>`__. The model grid for this problem
consists of 40 layers, 12 rows, and 30 columns. The top for layer 1 is
set to zero, and flat bottoms are assigned to all layers based on a
uniform layer thickness of 0.05 :math:`m`. DELR is set to 3.0 :math:`m`
and DELC is specified with a constant value of 0.5 :math:`m`. The
simulation consists of one stress period that is 400 :math:`d` in
length, and the stress period is divided into 400 equally sized time
steps. Velocity is specified to be 0.1 :math:`m/d` in the x direction
and zero in the y and z directions. The uniform flow field is
represented by specifying a constant inflow rate into all of the cells
in column 1 and by specifying a constant head condition to all of the
cells in column 30. A specified solute flux of 10 grams per day is
specified to the cell in layer 1, row 12, and column 8. Any water that
leaves through the constant-head cell leaves with the simulated
concentration of the water in that last cell. Advection is solved using
the TVD scheme to reduce numerical dispersion. In addition to the
longitudinal dispersion, transverse dispersion is represented with a
different value in the horizontal direction than in the vertical
direction. Because the velocity field is perfectly aligned with the
model grid, there are no cross-dispersion terms and the problem can be
simulated accurately without the need for XT3D.

.. container::
   :name: tab:ex-gwt-moc3d-p02-01

   .. table:: Model Parameters

      ============================================== ==========
      **Parameter**                                  **Value**
      ============================================== ==========
      Number of periods                              1
      Number of layers                               40
      Number of rows                                 12
      Number of columns                              30
      Column width (:math:`m`)                       3
      Row width (:math:`m`)                          0.5
      Layer thickness (:math:`m`)                    0.05
      Top of the model (:math:`m`)                   0.0
      Model bottom elevation (:math:`m`)             -2.0
      Velocity in x-direction (:math:`m d^{-1}`)     0.1
      Hydraulic conductivity (:math:`m d^{-1}`)      0.0125
      Porosity of mobile domain (unitless)           0.25
      Longitudinal dispersivity (:math:`m`)          0.6
      Transverse horizontal dispersivity (:math:`m`) 0.03
      Transverse vertical dispersivity (:math:`m`)   0.006
      Simulation time (:math:`d`)                    400.0
      Solute mass flux (:math:`g d^{-1}`)            2.5
      Source location (layer, row, column)           (1, 12, 8)
      ============================================== ==========

.. _example-results-1:

Example Results
---------------

A comparison of the MODFLOW 6 results with the analytical solution of
(Wexler 1992) is shown for layer 1 in
figure `13 <#fig:ex-gwt-moc3d-p02-map>`__.

.. figure:: _images/ex-gwt-moc3d-p02-map.png
   :name: fig:ex-gwt-moc3d-p02-map

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for three-dimensional flow with transport.
   Results are for the end of the simulation (time=400 :math:`d`) and
   for layer 1. Black lines represent solute concentration contours from
   the analytical solution (Wexler 1992); blue lines represent solute
   concentration contours simulated by MODFLOW 6. An aspect ratio of 4.0
   is specified to enhance the comparison.

Three-Dimensional Steady Flow with Transport (MOC3D Problem 2 with a Triangular Grid)
=====================================================================================

This problem corresponds to the second problem presented in the MOC3D
report (Konikow, Goode, and Hornberger 1996), which involves the
transport of a dissolved constituent in a steady, three-dimensional flow
field. An analytical solution for this problem is given by (Wexler
1992). As for the previous example, this example is simulated with the
GWT Model in MODFLOW 6, which receives flow information from a separate
simulation with the GWF Model in MODFLOW 6. In this example, however, a
triangular grid is used for the flow and transport simulation. Results
from the GWT Model are compared with the results from the (Wexler 1992)
analytical solution.

.. figure:: _images/ex-gwt-moc3d-p02tg-grid.png
   :name: fig:ex-gwt-moc3d-p02tg-grid

   Triangular model grid used for the MODFLOW 6 simulation. Model grid
   is shown using an aspect ratio of 4.

.. _example-description-4:

Example description
-------------------

(Wexler 1992) presents an analytical solution for three dimensional
solute transport from a point source in a one-dimensional flow field. As
described by (Konikow, Goode, and Hornberger 1996), only one quadrant of
the three-dimensional domain is represented by the numerical model.
Thus, the solute mass flux specified for the model is one quarter of the
solute mass flux used in the analytical solution.

The parameters used for this problem are listed in
table `6 <#tab:ex-gwt-moc3d-p02tg-01>`__. The model grid for this
problem consists of 40 layers, 695 cells per layer, and 403 vertices in
a layer. The top for layer 1 is set to zero, and flat bottoms are
assigned to all layers based on a uniform layer thickness of 0.05
:math:`m`. The remaining parameters are set similarly to the previous
simulation with a regular grid, except for in this simulation, the XT3D
method is used for flow and dispersive transport.

.. container::
   :name: tab:ex-gwt-moc3d-p02tg-01

   .. table:: Model Parameters

      ============================================== ==========
      **Parameter**                                  **Value**
      ============================================== ==========
      Number of periods                              1
      Number of layers                               40
      Number of rows                                 12
      Number of columns                              30
      Column width (:math:`m`)                       3
      Row width (:math:`m`)                          0.5
      Layer thickness (:math:`m`)                    0.05
      Top of the model (:math:`m`)                   0.0
      Model bottom elevation (:math:`m`)             -2.0
      Velocity in x-direction (:math:`m d^{-1}`)     0.1
      Hydraulic conductivity (:math:`m d^{-1}`)      0.0125
      Porosity of mobile domain (unitless)           0.25
      Longitudinal dispersivity (:math:`m`)          0.6
      Transverse horizontal dispersivity (:math:`m`) 0.03
      Transverse vertical dispersivity (:math:`m`)   0.006
      Simulation time (:math:`d`)                    400.0
      Solute mass flux (:math:`g d^{-1}`)            2.5
      Source location (layer, row, column)           (1, 12, 8)
      ============================================== ==========

.. _example-results-2:

Example Results
---------------

A comparison of the MODFLOW 6 results with the analytical solution of
(Wexler 1992) is shown for layer 1 in
figure `15 <#fig:ex-gwt-moc3d-p02tg-map>`__.

.. figure:: _images/ex-gwt-moc3d-p02tg-map.png
   :name: fig:ex-gwt-moc3d-p02tg-map

   Concentrations simulated by the MODFLOW 6 GWT Model and calculated by
   the analytical solution for three-dimensional flow with transport.
   Results are for the end of the simulation (time=400 :math:`d`) and
   for layer 1. Black lines represent solute concentration contours from
   the analytical solution (Wexler 1992); blue lines represent solute
   concentration contours simulated by MODFLOW 6. An aspect ratio of 4.0
   is specified to enhance the comparison.

Zero-Order Growth in a Uniform Flow Field (MT3DMS Supplemental Guide Problem 6.3.1)
===================================================================================

This example is for zero-order production in a uniform flow field. It is
based on example problem 6.3.1 described in (Zheng 2010). The problem
consists of a one-dimensional model grid with inflow into the first cell
and outflow through the last cell. This example is simulated with the
GWT Model in MODFLOW 6, which receives flow information from a separate
simulation with the GWF Model in MODFLOW 6. Results from the GWT Model
are compared with the results from a MT3DMS simulation (Zheng 1990) that
uses flows from a separate MODFLOW-2005 simulation (Harbaugh 2005).

.. _example-description-5:

Example description
-------------------

The parameters used for this problem are listed in
table `7 <#tab:ex-gwt-mt3dsupp631-01>`__. The model grid consists of 101
columns, 1 row, and 1 layer. The flow problem is confined and steady
state with an initial head set to the model top. The solute transport
simulation represents transient conditions, which begin with an initial
concentration specified as zero everywhere within the model domain. A
specified flow condition is assigned to the first model cell. For the
source pulse duration, the inflow concentration is specified as one.
Following the source pulse duration the inflowing water is assigned a
concentration of zero. A specified head condition is assigned to the
last model cell. Water exiting the model through the specified head cell
leaves with the simulated concentration of that cell.

.. container::
   :name: tab:ex-gwt-mt3dsupp631-01

   .. table:: Model Parameters

      ================================================ =========
      **Parameter**                                    **Value**
      ================================================ =========
      Number of periods                                2
      Number of layers                                 1
      Number of rows                                   1
      Number of columns                                101
      Column width (:math:`m`)                         0.16
      Row width (:math:`m`)                            1.0
      Top of the model (:math:`m`)                     1.0
      Layer bottom elevation (:math:`m`)               0
      Specific discharge (:math:`md^{-1}`)             0.1
      Longitudinal dispersivity (:math:`m`)            1.0
      Porosity of mobile domain (unitless)             0.37
      Zero-order production rate (:math:`mg/L d^{-1}`) -2.0e-3
      Source duration (:math:`d`)                      160.0
      Simulation time (:math:`t`)                      840.0
      Observation x location (:math:`m`)               8.0
      ================================================ =========

.. _example-results-3:

Example Results
---------------

Simulated concentrations from the MODFLOW 6 GWT Model and MT3DMS are
shown in figure `16 <#fig:ex-gwt-mt3dsupp631>`__. The close agreement
between the simulated concentrations demonstrate the
zero-order-production capabilities implemented in the GWT Model.

.. figure:: _images/ex-gwt-mt3dsupp631.png
   :name: fig:ex-gwt-mt3dsupp631

   Concentrations simulated by the MODFLOW 6 GWT Model and MT3DMS for
   zero-order growth in a uniform flow field.

Zero-Order Growth in a Dual-Domain System (MT3DMS Supplemental Guide Problem 6.3.2)
===================================================================================

This example is for zero-order production in a dual-domain system. It is
based on example problem 6.3.2 described in (Zheng 2010). The problem
consists of a one-dimensional model grid with inflow into the first cell
and outflow through the last cell. This example is simulated with the
GWT Model in MODFLOW 6, which receives flow information from a separate
simulation with the GWF Model in MODFLOW 6. This example is designed to
test the capabilities of the GWT Model to simulate zero-order production
in a dual-domain system with and without sorption. Results from the GWT
Model are compared with the results from a MT3DMS simulation (Zheng
1990) that uses flows from a separate MODFLOW-2005 simulation (Harbaugh
2005). This example was described by (Zheng 2010) who showed that the
results from MT3DMS were in good agreement with an analytical solution.

.. _example-description-6:

Example description
-------------------

The parameters used for this problem are listed in
table `8 <#tab:ex-gwt-mt3dsupp632-01>`__. The model grid consists of 401
columns, 1 row, and 1 layer. The flow problem is confined and steady
state with an initial head set to the model top. The solute transport
simulation represents transient conditions, which begin with an initial
concentration specified as zero everywhere within the model domain. A
specified flow condition is assigned to the first model cell. For the
source pulse duration, a specified concentration with a value of one is
assigned to the first model cell. Following the source pulse duration
the specified concentration in the first cell is zero. A specified head
condition is assigned to the last model cell. Water exiting the model
through the specified head cell leaves with the simulated concentration
of that cell.

.. container::
   :name: tab:ex-gwt-mt3dsupp632-01

   .. table:: Model Parameters

      +---------------------------------------------------------+-----------+
      | **Parameter**                                           | **Value** |
      +=========================================================+===========+
      | Number of periods                                       | 2         |
      +---------------------------------------------------------+-----------+
      | Number of layers                                        | 1         |
      +---------------------------------------------------------+-----------+
      | Number of rows                                          | 1         |
      +---------------------------------------------------------+-----------+
      | Number of columns                                       | 401       |
      +---------------------------------------------------------+-----------+
      | Column width (:math:`m`)                                | 2.5       |
      +---------------------------------------------------------+-----------+
      | Row width (:math:`m`)                                   | 1.0       |
      +---------------------------------------------------------+-----------+
      | Top of the model (:math:`m`)                            | 1.0       |
      +---------------------------------------------------------+-----------+
      | Layer bottom elevation (:math:`m`)                      | 0         |
      +---------------------------------------------------------+-----------+
      | Specific discharge (:math:`md^{-1}`)                    | 0.06      |
      +---------------------------------------------------------+-----------+
      | Longitudinal dispersivity (:math:`m`)                   | 10        |
      +---------------------------------------------------------+-----------+
      | Porosity of mobile domain (unitless)                    | 0.2       |
      +---------------------------------------------------------+-----------+
      | Porosity of immobile domain (unitless)                  | 0.05      |
      +---------------------------------------------------------+-----------+
      | Bulk density (:math:`gL^{-1})`                          | 4.0       |
      +---------------------------------------------------------+-----------+
      | First-order mass transfer rate between the mobile and   | 1.0e-3    |
      | immobile domains (:math:`d^{-1}`)                       |           |
      +---------------------------------------------------------+-----------+
      | Fraction of sorption sites in contact with mobile water | 0.8       |
      | (unitless)                                              |           |
      +---------------------------------------------------------+-----------+
      | Source duration (:math:`d`)                             | 1000      |
      +---------------------------------------------------------+-----------+
      | Simulation time (:math:`t`)                             | 10000     |
      +---------------------------------------------------------+-----------+
      | Observation x location (:math:`m`)                      | 200.0     |
      +---------------------------------------------------------+-----------+

.. _example-scenarios-1:

Example Scenarios
-----------------

This example problem consists of several different scenarios, as listed
in table `9 <#tab:ex-gwt-mt3dsupp632-scenario>`__. The first two
scenarios represent zero-order growth when sorbtion is active. Sorbtion
is not active in the last scenario. For all three scenarios, there is
mass transfer between the mobile domain and the immobile domain.

.. container::
   :name: tab:ex-gwt-mt3dsupp632-scenario

   .. table:: Model Scenario Parameters

      +--------------+-------------------+-------------------+-----------+
      | **Scenario** | **Scenario Name** | **Parameter**     | **Value** |
      +==============+===================+===================+===========+
      | 1            | ex                | distrib           | 0.25      |
      |              | -gwt-mt3dsupp632a | ution coefficient |           |
      |              |                   | (:                |           |
      |              |                   | math:`mL g^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay             | 0.0       |
      |              |                   | (:ma              |           |
      |              |                   | th:`g/mL d^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay sorbed      | -0.001    |
      |              |                   | (:ma              |           |
      |              |                   | th:`g/mL d^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      | 2            | ex                | distrib           | 0.25      |
      |              | -gwt-mt3dsupp632b | ution coefficient |           |
      |              |                   | (:                |           |
      |              |                   | math:`mL g^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay             | -0.0005   |
      |              |                   | (:ma              |           |
      |              |                   | th:`g/mL d^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay sorbed      | -0.0005   |
      |              |                   | (:ma              |           |
      |              |                   | th:`g/mL d^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      | 3            | ex                | distrib           | 0.0       |
      |              | -gwt-mt3dsupp632c | ution coefficient |           |
      |              |                   | (:                |           |
      |              |                   | math:`mL g^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay             | -0.001    |
      |              |                   | (:ma              |           |
      |              |                   | th:`g/mL d^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+
      |              |                   | decay sorbed      | 0.0       |
      |              |                   | (:ma              |           |
      |              |                   | th:`g/mL d^{-1}`) |           |
      +--------------+-------------------+-------------------+-----------+

.. _scenario-results-1:

Scenario Results
~~~~~~~~~~~~~~~~

Results from the three scenarios are shown in
figure `17 <#fig:ex-gwt-mt3dsupp632>`__. The close agreement between the
simulated concentrations for the MODFLOW 6 GWT Model and MT3DMS
demonstrate the zero-order growth and immobile-domain transfer
capabilities for MODFLOW 6.

.. figure:: _images/ex-gwt-mt3dsupp632.png
   :name: fig:ex-gwt-mt3dsupp632

   Concentrations simulated by the MODFLOW 6 GWT Model and MT3DMS for
   zero-order growth in a dual-domain system. Circles are for the GWT
   Model results; the lines represent simulated concentrations for
   MT3DMS.

Simulating the Effect of a Recirculating Well (MT3DMS Supplemental Guide Problem 8.2)
=====================================================================================

This example is for a recirculating well. It is based on example problem
8.2 described in (Zheng 2010). The problem consists of a
two-dimensional, one-layer model with flow from left to right. A solute
is introduced into the flow field by an injection well. Downgradient, an
extraction well pumps at the same rate as the injection well. This
extracted water is then injected into two other injection wells. This
example is simulated with the GWT Model in MODFLOW 6, which receives
flow information from a separate simulation with the GWF Model in
MODFLOW 6. Results from the GWT Model are compared with the results from
a MT3DMS simulation (Zheng 1990) that uses flows from a separate
MODFLOW-2005 simulation (Harbaugh 2005).

.. _example-description-7:

Example description
-------------------

The parameters used for this problem are listed in
table `10 <#tab:ex-gwt-mt3dsupp82-01>`__. The model grid consists of 31
rows, 46 columns, and 1 layer. The flow problem is confined and steady
state. The solute transport simulation represents transient conditions,
which begin with an initial concentration specified as zero everywhere
within the model domain.

For the MT3DMS representation of this problem, the Well Package is used
to inject water at a rate of 1 :math:`m^3/d` into model cell (1, 16,
16). Water is extracted at a rate of -1 :math:`m^3/d` from model cell
(1, 16, 21). Two additional wells, located in cells (1, 5, 16) and (1,
27, 16), reinject water at the concentration of the extracted water from
cell (1, 16, 21). The injection rate for each of these reinjection wells
is 0.5 :math:`m^3/d`.

For the MODFLOW 6 representation of this problem, the Multi-Aquifer Well
(MAW) Package is used for these injection and extraction wells, although
because the model is only a single layer, the MAW Package behaves just
like the Well Package. The Water Mover (MVR) Package is used to send
half of the extracted water into each of the reinjection wells. For the
MODFLOW 6 transport simulation, the Multi-Aquifer Transport (MWT)
Package is used to calculate the concentration in each of the well
bores. The Mover Transport (MVT) Package is used to move the solute from
the extraction well to the two reinjection wells based on the simulated
flows. Note that this approach used in MODFLOW 6 is slightly different
than the approach used in MT3DMS, because MODFLOW 6 is calculating the
concentrations in the well boreholes rather than using concentrations
directly from the model cells. By specifying a small well radius for the
MAW Package, the approaches are similar.

.. container::
   :name: tab:ex-gwt-mt3dsupp82-01

   .. table:: Model Parameters

      ============================================== =========
      **Parameter**                                  **Value**
      ============================================== =========
      Number of periods                              1
      Number of layers                               1
      Number of rows                                 31
      Number of columns                              46
      Column width (:math:`m`)                       10.0
      Row width (:math:`m`)                          10.0
      Top of the model (:math:`m`)                   10.0
      Layer bottom elevation (:math:`m`)             0.0
      Hydraulic conductivity (:math:`md^{-1}`)       10.0
      Longitudinal dispersivity (:math:`m`)          10.0
      Transverse horizontal dispersivity (:math:`m`) 3.0
      Transverse vertical dispersivity (:math:`m`)   0.3
      Simulation time (:math:`d`)                    365.0
      Porosity of mobile domain (unitless)           0.3
      ============================================== =========

.. _example-results-4:

Example Results
---------------

Simulated concentrations from MODFLOW 6 and MT3DMS are shown in
figure `18 <#fig:ex-gwt-mt3dsupp82-map>`__. The close agreement between
the simulated concentrations demonstrate the ability of MODFLOW 6 to
simulate the transfer of water and solute using the mover package
capability.

.. figure:: _images/ex-gwt-mt3dsupp82-map.png
   :name: fig:ex-gwt-mt3dsupp82-map

   Concentrations simulated by MODFLOW 6 and MT3DMS for a problem
   involving a recirculating well. This figure can be compared to figure
   8.2 in (Zheng 2010).

One-Dimensional Transport in a Uniform Flow Field (MT3DMS Example Problem 1)
============================================================================

Section 7 of (Zheng and Wang 1999) details a number of test problems
that verify the accuracy of MT3DMS. The first problem presented, titled
"one-dimensional transport in a uniform flow field," compared MT3DMS
solutions to analytical solutions given by (Van Genuchten and Alves
1982). For verifying the accuracy of transport calculations within
MODFLOW6, the transport solutions calculated by MT3DMS serve as the
benchmark to which the MODFLOW 6 solution is compared. The first
1-dimensional simulation solves for advection only. The second model
permutation uses both advection and dispersion to verify MODFLOW 6
results. Next, the accuracy of MODFLOW 6 is verified when advection,
dispersion, and some simple checmial reactions represented with sorption
processes are used. The fourth and final model permutation adds solute
decay to the previous model setup. As of the first release of MODFLOW 6
with transport capabilities, a linear isotherm is the only option
available for simulating sorption. Arbitrary values of bulk density and
distribution coefficient are uniformly applied to the entire model
domain to achieve the indicated retardation factor. The following table
summarizes how the four simulations incrementally increase model
complexity.

.. container::
   :name: tab:ex-gwt-mt3dms-p01-scenario

   .. table:: Model Scenario Parameters

      ============ ================== ======================== =========
      **Scenario** **Scenario Name**  **Parameter**            **Value**
      ============ ================== ======================== =========
      1            ex-gwt-mt3dms-p01a dispersivity (:math:`m`) 0.0
      \                               retardation (unitless)   1.0
      \                               decay (:math:`d^{-1}`)   0.0
      2            ex-gwt-mt3dms-p01b dispersivity (:math:`m`) 10.0
      \                               retardation (unitless)   1.0
      \                               decay (:math:`d^{-1}`)   0.0
      3            ex-gwt-mt3dms-p01c dispersivity (:math:`m`) 10.0
      \                               retardation (unitless)   5.0
      \                               decay (:math:`d^{-1}`)   0.0
      4            ex-gwt-mt3dms-p01d dispersivity (:math:`m`) 10.0
      \                               retardation (unitless)   5.0
      \                               decay (:math:`d^{-1}`)   0.002
      ============ ================== ======================== =========

.. _example-description-8:

Example description
-------------------

All four model scenarios have 101 columns, 1 row, and 1 layer. The first
and last columns use constant-head boundaries to simulate steady flow in
confined conditions. Because the analytical solution assumes an infinite
1-dimentional flow field, the last column is set far enough from the
source to avoid interfering with the final solution after 2,000 days.
Initially, the model domain is devoid of solute; however, the first
column uses a constant concentration boundary condition to ensure that
water entering the simulation has a unit concentration of 1. Additional
model parameters are shown in table `12 <#tab:ex-gwt-mt3dms-p01-01>`__.

.. container::
   :name: tab:ex-gwt-mt3dms-p01-01

   .. table:: Model Parameters

      =============================================== =========
      **Parameter**                                   **Value**
      =============================================== =========
      Number of periods                               1
      Number of layers                                1
      Number of columns                               101
      Number of rows                                  1
      Column width (:math:`m`)                        10.0
      Row width (:math:`m`)                           1.0
      Top of the model (:math:`m`)                    0.0
      Layer bottom elevations (:math:`m`)             -1.0
      Porosity                                        0.25
      Simulation time (:math:`days`)                  2000
      Horizontal hydraulic conductivity (:math:`m/d`) 1.0
      =============================================== =========

.. _example-results-5:

Example Results
---------------

Currently no options are available with MODFLOW 6 for simulating solute
transport using particle tracking methods [referred to as Method of
Characteristics (MOC) in the MT3DMS manual (Zheng and Wang 1999). Thus,
the MODFLOW 6 solution is compared to an MT3DMS solution that uses the
third-order total variation diminishing (TVD) option for solving the
advection-only problem rather than invoking one of the MOC options
available within MT3DMS. Owing to different approaches between the two
codes, namely TVD scheme of MT3DMS and the second-order approach of
MODFLOW 6, differences between the two solutions and reflected in
figure `[fig:ex-gwt-mt3dms-p01a] <#fig:ex-gwt-mt3dms-p01a>`__ are
expected. However, the differences are within acceptable tolerances.

The comparison of the MT3DMS andMODFLOW 6 solutions for problem 1a, an
advection dominated problem, represents an end-member test as the
migrating concentration front is sharp (i.e., discontinuous). In
technical terms, the grid Peclet number is infinity for this problem
(:math:`P_e` = :math:`v\Delta x/D_{xx}` =
:math:`\Delta x`/:math:`\alpha_L` = :math:`\infty`).

.. figure:: _images/ex-gwt-mt3dms-p01a.png
   :name: fig:ex-gwf-mt3dms-p01a

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for a
   one-dimensional advection dominated test problem. The analytical
   solution for this problem was originally given in (Van Genuchten and
   Alves 1982) and is not shown here

A comparison of MT3DMS and MODFLOW 6 for scenario 2 in the MT3DMS manual
represents a more common situation whereby dispersion acts to spread or
smooth the advancing concentration front. For this problem, the
dispersion term :math:`\alpha_L` is set equal to the length of the grid
cell in the direction of flow, 10 cm, resulting in a Peclet number equal
to one (:math:`P_e` = :math:`v\Delta x/D_{xx} = 10/10 = 1`). Owing to
the presense of dispersion, the finite-difference solutions employed by
both MT3DMS and MODFLOW 6 for this problem are more accurate, and as a
result are in closer agreement (figure
`20 <#fig:ex-gwt-mt3dms-p01b>`__).

.. figure:: _images/ex-gwt-mt3dms-p01b.png
   :name: fig:ex-gwt-mt3dms-p01b

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for a
   one-dimensional test problem with dispersion (:math:`\alpha_L` = 10
   :math:`cm`). The analytical solution for this problem was originally
   given in (Van Genuchten and Alves 1982) and is not shown here

The third comparison for the one-dimensional transport in a steady flow
field includes uses the same dispersion specified for the second
scenario, but adds retardation. For this problem, retardation slows the
advance of the migrating concentration front by simulating sorption of
the dissolved solute onto the matrix material through which the fluid is
moving. In this way, dissolved mass is transferred from the aqueous
phase to the solid phase. Appropriate reaction package parameter values
are determined for obtaining the specified retardation (5.0) within the
code. Figure `21 <#fig:ex-gwt-mt3dms-p01c>`__ shows a close match
between MT3DMS and MODFLOW 6. In addition,
figure `21 <#fig:ex-gwt-mt3dms-p01c>`__ also shows that after 2,000
days, the concentration front did not advance as far as shown in
figure `20 <#fig:ex-gwt-mt3dms-p01b>`__.

.. figure:: _images/ex-gwt-mt3dms-p01c.png
   :name: fig:ex-gwt-mt3dms-p01c

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for a
   one-dimensional test problem with dispersion (:math:`\alpha_L` = 10
   :math:`cm`) and retardation (:math:`R` = :math:`5.0`). The analytical
   solution for this problem was originally given in (Van Genuchten and
   Alves 1982) and is not shown here

The final comparison for the fourth scenario of problem 1 adds decay to
the dispersion and retardation simulated in the third scenario. For this
case, decay represents the irreversible loss of mass from both the
aqueous and sorbed phases, further stunting the advance of the migrating
concentration front `22 <#fig:ex-gwt-mt3dms-p01d>`__.

.. figure:: _images/ex-gwt-mt3dms-p01d.png
   :name: fig:ex-gwt-mt3dms-p01d

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for a
   one-dimensional test problem with dispersion (:math:`\alpha_L` = 10
   :math:`cm`), retardation (:math:`R` = :math:`5.0`) and decay
   (:math:`\lambda` = 0.002 :math:`d^{-1}`). The analytical solution for
   this problem was originally given in (Van Genuchten and Alves 1982)
   and is not shown here

Two-Dimensional Transport in a Uniform Flow Field (MT3DMS Example Problem 3)
============================================================================

The second example problem in (Zheng and Wang 1999) verifies accurate
simulation of nonlinear and nonequilibrium sorption by comparing results
to corresponding analytical solutions. As neither capability is included
in the initial release of the transport process written forMODFLOW 6,
the next MT3DMS-MODFLOW 6 transport comparison is the third problem
appearing in (Zheng and Wang 1999), titled, "two-dimensional transport
in a uniform flow field." In contrast to the first demonstrated test
problem, transport is simulated in two dimensions with dispersion but no
reactions. An analytical solution for this problem was originally
published in (Wilson and Miller 1978). Two assumptions that make the
analytical solution possible are that (1) the aquifer is areally
infinite and relatively thin to support the assumption that
instantaneous mixing occurs in the vertical direction, and (2) that
compared to the ambient flow field, the injection rate is insignificant.

.. _example-description-9:

Example description
-------------------

Steady uniform flow enters the left edge of a numerical grid with 31
rows, 46 columns, and 1 layer through a constant head boundary and exits
along the right edge. Constant heads are selected to ensure the
hydraulic gradient matches with the analytical solution. The other
boundaries are all no flow. Boundaries are sufficiently far away from
the injection well where the contaminant is released so as not to
interfere with the final solution after 365 days.
Table `13 <#tab:ex-gwt-mt3dms-p03-01>`__ summarizes model setup:

.. container::
   :name: tab:ex-gwt-mt3dms-p03-01

   .. table:: Model Parameters

      ================================================ =========
      **Parameter**                                    **Value**
      ================================================ =========
      Number of layers                                 1
      Number of rows                                   31
      Number of columns                                46
      Column width (:math:`m`)                         10.0
      Row width (:math:`m`)                            10.0
      Layer thickness (:math:`m`)                      10.0
      Top of the model (:math:`m`)                     0.0
      Porosity                                         0.3
      Simulation time (:math:`days`)                   365
      Horizontal hydraulic conductivity (:math:`m/d`)  1.0
      Volumetric injection rate (:math:`m^3/d`)        1.0
      Concentration of injected water (:math:`mg/L`)   1000.0
      Longitudinal dispersivity (:math:`m`)            10.0
      Ratio of transverse to longitudinal dispersivity 0.3
      ================================================ =========

After 365 days, the MODFLOW 6 solution aligns well with the MT3DMS
solution `[fig:mt3dms_p03] <#fig:mt3dms_p03>`__. In addition to the good
agreement that was seen in the first MT3DMS test problem, the current
comparison confirms that the lateral dispersion is accurately simulated
within MODFLOW 6.

.. figure:: _images/ex-gwt-mt3dms-p03.png
   :name: fig:ex-gwf-mt3d-p03

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for a
   two-dimensional advection-dispersion test problem. The analytical
   solution for this problem was originally given in (Wilson and Miller
   1978) and is not shown here

Two-Dimensional Transport in a Diagonal Flow Field (MT3DMS Example Problem 4)
=============================================================================

The third demonstrated MT3DMS-MODFLOW 6 transport comparson is for a
two-dimensional transport in a diagonal flow field. This problem is
similar to the preceding problem with two important changes. First, the
flow direction is now oriented at a 45-degree angle relative to rows and
columns of the numerical grid. Owing the use of MT3DMS for comparison,
the MODFLOW 6 solution uses the traditional DIS package (not DISU or
DISV). The second notable change is that the number of rows and columns
has been expanded in order to accomodate a longer simulation period of
1,000 days. Because of the orientation of the flow field relative to the
model grid, and the sharpness of the migrating concentration front, this
test problem presents a challenging set of conditions to simulate. Three
scenarios test alternative advection formulations, as summarized in
table `14 <#tab:ex-gwt-mt3dms-p04-scenario>`__

.. container::
   :name: tab:ex-gwt-mt3dms-p04-scenario

   .. table:: Model Scenario Parameters

      ============ ================== ================= =========
      **Scenario** **Scenario Name**  **Parameter**     **Value**
      ============ ================== ================= =========
      1            ex-gwt-mt3dms-p04a xt3d (unitless)   True
      \                               mixelm (unitless) 0
      2            ex-gwt-mt3dms-p04b xt3d (unitless)   True
      \                               mixelm (unitless) -1
      3            ex-gwt-mt3dms-p04c xt3d (unitless)   True
      \                               mixelm (unitless) 1
      ============ ================== ================= =========

Model parameter values for this problem are provided in
table `15 <#tab:ex-gwt-mt3dms-p04-01>`__.

.. container::
   :name: tab:ex-gwt-mt3dms-p04-01

   .. table:: Model Parameters

      ================================================== =========
      **Parameter**                                      **Value**
      ================================================== =========
      Number of layers                                   1
      Number of rows                                     100
      Number of columns                                  100
      Column width (:math:`m`)                           10.0
      Row width (:math:`m`)                              10.0
      Layer thickness (:math:`m`)                        1.0
      Top of the model (:math:`m`)                       0.0
      Porosity                                           0.14
      Simulation time (:math:`days`)                     365
      Horizontal hydraulic conductivity (:math:`m/d`)    1.0
      Volumetric injection rate (:math:`m^3/d`)          0.01
      Concentration of injected water (:math:`mg/L`)     1000.0
      Longitudinal dispersivity (:math:`m`)              2.0
      Ratio of transverse to longitudinal dispersitivity 0.1
      Molecular diffusion coefficient (:math:`m^2/d`)    1.0e-9
      ================================================== =========

The same analytical solution used in the previous problem can be used
for this problem after applying the necessary updates to select
parameters - most noteably the dispersion and porosity terms and that an
inter-model comparison is drawn after 1,000 days instead of 365 days.
Figure 36 in (Zheng and Wang 1999) shows four different solutions for
this problem: (1) analytical, (2) Method of Characteristics (MOC), (3)
upstream finite difference (FD), and (4) Total Variation Dimishing (TVD)
or “ULTIMATE” scheme. Both the MOC and TVD solutions demonstrate a
reasonable agreement with the analytical solution. However, the upstream
finite difference solution reflects considerably more spread from
simulation of too much dispersion - in this case numerical dispersion
instead of hydrodynamic dispersion.

The MODFLOW 6 transport solution is compared to all three numerical
solutions (FD, TVD, and MOC) presented in (Zheng and Wang 1999). The
first comparison shows complete agreement between MT3DMS and the
MODFLOW 6 transport solution when the finite difference approach is
applied (figure `24 <#fig:ex-gwt-mt3dms-p04a>`__).

.. figure:: _images/ex-gwt-mt3dms-p04a.png
   :name: fig:ex-gwt-mt3dms-p04a

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for
   two-dimensional transport in a diagonal flow field. Both models are
   using their respective finite difference solutions without the TVD
   option

Figure `25 <#fig:ex-gwt-mt3dms-p04b>`__ shows a comparison between the
MT3DMS and MODFLOW 6 solution with the respective TVD options for each
model activated. Owing to the fact that MT3DMS uses a third-order TVD
scheme while MODFLOW 6 uses a second-order scheme, differences between
the two solutions are expected.

.. figure:: _images/ex-gwt-mt3dms-p04b.png
   :name: fig:ex-gwt-mt3dms-p04b

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for
   two-dimensional transport in a diagonal flow field. Both models are
   using their respective finite difference solutions with the use of
   the TVD option, which serves as the main difference with results
   displayed in figure `25 <#fig:ex-gwt-mt3dms-p04b>`__

The third model comparison shows the largest difference between the two
solutions (figure `26 <#fig:ex-gwt-mt3dms-p04c>`__). Because the MOC
solution is the closest facsimile of the analytical solution, comparison
of MODFLOW 6 with the MT3DMS MOC solution is as close to a comparison
with the analytical solution as will be shown for the current set of
model runs.

.. figure:: _images/ex-gwt-mt3dms-p04c.png
   :name: fig:ex-gwt-mt3dms-p04c

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for
   two-dimensional transport in a diagonal flow field. Here, MT3DMS is
   using a MOC technique to find a solution while MODFLOW 6 uses finite
   difference without TVD activated.

Two-Dimensional Transport in a Radial Flow Field (MT3DMS Example Problem 5)
===========================================================================

The next example problem tests two-dimensional transport in a radial
flow field. The radial flow field is established by injecting water in
the center of the model domain (row 16, column 16) and allowing it to
flow outward toward the model perimeter. No regional groundwater flow
gradient exists as in some of the previous comparisons with MT3DMS.
Constant head cells located around the perimeter of the model drain
water and solute from the simulation domain. Solute enters the model
domain through the injection well with a unit concentration. The
starting concentration is zero across the entire domain. Flow remains
steady and confined throughout the 27 day simulation period. The aquifer
is homogenous, isotropic, and boundaries are sufficiently far from the
injection well to avoid solute reaching the boundary during the
simulation interval. Table `16 <#tab:ex-gwt-mt3dms-p05-01>`__ summarizes
many of the model inputs:

.. container::
   :name: tab:ex-gwt-mt3dms-p05-01

   .. table:: Model Parameters

      ================================================== =========
      **Parameter**                                      **Value**
      ================================================== =========
      Number of layers                                   1
      Number of rows                                     31
      Number of columns                                  31
      Column width (:math:`m`)                           10.0
      Row width (:math:`m`)                              10.0
      Layer thickness (:math:`m`)                        1.0
      Top of the model (:math:`m`)                       0.0
      Porosity                                           0.3
      Simulation time (:math:`days`)                     27
      Horizontal hydraulic conductivity (:math:`m/d`)    1.0
      Volumetric injection rate (:math:`m^3/d`)          100.0
      Concentration of injected water (:math:`mg/L`)     1.0
      Longitudinal dispersivity (:math:`m`)              10.0
      Ratio of transverse to longitudinal dispersitivity 1.0
      Molecular diffusion coefficient (:math:`m^2/d`)    1.0e-9
      ================================================== =========

An analytical solution for this problem was originally given in (Moench
and Ogata 1981). The MT3DMS solution with the TVD option activated most
closely matched the analytical solution. Therefore the TVD option is
activated in both MT3DMS and MODFLOW 6 for verifying the transport
solution. Figure `27 <#fig:ex-gwt-mt3dms-p05-xsec>`__ shows a slight
under simulation of the outward spread of solute in the MODFLOW 6
solution compared to MT3DMS. Figure
`28 <#fig:ex-gwt-mt3dms-p05-planView>`__ shows close agreement among the
MT3DMS and MODFLOW 6 isoconcentration contours with the TVD advection
scheme activated.

.. figure:: _images/ex-gwt-mt3dms-p05-xsec.png
   :name: fig:ex-gwt-mt3dms-p05-xsec

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for a
   point source in a two-dimensional radial flow field simulation. The
   thick black line in figure `28 <#fig:ex-gwt-mt3dms-p05-planView>`__
   shows the location of this profile view of concentrations. The
   analytical solution for this problem was originally given in (Moench
   and Ogata 1981) and is not shown here

.. figure:: _images/ex-gwt-mt3dms-p05-planView.png
   :name: fig:ex-gwt-mt3dms-p05-planView

   Comparison of the MT3DMS and MODFLOW 6 numerical solutions for
   two-dimensional transport in a radial flow field. The thick black
   line shows the location of the concentration profile shown in
   figure `27 <#fig:ex-gwt-mt3dms-p05-xsec>`__. Both models are using
   their respective finite difference solutions with the use of the TVD
   option.

Concentration at an Injection/Extraction Well (MT3DMS Example Problem 6)
========================================================================

In this example problem, concentrations are compared between MODFLOW 6
and MT3D-USGS at an injection/extraction well. The well is fully
penetrating in a confined aquifer and injects contaminated water for a
period of 2.5 years at a rate of 1 :math:`ft^3/sec`. At the end of 2.5
years, the injection well is reversed and begins pumping (extracting)
contaminated groundwater for a period of 7.5 years, also at the rate of
1 :math:`ft^3/sec`. (El‐Kadi 1988) was the first to develop the test
problem which was later used by (Zheng 1993) to test ongoing
method-of-characteristics (MOC) developments. The model boundary is
placed far enough away from the injection/extraction well to ensure no
solute exits the model domain during the injection period. Moreover,
steady flow conditions are reached immediately during the injection and
extraction stress periods. Problem specifics are provided in
table `17 <#tab:ex-gwt-mt3dms-p06-01>`__.

.. container::
   :name: tab:ex-gwt-mt3dms-p06-01

   .. table:: Model Parameters

      ===================================================== =========
      **Parameter**                                         **Value**
      ===================================================== =========
      Number of layers                                      1
      Number of rows                                        31
      Number of columns                                     31
      Column width (:math:`ft`)                             900.0
      Row width (:math:`ft`)                                900.0
      Layer thickness (:math:`ft`)                          20.0
      Top of the model (:math:`ft`)                         0.0
      Porosity                                              0.35
      Length of the injection period (:math:`years`)        2.5
      Length of the extraction period (:math:`years`)       7.5
      Horizontal hydraulic conductivity (:math:`ft/d`)      432.0
      Volumetric injection rate (:math:`ft^3/d`)            1.0
      Relative concentration of injected water (:math:`\%`) 100.0
      Longitudinal dispersivity (:math:`ft`)                100.0
      Ratio of transverse to longitudinal dispersitivity    1.0
      ===================================================== =========

An analytical solution for this problem was originally given in (Gelhar
and Collins 1971). Because this is an advection dominated problem, both
numerical solutions invoke their TVD schemes. The MODFLOW 6 solution
shows a quicker rise in concentration at the well site than does the
MT3D-USGS solution (figure `29 <#fig:ex-gwt-mt3dms-p06>`__).

.. figure:: _images/ex-gwt-mt3dms-p06.png
   :name: fig:ex-gwt-mt3dms-p06

   Comparison of the MT3D-USGS and MODFLOW 6 numerical solutions at an
   injection/extraction well. The analytical solution for this problem
   was originally given in (Gelhar and Collins 1971) and is not shown
   here

Three-Dimensional Transport in a Uniform Flow Field (MT3DMS Example Problem 7)
==============================================================================

In a previous problem titled “two-dimensional transport in a uniform
flow field” concentrations were compared between MODFLOW 6 and MT3D-USGS
for a relatively thin aquifer (10 :math:`m`) wherein instantaneous
vertical mixing was assumed. In order to test transport simulation in a
thicker aquifer, where all three spatial dimensions are required to
adequately simulate the movement of solute, the current problem was
devised. (Hunt 1978) provides an analytical solution, which is used to
verify MT3DMS in (Zheng and Wang 1999), but not shown here. Instead,
only the MODFLOW 6 and MT3DMS solutions are compared here.

Problem dimensions and aquifer properties are given in
table `18 <#tab:ex-gwt-mt3dms-p07-01>`__. The point source is located in
layer 7, row 8, and column 3

.. container::
   :name: tab:ex-gwt-mt3dms-p07-01

   .. table:: Model Parameters

      ================================================== =========
      **Parameter**                                      **Value**
      ================================================== =========
      Number of layers                                   8
      Number of rows                                     15
      Number of columns                                  21
      Column width (:math:`m`)                           10.0
      Row width (:math:`m`)                              10.0
      Layer thickness (:math:`m`)                        10.0
      Top of the model (:math:`m`)                       0.0
      Porosity                                           0.2
      Horizontal hydraulic conductivity (:math:`m/d`)    0.5
      Volumetric injection rate (:math:`m^3/d`)          0.5
      Longitudinal dispersivity (:math:`m`)              10.0
      Ratio of transverse to longitudinal dispersitivity 0.3
      Ratio of vertical to longitudinal dispersitivity   0.3
      Simulation time (:math:`days`)                     100.0
      ================================================== =========

An analytical solution for this problem was originally given in (Hunt
1978). Both numerical solutions invoke their respective TVD schemes.
Moreover, MODFLOW 6 is using the XT3D package for simulating dispersion.
The MODFLOW 6 solution shows great agreement with the MT3DMS calculated
concentrations for the three layers displayed in
figure `29 <#fig:ex-gwt-mt3dms-p06>`__.

.. figure:: _images/ex-gwt-mt3dms-p07.png
   :name: fig:ex-gwt-mt3dms-p07

   Comparison of the MT3D-USGS and MODFLOW 6 numerical solutions for
   three-dimensional transport in a uniform flow field. The analytical
   solution for this problem was originally given in (Hunt 1978) and is
   not shown here

Stream-Lake Interaction with Solute Transport (SFR1 Manual Test Problem 2)
==========================================================================

This problem is based on the stream-aquifer interaction problem
described as test 2 by (Prudic, Konikow, and Banta 2004). (Prudic,
Konikow, and Banta 2004) designed their test 2 problem by modifying a
variant originally described by (Merritt and Konikow 2000). The
description in the text and the figures presented here are largely based
on the text and figures presented by (Prudic, Konikow, and Banta 2004).
The purpose for including this problem here is to demonstrate the use of
MODFLOW 6 to simulate solute transport through a coupled system
consisting of an aquifer, streams, and lakes. The example requires
accurate simulation of transport within the streams and lakes and also
between the surface water features and the underlying aquifer.

.. _example-description-10:

Example description
-------------------

The example problem consists of two lakes and a stream network. Figure
`31 <#fig:ex-gwt-prudic2004t2-bcmap>`__ shows the configuration of the
hypothetical (but realistic) problem that was used by (Prudic, Konikow,
and Banta 2004) to demonstrate integration of the SFR1 Package with the
LAK3 Package and the MODFLOW-2000 GWT Process.

.. figure:: _images/ex-gwt-prudic2004t2-bcmap.png
   :name: fig:ex-gwt-prudic2004t2-bcmap

   Model grid, boundary conditions and locations of lakes and streams
   used for the stream-lake interaction with solute transport problem.
   From (Prudic, Konikow, and Banta 2004).

As described in (Prudic, Konikow, and Banta 2004), the aquifer is
moderately permeable and has homogeneous properties and uniform
thickness (table `19 <#tab:ex-gwt-prudic2004t2-01>`__). The aquifer was
discretized into 8 layers (each 15 ft thick), 36 rows (at equal spacing
of 405.7 ft), and 23 columns (at equal spacing of 403.7 ft). The flow
field is represented as steady state. Uniform recharge was applied at a
rate of 21 in/yr to layer 1 of the model. Two lakes are located within
the model domain. Lake 1 has inflow from stream segment 1 and has
outflow to stream segment 2. Lake 2 has no stream inflows or outflows.
Both streams are in contact with the aquifer. Lakes are in horizontal
contact with the aquifer in layer 1 around their edges and in vertical
contact with the aquifer in layer 2.

.. container::
   :name: tab:ex-gwt-prudic2004t2-01

   .. table:: Model Parameters

      ===================================================== =========
      **Parameter**                                         **Value**
      ===================================================== =========
      Horizontal hydraulic conductivity (:math:`ft d^{-1}`) 250.0
      Vertical hydraulic conductivity (:math:`ft d^{-1}`)   125.0
      Storage coefficient (unitless)                        0.0
      Aquifer thickness (:math:`ft`)                        120.0
      Porosity of mobile domain (unitless)                  0.30
      Recharge rate (:math:`ft d^{-1}`)                     4.79e-3
      Lakebed leakance (:math:`ft^{-1}`)                    1.0
      Streambed hydraulic conductivity (:math:`ft d^{-1}`)  100.0
      Streambed thickness (:math:`ft`)                      1.0
      Stream width (:math:`ft`)                             5.0
      Manning’s roughness coefficient (unitless)            0.03
      Longitudinal dispersivity (:math:`ft`)                20.0
      Transverse horizontal dispersivity (:math:`ft`)       2.0
      Transverse vertical dispersivity (:math:`ft`)         0.2
      Diffusion coefficient (:math:`ft^2 d^{-1}`)           0.0
      Initial concentration (micrograms per liter)          0.0
      Source concentration (micrograms per liter)           500.0
      Number of layers                                      8
      Number of rows                                        36
      Number of columns                                     23
      Column width (:math:`ft`)                             405.665
      Row width (:math:`ft`)                                403.717
      Layer thickness (:math:`ft`)                          15.0
      Top of the model (:math:`ft`)                         100.0
      Total simulation time (:math:`d`)                     9131.0
      ===================================================== =========

The boundary conditions are illustrated in figure
`31 <#fig:ex-gwt-prudic2004t2-bcmap>`__ and were designed to produce
flow that is generally from north to south. For the numerical model,
constant-head conditions were specified along the northern and southern
edges of the model domain, and no-flow boundaries were set along the
east and west edges of the grid. In MODFLOW 6 lakes can either sit on
top of the model grid, or they can be incised into the model grid as is
done with previous MODFLOW versions. For this example, aquifer cells in
layer 1 (fig. `31 <#fig:ex-gwt-prudic2004t2-bcmap>`__) that share the
same space as the lake are made inactive in the model grid by setting
their IDOMAIN values to zero. The constant-head boundaries were placed
in all 8 layers at the map locations shown in figure
`31 <#fig:ex-gwt-prudic2004t2-bcmap>`__, with two exceptions. The first
exception is related to the last downstream reach of stream segment 4.
The grid cell in which the reach is located and the cell underlying it
in layer 2 are both specified as active aquifer cells rather than as
constant-head cells because a constant head in that cell with the stream
reach would have set the gradient across the streambed. The second
exception is related to the contaminant source (from treated sewage
effluent), which is only introduced into the upper two model layers
(that is, to a total of 8 cells). The four cells in layer 3 that
underlie the contaminant source are specified as active aquifer cells
rather than as constant-head cells to simulate transport beneath the
source. Constant-head elevations were specified as 50.0 ft along the
north boundary, except at the 4 cells in each of the upper two model
layers that represent inflow from the contaminant source, where the
fixed heads were 50.15 ft in the two middle cells and 50.10 ft in the
two outer cells. The constant-head elevations were set to 28.0 ft along
the south boundary.

The stream network consists of four segments and 38 reaches. (In
MODFLOW 6 there is no concept of a segment, however the reaches are
assigned names based on the segment number so that their combined flows
can be compared with the results from MODFLOW-GWT.) The stream depths
are calculated using Manning’s equation assuming a wide rectangular
channel. For all stream reaches, the channel width was assumed constant
at 5.0 ft and the roughness coefficient for the channel was 0.03.
Inflows to segments 1 and 3 were specified (86,400 and 8,640 ft3/d,
respectively). The inflow to stream segment 2 was equal to the outflow
from lake 1, and it was calculated using Manning’s equation assuming a
wide rectangular channel with a depth based on the difference between
the calculated lake stage and the elevation of the top of the streambed
(see (Merritt and Konikow 2000), p. 11). The inflow to segment 4 is
calculated as the sum of outflows from tributary segments 2 and 3. In
MODFLOW 6 the Water Mover Package was used to route the water from the
end of segment 1 into lake 1, and from the southern outlet of lake 1
into stream segment 2.

As reported by (Prudic, Konikow, and Banta 2004) the test problem
focuses on simulation of a boron plume, which results from sewage
effluent. Variables related to the transport simulation are listed in
table `19 <#tab:ex-gwt-prudic2004t2-01>`__. For the purposes of this
test, boron was assumed nonreactive. Molecular diffusion in the aquifer
was assumed a negligible contributor to solute spreading at the scale of
the field problem, so that hydrodynamic dispersion was related solely to
mechanical dispersion, which was computed in MODFLOW 6 as a function of
the specified dispersivity of the medium and the velocity of the flow
field. The initial boron concentrations in the aquifer and in the lakes
were assumed to be zero, and the sewage effluent concentration was
assumed to be 500 micrograms/L. The source concentration in recharge was
assumed to be zero and the concentration in specified inflow to stream
segments 1 and 3 was also zero. The solute-transport model was run for a
period of 25 years.

.. _example-results-6:

Example Results
---------------

The calculated steady-state head distributions in layers 1 and 2 are
shown in figure `32 <#fig:ex-gwt-prudic2004t2-head>`__ for MODFLOW 6.
This figure can be compared to figure 14 in (Prudic, Konikow, and Banta
2004). The heads in layers 3 through 8 are almost identical to the heads
shown for layer 2 (fig. `32 <#fig:ex-gwt-prudic2004t2-head>`__). Flow is
generally from north to south and predominantly horizontal. Because of
recharge at the water table, however, there is a slight vertically
downward flow in most areas. The lakes and streams exert a strong
influence on the location and magnitude of vertical flow. In layer 2,
the good hydraulic connection with the lakebed results in an almost flat
horizontal hydraulic gradient in head beneath the lakes (fig.
`32 <#fig:ex-gwt-prudic2004t2-head>`__). In general the simulated water
table and aquifer heads simulated by MODFLOW 6 are in good agreement
with those simulated by MODFLOW-2000.

.. figure:: _images/ex-gwt-prudic2004t2-head.png
   :name: fig:ex-gwt-prudic2004t2-head

   Contours of head simulated by the MODFLOW 6 GWF Model for the
   stream-lake interaction example. Contours are for (A) layer 1, and
   (B) layer 2. This figure can be compared with figure 14 in (Prudic,
   Konikow, and Banta 2004), which shows head contours as simulated by
   MODFLOW-2000.

The MODFLOW 6 model calculated steady-state stage in lake 1 was 45.07 ft
(compared to 44.97 ft in MODFLOW-GWT) and in lake 2 was 37.15 ft
(compared to 37.14 ft in MODFLOW-2000). Stream segment 1 was mostly a
gaining stream (leakage across streambed was from ground water), and
stream segment 2 was losing such that outflow from the last reach in
segment 2 was only half of the inflow from lake 1. Stream segment 3 was
mostly gaining and outflow from this segment was only from groundwater
leakage. Lastly, stream segment 4 was a losing stream and leakage was
from the stream to the aquifer in every reach. Simulated flows between
MODFLOW 6 and MODFLOW-2000 are in good qualitative agreement, however
there are differences in individual flows, which can be attributed to
slight differences in the way MODFLOW 6 and MODFLOW-2000 simulate lakes,
streams, and groundwater flow.

For the MODFLOW 6 simulation, the 25-year period was divided into 300
time steps (compared with 1229 time steps used for the MODFLOW-GWT
simulation). Advective groundwater flow was solved using the
second-order implicit Total Variation Diminishing (TVD) scheme.
Dispersion was solved using the XT3D approach, which was originally
designed to represent full three-dimensional anisotropic groundwater
flow (Provost, Langevin, and Hughes 2017). Transport through the surface
water system was solved using the Lake Transport (LKT) Package, the
Streamflow Transport (SFT) Package, and the Mover Transport (MVT)
Package, which transfers solute between the lake and stream according to
simulated flows. Additional detail on the transport parameters for the
MODFLOW-GWT simulation are described in (Prudic, Konikow, and Banta
2004) and in (Merritt and Konikow 2000). The calculated concentration in
lake 1 and in the outflow from the last reach in stream segments 2, 3,
and 4 during the 25-year simulation period are shown in figure
`33 <#fig:ex-gwt-prudic2004t2-cvt>`__ for both the MODFLOW 6 simulation
and the MODFLOW-GWT simulation. Note that because stream segment 2 lost
flow to ground water in all reaches and its only source was inflow from
lake 1, solute concentration in the outflow from the last reach in
stream segment 2 was equal to the concentration in the discharge from
lake 1 (fig. `33 <#fig:ex-gwt-prudic2004t2-cvt>`__). The leading edge of
the plume reaches the upstream edge of lake 1 after about 4 years, at
which time the concentration in the lake begins to increase rapidly.
After about 22 years, the part of the plume close to the source and near
the lake has stabilized and the concentration in lake 1 reaches an
equilibrium concentration of 37.2 micrograms/L as simulated by MODFLOW 6
and 37.4 micrograms/L as simulated by MODFLOW-GWT. Although there are
differences in the surface water concentrations simulated by MODFLOW 6
and MODFLOW-GWT, the general pattern and behavior is quite similar.
Differences between the two models are generally attributed to slight
differences in simulated flows as well as slight differences in how
solute transport is represented. For example, MODFLOW-GWT uses the
method-of-characteristics to simulate advective flow, whereas MODFLOW 6
uses an implicit TVD approach. The methods-of-characteristics approach
implemented in MODFLOW-GWT is exceptional for reducing numerical
dispersion, whereas the second-order TVD approach implemented in
MODFLOW 6 is relatively fast and efficient, but it has more numerical
dispersion than MODFLOW-GWT.

.. figure:: _images/ex-gwt-prudic2004t2-cvt.png
   :name: fig:ex-gwt-prudic2004t2-cvt

   Concentration versus time simulated by MODFLOW 6 for the stream-lake
   interaction example. Solid lines represent the MODFLOW 6 simulated
   change in boron concentration in lake 1 and at the end of stream
   segments 2, 3, and 4. Dashed lines are results from the MODFLOW-GWT
   simulation (Prudic, Konikow, and Banta 2004).

As the lake concentration increases, it in turn acts as a source of
contamination to the aquifer in the areas where the lake is a source of
water to the aquifer. Although the lake significantly dilutes the
contaminants that enter it from the aquifer, the lake and the stream
segments downstream from it in effect provide a short circuit for the
relatively fast transmission of low levels of the contaminant. This is
evident in figure `34 <#fig:ex-gwt-prudic2004t2-conc>`__, which shows
the computed solute distributions in layers 1, 3, 5, and 8 after 25
years. The low-concentration part of the plume emanating from the
downgradient side of lake 1 has advanced farther, and is wider, than the
main plume that emanated directly from the source at the north edge of
the model. The influence of groundwater discharge to stream segment 3 is
most apparent in the concentration pattern shown in layer 1 of figure
`34 <#fig:ex-gwt-prudic2004t2-conc>`__ (that is, the 25 microgram/L
contour). Comparison of concentration levels at different depths in the
system indicates that in the southern part of the area, concentrations
generally increase with depth. In contrast, in the northern part
downgradient from the source, the highest concentrations occur in layer
3. These various patterns result from the dilution effect of recharge of
uncontaminated water at the water table coupled with the consequent
downward component of flow, which causes the solute to move slowly
downward as it migrates to the south. Solute concentrations simulated by
MODFLOW 6 and shown in figure `34 <#fig:ex-gwt-prudic2004t2-conc>`__ are
generally in good agreement with those simulated by MODFLOW-GWT (figure
17 in (Prudic, Konikow, and Banta 2004)). There are some differences,
which can be attributed to the slightly different flow field and the
differences in solute transport solution schemes.

.. figure:: _images/ex-gwt-prudic2004t2-conc.png
   :name: fig:ex-gwt-prudic2004t2-conc

   Contours of concentration simulated by the MODFLOW 6 GWT Model for
   the stream-lake interaction example. Contours are shown for (A) layer
   1, (B) layer 3, (C) layer 5, and (D) layer 8. This figure can be
   compared with figure 17 in (Prudic, Konikow, and Banta 2004), which
   shows an equivalent plot using model results simulated by
   MODFLOW-GWT.

.. container:: references hanging-indent
   :name: refs

   .. container::
      :name: ref-bakker2016

      Bakker, Mark, Vincent Post, Christian D Langevin, Joseph D Hughes,
      Jeremy T White, Jeff J Starn, and Michael N Fienen. 2016.
      “Scripting Modflow Model Developmentusing Python and Flopy.”
      *Groundwater* 54 (5): 733–39. https://doi.org/10.1111/gwat.12413.

   .. container::
      :name: ref-elkadi1988

      El‐Kadi, Aly L. 1988. “Applying the Usgs Mass‐transport Model
      (Moc) to Remedial Actions by Recovery Wells.” *Groundwater* 26
      (3): 281–88.
      https://ngwa.onlinelibrary.wiley.com/doi/abs/10.1111/j.1745-6584.1988.tb00391.x.

   .. container::
      :name: ref-gelhar1971

      Gelhar, Lynn W, and M A Collins. 1971. “General Analysis of
      Longitudinal Dispersion in Nonuniform Flow.” *Water Resources
      Research* 7 (6): 1511–21. https://doi.org/0.1029/WR007i006p01511.

   .. container::
      :name: ref-modflow2005

      Harbaugh, Arlen W. 2005. *MODFLOW-2005, the U.s. Geological Survey
      Modular Ground-Water Model—the Ground-Water Flow Process*. U.S.
      Geological Survey Techniques and Methods, book 6, chap. A16,
      variously paged. https://pubs.usgs.gov/tm/2005/tm6A16/.

   .. container::
      :name: ref-harbaugh1996user

      Harbaugh, Arlen W, and Michael G McDonald. 1996. *User’s
      Documentation for Modflow-96, an Update to the U.s. Geological
      Survey Modular Finite-Difference Ground-Water Flow Model*. U.S.
      Geological Survey Techniques of Water-Resources Investigations,
      book 6, chap. A1, 586 p.

   .. container::
      :name: ref-hunt1978

      Hunt, Bruce. 1978. “Dispersive Sources in Uniform Ground-Water
      Flow.” *Journal of the Hydraulics Division* 104 (ASCE 13467
      Proceeding). https://trid.trb.org/view/72451.

   .. container::
      :name: ref-konikow1996three

      Konikow, Leonard F, Daniel J Goode, and George Z Hornberger. 1996.
      *A Three-Dimensional Method-of-Characteristics Solute-Transport
      Model (Moc3d)*. U.S. Geological Survey Water-Resources
      Investigations Report 96–4267, 87 p.
      https://pubs.er.usgs.gov/publication/wri964267.

   .. container::
      :name: ref-modflow6software

      Langevin, Christian D, Joseph D Hughes, Alden M Provost, Edward R
      Banta, Richard G Niswonger, and Sorab Panday. 2017. *MODFLOW 6,
      the U.s. Geological Survey Modular Hydrologic Model*. U.S.
      Geological Survey Computer software Release.
      https://doi.org/10.5066/F76Q1VQV.

   .. container::
      :name: ref-modflow88

      McDonald, Michael G, and Arlen W Harbaugh. 1988. *A Modular
      Three-Dimensional Finite-Difference Ground-Water Flow Model*. U.S.
      Geological Survey Techniques of Water-Resources Investigations,
      book 6, chap. A1, 586 p.
      https://pubs.er.usgs.gov/publication/twri06A1.

   .. container::
      :name: ref-mehl2013

      Mehl, Steffen W, and Mary C Hill. 2013. *MODFLOW–Lgr Documentation
      of Ghost Node Local Grid Refinement (Lgr2) for Multiple Areas and
      the Boundary Flow and Head (Bfh2) Package*. U.S. Geological Survey
      Techniques and Methods, book 6, chap. A44, 43 p.
      https://doi.org/10.3133/tm6A44.

   .. container::
      :name: ref-modflowlak3pack

      Merritt, Michael L, and Leonard F Konikow. 2000. *Documentation of
      a Computer Program to Simulate Lake-Aquifer Interaction Using the
      Modflow Ground-Water Flow Model and the Moc3d Solute-Transport
      Model*. U.S. Geological Survey Water-Resources Investigations
      Report 00–4167, 146 p.
      https://pubs.er.usgs.gov/publication/wri004167.

   .. container::
      :name: ref-moench1981

      Moench, AF, and A Ogata. 1981. “A Numerical Inversion of the
      Laplace Transform Solution to Radial Dispersion in a Porous
      Medium.” *Water Resources Research* 17 (1): 250–52.
      https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/WR017i001p00250.

   .. container::
      :name: ref-modflow6xt3d

      Provost, Alden M, Christian D Langevin, and Joseph D Hughes. 2017.
      *Documentation for the “Xt3d” Option in the Node Property Flow
      (Npf) Package of Modflow 6*. U.S. Geological Survey Techniques and
      Methods, book 6, chap. A56, 46 p. https://doi.org/10.3133/tm6A56.

   .. container::
      :name: ref-modflowsfr1pack

      Prudic, David E., Leonard F. Konikow, and Edward R. Banta. 2004.
      *A New Streamflow-Routing (Sfr1) Package to Simulate
      Stream-Aquifer Interaction with Modflow-2000*. U.S. Geological
      Survey Open File Report 2004–1042, 104 p.
      https://pubs.er.usgs.gov/publication/ofr20041042.

   .. container::
      :name: ref-vanGenuchtenAlves1982

      Van Genuchten, Martinus Theodorus, and W J Alves. 1982.
      *Analytical Solutions of the One-Dimensional Convective-Dispersive
      Solute Transport Equation*. U.S. Department of Agriculture
      Technical Bulletin No. 1661.
      https://naldc.nal.usda.gov/download/CAT82780278/PDF.

   .. container::
      :name: ref-wexler1992

      Wexler, E. J. 1992. *Analytical Solutions for One-, Two-, and
      Three-Dimensional Solute Transport in Ground-Water Systems with
      Uniform Flow*. U.S. Geological Survey Techniques of
      Water-Resources Investigations, Book 3, Chapter B7, 190 p.

   .. container::
      :name: ref-wilson1978

      Wilson, John L, and Paul J Miller. 1978. “Two-Dimensional Plume in
      Uniform Ground-Water Flow.” *Journal of the Hydraulics Division*
      104 (4): 503–14.
      https://cedb.asce.org/CEDBsearch/record.jsp?dockey=0007986.

   .. container::
      :name: ref-zheng1993

      Zheng, C. 1993. “Extension of the Method of Characteristics for
      Simulation of Solute Transport in Three Dimensions.” *Groundwater*
      31 (3): 456–65.
      https://ngwa.onlinelibrary.wiley.com/doi/abs/10.1111/j.1745-6584.1993.tb01848.x.

   .. container::
      :name: ref-zheng1990mt3d

      Zheng, Chunmiao. 1990. *MT3D, a Modular Three-Dimensional
      Transport Model for Simulation of Advection, Dispersion and
      Chemical Reactions of Contaminants in Groundwater Systems*. Report
      to the U.S. Environmental Protection Agency, 170 p.

   .. container::
      :name: ref-zheng2010mt3dmsv5.3

      ———. 2010. *MT3DMS V5.3—a Modular Three-Dimensional Multi-Species
      Transport Model for Simulation of Advection, Dispersion and
      Chemical Reactions of Contaminants in Groundwater Systems;
      Supplemental User’s Guide*. Technical Report, Department of
      Geological Sciences, The University of Alabama, 51 p.

   .. container::
      :name: ref-zheng1999mt3dms

      Zheng, Chunmiao, and P Patrick Wang. 1999. *MT3DMS—a Modular
      Three-Dimensional Multi-Species Transport Model for Simulation of
      Advection, Dispersion and Chemical Reactions of Contaminants in
      Groundwater Systems; Documentation and User’s Guide*. Contract
      report SERDP–99–1: Vicksburg, Miss., U.S. Army Engineer Research
      and Development Center, 169 p.
