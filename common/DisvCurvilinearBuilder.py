import numpy as np

from DisvPropertyContainer import DisvPropertyContainer

__all__ = ["DisvCurvilinearBuilder"]


class DisvCurvilinearBuilder(DisvPropertyContainer):
    """
    A class for generating a curvilinear MODFLOW 6 DISV grid. A curvilinear
    grid is similar to a radial grid, composed of radial bands, but includes
    ncol discretization within a radial band and does not have to form an
    entire circle (such as, a discretized wedge).

    This class inherits from the `DisvPropertyContainer` class and provides
    methods to generate a curvilinear grid using radial and angular parameters.

    All indices are zero-based, but translated to one-base for the figures and
    by flopy for use with MODFLOW 6. Angles are in degrees, with ``0`` being in
    the positive x-axis direction and ``90`` in the positive y-axis direction.

    If no arguments are provided then an empty object is returned.

    Parameters
    ----------
    nlay : int
        Number of layers
    radii : array_like
        List of radial distances that describe the radial bands.
        The first radius is the innermost radius, and then the rest are the
        outer radius of each radial band. Note that the number of radial bands
        is equal to ``len(radii) - 1``.
    angle_start : float
        Starting angle in degrees for the curvilinear grid.
    angle_stop : float
        Stopping angle in degrees for the curvilinear grid.
    angle_step : float
        Column discretization of each radial band.
        If positive, then represents the angle step in degrees for each column
        in a radial band. That is, the number of columns (`ncol`) is:
           ``ncol = (angle_stop - angle_start)/angle_step``
        If negative, then the absolute value is the number of columns (ncol).
    surface_elevation : float or array_like
        Surface elevation for the top layer. Can either be a single float
        for the entire `top`, or array_like of length `nradial`, or
        array_like of length `ncpl`.
    layer_thickness : float or array_like
        Thickness of each layer. Can either be a single float
        for model cells, or array_like of length `nlay`, or
        array_like of length `ncpl`.
    single_center_cell : bool, default=False
        If True, include a single center cell. If true, then innermost `radii`
        must be **zero**. That is, the innermost, radial band has ``ncol=1``.
    origin_x : float, default=0.0
        X-coordinate reference point for the `radii` distance.
    origin_y : float, default=0.0
        Y-coordinate reference point for the `radii` distance.

    Attributes
    ----------
    nradial : int
        Number of radial bands in the grid.
    ncol : int
        Number of columns in each radial band.
    inner_vertex_count : int
        Number of vertices in the innermost radial band.
    single_center_cell : bool
        Whether a single center cell is included.
    full_circle : bool
        Whether the grid spans a full circle. That is,
         full_circle = `angle_start`==`angle_stop`==``0``).
    radii : numpy.ndarray
        Array of radial distances from (origin_x, origin_y) for each radial
        band. The first value is the innermost radius and the remaining are
        each radial bands outer radius.
    angle_start : float
        Starting angle in degrees for the curvilinear grid.
    angle_stop : float
        Stopping angle in degrees for the curvilinear grid.
    angle_step : float
        Angle step in degrees for each column in a radial band.
    angle_span : float
        Span of the angle range in degrees for the curvilinear grid.

    Methods
    -------
    get_disv_kwargs()
        Get the keyword arguments for creating a MODFLOW-6 DISV package.
    plot_grid(...)
        Plot the model grid from `vertices` and `cell2d` attributes.
    get_cellid(rad, col, col_check=True)
        Get the cellid given the radial and column indices.
    get_rad_col(cellid)
        Get the radial and column indices given the cellid.
    get_vertices(rad, col)
        Get the vertex indices for a cell given the radial and column indices.
    calc_curvilinear_ncol(angle_start, angle_stop, angle_step)
        Calculate the number of columns in the curvilinear grid based on
        the given angle parameters. It will adjust `angle_step` to ensure
        that the number of columns is an integer value.
    iter_rad_col()
        Iterate through the radial band columns, then bands.
    iter_radial_cellid(rad)
        Iterate through the cellid within a radial band.
    iter_column_cellid(col)
        Iterate through the cellid along a column across all radial bands.
    """

    nradial: int
    ncol: int
    inner_vertex_count: int
    single_center_cell: bool
    full_circle: bool
    radii: np.ndarray
    angle_start: float
    angle_stop: float
    angle_step: float
    angle_span: float

    def __init__(
        self,
        nlay=-1,
        radii=np.array((0.0, 1.0)),
        angle_start=0.0,
        angle_stop=90.0,
        angle_step=-1,
        surface_elevation=100.0,
        layer_thickness=100.0,
        single_center_cell=False,
        origin_x=0.0,
        origin_y=0.0,
    ):
        if nlay is None or nlay < 1:
            self._init_empty()
            return

        if angle_start < 0.0:
            angle_start += 360.0
        if angle_stop < 0.0:
            angle_stop += 360.0
        if abs(angle_step) < 1.0e-30:
            raise RuntimeError(
                "DisvCurvilinearBuilder: angle_step is near zero"
            )

        angle_span = self._get_angle_span(angle_start, angle_stop)

        ncol, angle_step = self.calc_curvilinear_ncol(
            angle_start, angle_stop, angle_step
        )

        if angle_step > 90.0:
            angle_step = 90.0
            ncol, angle_step = self.calc_curvilinear_ncol(
                angle_start, angle_stop, angle_step
            )

        if angle_span < angle_step:
            raise RuntimeError(
                "DisvCurvilinearBuilder: angle_step is greater than "
                "the total angel, that is:\n"
                "angle_step > |angle_stop - angle_start|\n"
                f"{angle_step} > {angle_span}"
            )

        try:
            nradial = len(radii) - 1
        except TypeError:
            raise RuntimeError(
                "DisvCurvilinearBuilder: radii must be list-like type"
            )

        if nradial < 1:
            raise RuntimeError(
                "DisvCurvilinearBuilder: len(radii) must be greater than 1"
            )

        if single_center_cell and radii[0] > 1.0e-100:
            raise RuntimeError(
                "DisvCurvilinearBuilder: single_center_cell=True must "
                "have the first radii be zero, that is: radii[0] = 0.0\n"
                f"Input received radii[0]={radii[0]}"
            )

        full_circle = 359.999 < angle_span
        nver = ncol if full_circle else ncol + 1

        ncpl = ncol * nradial  # Nodes per layer
        if single_center_cell:
            ncpl = (ncol * nradial) - ncol + 1

        self.radii = np.array(radii, dtype=np.float64)
        self.nradial = nradial
        self.ncol = ncol

        self.single_center_cell = single_center_cell
        self.full_circle = full_circle

        self.angle_start = angle_start
        self.angle_stop = angle_stop
        self.angle_step = angle_step
        self.angle_span = angle_span

        cls_name = "DisvCurvilinearBuilder"
        top = self._get_array(cls_name, surface_elevation, ncpl, nradial)
        thick = self._get_array(cls_name, layer_thickness, nlay, ncpl * nlay)

        if top.size == nradial and nradial != ncpl:
            tmp = []
            for it, rad in top:
                if it == 0 and single_center_cell:
                    tmp.append(rad)
                else:
                    tmp += ncol * [rad]
            top = np.array(tmp)
            del tmp

        bot = []

        if thick.size == nlay:
            for lay in range(nlay):
                bot.append(top - thick[: lay + 1].sum())
        else:
            st = 0
            sp = ncpl
            bt = top.copy()
            for lay in range(nlay):
                bt -= thick[st:sp]
                st, sp = sp, sp + ncpl
                bot.append(bt)

        if single_center_cell and full_circle:
            # Full, filled circle - No vertex at center
            inner_vertex_count = 0
        elif self.radii[0] < 1.0e-100:
            # Single point at circle center
            inner_vertex_count = 1
        else:
            # Innermost vertices are the same as outer bands
            inner_vertex_count = nver

        self.inner_vertex_count = inner_vertex_count

        # Build the grid

        vertices = []
        iv = 0
        stp = np.radians(angle_step)  # angle step in radians

        # Setup center vertex
        if inner_vertex_count == 1:
            vertices.append([iv, 0.0, 0.0])  # Single vertex at center
            iv += 1

        # Setup vertices
        st = 0 if inner_vertex_count > 1 else 1
        for rad in self.radii[st:]:
            ang = np.radians(angle_start)  # angle start in radians
            for it in range(nver):
                xv = rad * np.cos(ang)
                yv = rad * np.sin(ang)
                vertices.append([iv, xv, yv])
                iv += 1
                ang += stp

        # cell2d: [icell2d, xc, yc, ncvert, icvert]
        cell2d = []
        ic = 0
        for rad in range(nradial):
            single_cell_rad0 = self.single_center_cell and rad == 0
            for col in range(ncol):
                icvert = self.get_vertices(rad, col)
                # xc, yc = get_cell_center(rad, col)
                if single_cell_rad0:
                    xc, yc = 0.0, 0.0
                else:
                    xc, yc = self.get_centroid(icvert, vertices)
                cell2d.append([ic, xc, yc, len(icvert), *icvert])
                ic += 1
                if single_cell_rad0:
                    break

        super().__init__(nlay, vertices, cell2d, top, bot, origin_x, origin_y)

    def __repr__(self):
        return super().__repr__("DisvCurvilinearBuilder")

    def _init_empty(self):
        super()._init_empty()
        nul = np.array([])
        self.nradial = 0
        self.ncol = 0
        self.inner_vertex_count = 0
        self.single_center_cell = False
        self.full_circle = False
        self.radii = nul
        self.angle_start = 0
        self.angle_stop = 0
        self.angle_step = 0
        self.angle_span = 0

    def property_copy_to(self, DisvCurvilinearBuilderType):
        if isinstance(DisvCurvilinearBuilderType, DisvCurvilinearBuilder):
            super().property_copy_to(DisvCurvilinearBuilderType)
            DisvCurvilinearBuilderType.nradial = self.nradial
            DisvCurvilinearBuilderType.ncol = self.ncol
            DisvCurvilinearBuilderType.full_circle = self.full_circle
            DisvCurvilinearBuilderType.radii = self.radii
            DisvCurvilinearBuilderType.angle_start = self.angle_start
            DisvCurvilinearBuilderType.angle_stop = self.angle_stop
            DisvCurvilinearBuilderType.angle_step = self.angle_step
            DisvCurvilinearBuilderType.angle_span = self.angle_span
            DisvCurvilinearBuilderType.inner_vertex_count = (
                self.inner_vertex_count
            )
            DisvCurvilinearBuilderType.single_center_cell = (
                self.single_center_cell
            )
        else:
            raise RuntimeError(
                "DisvCurvilinearBuilder.property_copy_to "
                "can only copy to objects that inherit "
                "properties from DisvCurvilinearBuilder"
            )

    def copy(self):
        cp = DisvCurvilinearBuilder()
        self.property_copy_to(cp)
        return cp

    def get_cellid(self, rad, col, col_check=True):
        """
        Get the cellid given the radial and column indices.

        Parameters
        ----------
        rad : int
            Radial index.
        col : int
            Column index.
        col_check : bool, default=True
            If True, than a RuntimeError error is raised for single_center_cell
            grids with ``rad==0`` and ``col>0``. Otherwise, assumes ``col=0``.

        Returns
        -------
        int
            cellid index
        """
        ncol = self.ncol
        if self.single_center_cell:
            # Have to account for only one cell at the center
            if rad == 0 and col > 0:
                if col_check:
                    raise RuntimeError(
                        "DisvCurvilinearBuilder: Bad rad and col given"
                    )
                return 0
            # if rad == 0, then first cell and pos =  0
            # else account for inner cell, plus each ncol band
            pos = 1 + ncol * (rad - 1) + col if rad > 0 else 0
        else:
            pos = rad * ncol + col

        return pos

    def get_rad_col(self, cellid):
        """
        Get the radial and column indices given the cellid.

        Parameters
        ----------
        cellid : int
            cellid index

        Returns
        -------
        (int, int)
            Radial index, Column index
        """
        ncol = self.ncol

        if cellid < 1:
            rad, col = 0, 0
        elif self.single_center_cell:
            cellid -= 1  # drop out first radial band (single cell)
            rad = cellid // ncol + 1
            col = cellid - ncol * (rad - 1)
        else:
            rad = cellid // ncol
            col = cellid - ncol * rad

        return rad, col

    def get_vertices(self, rad, col):
        """
        Get the vertex indices for a cell given the radial and column indices.

        Parameters
        ----------
        rad : int
            Radial index.
        col : int
            Column index.

        Returns
        -------
        list[int]
            List of vertex indices that define the cell at (rad, col).
        """
        ivc = self.inner_vertex_count
        full_circle = self.full_circle
        ncol = self.ncol
        nver = ncol if full_circle else ncol + 1

        if rad == 0:  # Case with no center point or single center point
            if self.single_center_cell:
                return [iv for iv in range(nver + ivc)][::-1]
            elif ivc == 1:  # Single center point
                if full_circle and col == ncol - 1:
                    return [1, col + 1, 0]  # [col+2-nver, col+1, 0]
                return [col + 2, col + 1, 0]
            elif full_circle and col == ncol - 1:
                return [col + 1, nver + col, col, col + 1 - nver]
            else:  # Normal inner band
                return [nver + col + 1, nver + col, col, col + 1]

        n = (rad - 1) * nver + ivc

        if full_circle and col == ncol - 1:
            return [n + col + 1, n + nver + col, n + col, n + col + 1 - nver]

        return [n + nver + col + 1, n + nver + col, n + col, n + col + 1]

    def iter_rad_col(self):
        """Generator that iterates through the radial band columns, then bands.

        Yields
        -------
        (int, int)
            radial band index, column index
        """
        for cellid in range(self.ncpl):
            yield self.get_rad_col(cellid)

    def iter_radial_cellid(self, rad):
        """Generator that iterates through the cellid within a radial band.

        Parameters
        ----------
        rad : int
            Radial index.

        Yields
        -------
        int
            cellid index
        """
        st = self.get_cellid(rad, 0)
        if self.single_center_cell and rad == 0:
            return iter([st])
        sp = self.get_cellid(rad, self.ncol - 1) + 1
        return iter(range(st, sp))

    def iter_column_cellid(self, col):
        """Generator that iterates through the cellid along a column across
        all radial bands.

        Parameters
        ----------
        col : int
            Column index.

        Yields
        -------
        int
            cellid index
        """
        rad = 0
        while rad < self.nradial:
            yield self.get_cellid(rad, col)
            rad += 1

    def iter_columns(self, rad):
        """Generator that iterates through the columns within a radial band.

        Parameters
        ----------
        rad : int
            Radial index.

        Yields
        -------
        int
            column index
        """
        if self.single_center_cell and rad == 0:
            return iter([0])
        return iter(range(0, self.ncol))

    @staticmethod
    def _get_angle_span(angle_start, angle_stop):
        # assumes angles are between 0 and 360
        if abs(angle_stop - angle_start) < 0.001:  # angle_stop == angle_start
            return 360.0
        if angle_start < angle_stop:
            return angle_stop - angle_start
        return 360.0 - angle_start + angle_stop

    @staticmethod
    def calc_curvilinear_ncol(angle_start, angle_stop, angle_step):
        """
        Calculate the number of columns in the curvilinear grid based on
        the given angle parameters. It will adjust `angle_step` to ensure
        that the number of columns is an integer value.

        Parameters
        ----------
        angle_start : float
            Starting angle in degrees for the curvilinear grid.
        angle_stop : float
            Stopping angle in degrees for the curvilinear grid.
        angle_step : float
            If positive, then represents the largest angle step in degrees
            for each column in a radial band. It may be reduced to make
            the number of columns be a positive, integer.
            If negative, then the absolute value is the number of columns
            (ncol) and angle_step is calculated based on it.

        Returns
        -------
        (int, float)
            The number of columns in the curvilinear grid and the angle_step
            that can reproduce the exact integer number.
        """
        angle_span = DisvCurvilinearBuilder._get_angle_span(
            angle_start, angle_stop
        )

        if angle_step > 0.0:
            ncol = int(angle_span // angle_step)
            if (angle_span / angle_step) - ncol > 0.1:  # error towards larger
                ncol += 1
        else:
            ncol = int(round(-1 * angle_step))
        angle_step = angle_span / ncol
        return ncol, angle_step


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    test1 = DisvCurvilinearBuilder(
        1,
        [0, 1, 2, 3],
        angle_start=0,
        angle_stop=90,
        angle_step=10,
        single_center_cell=False,
    )

    tst = test1.copy()
    print(test1)
    print(tst)

    def tmp(**kwargs):
        assert "nlay" in kwargs

    tmp(**test1)
    del test1

    def check(
        nlay,
        radii,
        a_st,
        a_sp,
        a_stp,
        single,
        plot_time=0.0,
        title="",
        dpi=150,
    ):
        ob = DisvCurvilinearBuilder(
            nlay,
            radii,
            angle_start=a_st,
            angle_stop=a_sp,
            angle_step=a_stp,
            single_center_cell=single,
        )

        assert nlay == ob.nlay
        assert single == ob.single_center_cell
        nradial = ob.nradial
        ncol = ob.ncol
        for rad in range(nradial):
            for col in range(ncol):
                node = ob.get_cellid(rad, col)
                assert (rad, col) == ob.get_rad_col(node)

                if single and rad == 0:
                    break

        if title == "":
            title = (
                f"({ob.angle_step}°, {ob.angle_span}°) and "
                f"({ncol}, {nradial}) and SingleCenter={single}"
                f""
            )

        if plot_time < 0:
            fig, _ = ob.plot_grid(show=False)
            plt.close(fig)
        else:
            ob.plot_grid(title=title, plot_time=plot_time, dpi=dpi)

    # Non-Plot Checks Angle -------------------------------------------------

    check(3, [0, 10], 0, 90, 15, False, plot_time=-1)
    check(3, [0, 10], 0, 180, 15, False, plot_time=-1)
    check(3, [0, 10], 0, 270, 15, False, plot_time=-1)
    check(3, [0, 10], 0, 360, 15, False, plot_time=-1)
    check(3, [0, 10], 0, 0, 15, False, plot_time=-1)

    check(3, [10, 20], 0, 90, 15, False, plot_time=-1)
    check(3, [10, 20], 0, 180, 15, False, plot_time=-1)
    check(3, [10, 20], 0, 270, 15, False, plot_time=-1)
    check(3, [10, 20], 0, 360, 15, False, plot_time=-1)
    check(3, [10, 20], 0, 0, 15, False, plot_time=-1)

    check(3, [0, 10], 0, 90, 15, True, plot_time=-1)
    check(3, [0, 10], 0, 180, 15, True, plot_time=-1)
    check(3, [0, 10], 0, 270, 15, True, plot_time=-1)
    check(3, [0, 10], 0, 360, 15, True, plot_time=-1)
    check(3, [0, 10], 0, 0, 15, True, plot_time=-1)

    check(3, [0, 10, 20], 0, 90, 15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 180, 15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 270, 15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 360, 15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 0, 15, False, plot_time=-1)

    check(3, [10, 20, 30], 0, 90, 15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 180, 15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 270, 15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 360, 15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 0, 15, False, plot_time=-1)

    check(3, [0, 10, 20], 0, 90, 15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 180, 15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 270, 15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 360, 15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 0, 15, True, plot_time=-1)

    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 90, 15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 180, 15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 270, 15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 360, 15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 0, 15, False, plot_time=-1)

    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 90, 15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 180, 15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 270, 15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 360, 15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 0, 15, False, plot_time=-1)

    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 90, 15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 180, 15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 270, 15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 360, 15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 0, 15, True, plot_time=-1)

    # Non-Plot Checks NCOL -------------------------------------------------

    check(3, [0, 10], 0, 90, -15, False, plot_time=-1)
    check(3, [0, 10], 0, 180, -15, False, plot_time=-1)
    check(3, [0, 10], 0, 270, -15, False, plot_time=-1)
    check(3, [0, 10], 0, 360, -15, False, plot_time=-1)
    check(3, [0, 10], 0, 0, -15, False, plot_time=-1)

    check(3, [10, 20], 0, 90, -15, False, plot_time=-1)
    check(3, [10, 20], 0, 180, -15, False, plot_time=-1)
    check(3, [10, 20], 0, 270, -15, False, plot_time=-1)
    check(3, [10, 20], 0, 360, -15, False, plot_time=-1)
    check(3, [10, 20], 0, 0, -15, False, plot_time=-1)

    check(3, [0, 10], 0, 90, -15, True, plot_time=-1)
    check(3, [0, 10], 0, 180, -15, True, plot_time=-1)
    check(3, [0, 10], 0, 270, -15, True, plot_time=-1)
    check(3, [0, 10], 0, 360, -15, True, plot_time=-1)
    check(3, [0, 10], 0, 0, -15, True, plot_time=-1)

    check(3, [0, 10, 20], 0, 90, -15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 180, -15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 270, -15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 360, -15, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 0, -15, False, plot_time=-1)

    check(3, [10, 20, 30], 0, 90, -15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 180, -15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 270, -15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 360, -15, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 0, -15, False, plot_time=-1)

    check(3, [0, 10, 20], 0, 90, -15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 180, -15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 270, -15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 360, -15, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 0, -15, True, plot_time=-1)

    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 90, -15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 180, -15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 270, -15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 360, -15, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 0, -15, False, plot_time=-1)

    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 90, -15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 180, -15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 270, -15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 360, -15, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 0, -15, False, plot_time=-1)

    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 90, -15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 180, -15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 270, -15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 360, -15, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 0, -15, True, plot_time=-1)

    # Non-Plot Checks Extreme Points ---------------------------------------

    check(3, [0, 10], 0, 90, -1, False, plot_time=-1)
    check(3, [0, 10], 0, 180, -1, False, plot_time=-1)
    check(3, [0, 10], 0, 270, -1, False, plot_time=-1)
    check(3, [0, 10], 0, 360, -1, False, plot_time=-1)
    check(3, [0, 10], 0, 0, -1, False, plot_time=-1)

    check(3, [10, 20], 0, 90, -1, False, plot_time=-1)
    check(3, [10, 20], 0, 180, -1, False, plot_time=-1)
    check(3, [10, 20], 0, 270, -1, False, plot_time=-1)
    check(3, [10, 20], 0, 360, -1, False, plot_time=-1)
    check(3, [10, 20], 0, 0, -1, False, plot_time=-1)

    check(3, [0, 10], 0, 90, -1, True, plot_time=-1)
    check(3, [0, 10], 0, 180, -1, True, plot_time=-1)
    check(3, [0, 10], 0, 270, -1, True, plot_time=-1)
    check(3, [0, 10], 0, 360, -1, True, plot_time=-1)
    check(3, [0, 10], 0, 0, -1, True, plot_time=-1)

    check(3, [0, 10, 20], 0, 90, -1, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 180, -1, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 270, -1, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 360, -1, False, plot_time=-1)
    check(3, [0, 10, 20], 0, 0, -1, False, plot_time=-1)

    check(3, [10, 20, 30], 0, 90, -1, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 180, -1, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 270, -1, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 360, -1, False, plot_time=-1)
    check(3, [10, 20, 30], 0, 0, -1, False, plot_time=-1)

    check(3, [0, 10, 20], 0, 90, -1, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 180, -1, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 270, -1, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 360, -1, True, plot_time=-1)
    check(3, [0, 10, 20], 0, 0, -1, True, plot_time=-1)

    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 90, 5, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 180, 5, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 270, 5, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 360, 5, False, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 0, 5, False, plot_time=-1)

    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 90, 5, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 180, 5, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 270, 5, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 360, 5, False, plot_time=-1)
    check(3, [10, 20, 30, 40, 50, 60, 70], 0, 0, 5, False, plot_time=-1)

    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 90, 5, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 180, 5, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 270, 5, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 360, 5, True, plot_time=-1)
    check(3, [0, 10, 20, 30, 40, 50, 60, 70], 0, 0, 5, True, plot_time=-1)

    # Plot Checking ---------------------------------------------------------

    check(3, [0, 10, 20, 30, 40, 50], 0, 90, 15, False, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 180, 15, False, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 270, 15, False, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 360, 15, False, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 0, 15, False, plot_time=-1)

    check(3, [10, 20, 30, 40, 50], 0, 90, 15, False, plot_time=0.1)
    check(3, [10, 20, 30, 40, 50], 0, 180, 15, False, plot_time=0.1)
    check(3, [10, 20, 30, 40, 50], 0, 270, 15, False, plot_time=0.1)
    check(3, [10, 20, 30, 40, 50], 0, 360, 15, False, plot_time=0.1)
    check(3, [10, 20, 30, 40, 50], 0, 0, 15, False, plot_time=-1)

    check(3, [0, 10, 20, 30, 40, 50], 0, 90, 15, True, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 180, 15, True, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 270, 15, True, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 360, 15, True, plot_time=0.1)
    check(3, [0, 10, 20, 30, 40, 50], 0, 0, 15, True, plot_time=-1)
