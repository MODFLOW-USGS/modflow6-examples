import numpy as np

__all__ = ["disv_curvilinear_builder"]


class disv_curvilinear_builder:
    nlay: int
    nradial: int
    ncol: int
    ncpl: int
    nodes: int
    inner_vertex_count: int
    single_center_cell: bool
    full_circle: bool
    radii: np.ndarray
    angle_start: float
    angle_stop: float
    angle_step: float
    angle_span: float
    disv_kw: dict

    def __init__(
        self,
        nlay,
        radii,
        angle_start,
        angle_stop,
        angle_step,
        surface_elevation=100.0,
        layer_thickness=100.0,
        single_center_cell=False,
    ):
        if angle_start < 0.0:
            angle_start += 360.0
        if angle_stop < 0.0:
            angle_stop += 360.0
        if abs(angle_step) < 1.0e-30:
            raise RuntimeError(
                "disv_curvilinear_builder: angle_step is near zero"
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
                "disv_curvilinear_builder: angle_step is greater than "
                "the total angel, that is:\n"
                "angle_step > |angle_stop - angle_start|\n"
                f"{angle_step} > {angle_span}"
            )

        try:
            nradial = len(radii) - 1
        except TypeError:
            raise RuntimeError(
                "disv_curvilinear_builder: radii must be list-like type"
            )

        if nradial < 1:
            raise RuntimeError(
                "disv_curvilinear_builder: len(radii) must be greater than 1"
            )

        if single_center_cell and radii[0] > 1.0e-100:
            raise RuntimeError(
                "disv_curvilinear_builder: single_center_cell=True must "
                "have the first radii be zero, that is: radii[0] = 0.0\n"
                f"Input received radii[0]={radii[0]}"
            )

        # if single_center_cell and 90.0 < angle_span < 360.0:
        #     raise RuntimeError(
        #         "disv_curvilinear_builder: single_center_cell=True must "
        #         "|angle_stop - angle_start| be either <90 or equal to 0 or 360"
        #     )

        full_circle = 359.999 < angle_span
        nver = ncol if full_circle else ncol + 1

        ncpl = ncol * nradial  # Nodes per layer
        if single_center_cell:
            ncpl = (ncol * nradial) - ncol + 1

        nodes = ncpl * nlay

        self.radii = np.array(radii, dtype=np.float64)
        self.nlay = nlay
        self.nradial = nradial
        self.ncol = ncol
        self.ncpl = ncpl
        self.nodes = nodes

        self.single_center_cell = single_center_cell
        self.full_circle = full_circle

        self.angle_start = angle_start
        self.angle_stop = angle_stop
        self.angle_step = angle_step
        self.angle_span = angle_span

        top = self._get_rad_array(surface_elevation, ncpl)
        thick = self._get_rad_array(layer_thickness, nlay)

        bot = []
        for lay in range(nlay):
            bot.append(top - thick[: lay + 1].sum())

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
                    xc, yc = self._find_Centroid(icvert, vertices)
                cell2d.append([ic, xc, yc, len(icvert), *icvert])
                ic += 1
                if single_cell_rad0:
                    break

        self.disv_kw = {}
        self.disv_kw["nlay"] = nlay
        self.disv_kw["ncpl"] = ncpl
        self.disv_kw["top"] = top
        self.disv_kw["botm"] = bot
        self.disv_kw["nvert"] = len(vertices)
        self.disv_kw["vertices"] = vertices
        self.disv_kw["cell2d"] = cell2d

    def keys(self):
        # Return keys in disv_kw
        return self.disv_kw.keys()

    def __getitem__(self, k):
        if k in self.disv_kw:
            return self.disv_kw[k]
        if hasattr(self, k):
            return getattr(self, k)
        raise KeyError(f"{k}")

    def get_disv_kw(self):
        return self.disv_kw

    def get_node(self, rad, col, col_check=True):
        ncol = self.ncol
        if self.single_center_cell:
            # Have to account for only one cell at the center
            if rad == 0 and col > 0:
                if col_check:
                    raise RuntimeError(
                        "get_curvilinear_node: Bad rad and col given"
                    )
                return 0
            # if rad == 0, then first cell and pos =  0
            # else account for inner cell, plus each ncol band
            pos = 1 + ncol * (rad - 1) + col if rad > 0 else 0
        else:
            pos = rad * ncol + col

        return pos

    def get_rad_col(self, node):
        ncol = self.ncol

        if node < 1:
            rad, col = 0, 0
        elif self.single_center_cell:
            node -= 1  # drop out first radial band (single cell)
            rad = node // ncol + 1
            col = node - ncol * (rad - 1)
        else:
            rad = node // ncol
            col = node - ncol * rad

        return rad, col

    def get_vertices(self, rad, col):
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
        for node in range(self.ncpl):
            yield self.get_rad_col(node)

    def iter_radial_nodes(self, rad):
        st = self.get_node(rad, 0)
        if self.single_center_cell and rad == 0:
            return iter([st])
        sp = self.get_node(rad, self.ncol - 1) + 1
        return iter(range(st, sp))

    # def iter_column_nodes(self, col):
    #     nradial = self.nradial
    #     rad = 0
    #
    #     def column_node_generator():
    #         nonlocal rad, col, self
    #         while rad < self.nradial:
    #             yield self.get_node(rad, col)
    #             rad += 1
    #     return column_node_generator

    def iter_column_nodes(self, col):
        rad = 0
        while rad < self.nradial:
            yield self.get_node(rad, col)
            rad += 1

    def iter_columns(self, rad):
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
        angle_span = disv_curvilinear_builder._get_angle_span(
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

    @staticmethod
    def _find_Centroid(icvert, vertices):
        # From https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
        nv = len(icvert)
        x = []
        y = []
        for iv in icvert:
            x.append(vertices[iv][1])
            y.append(vertices[iv][2])

        if nv < 3:
            raise RuntimeError("_find_Centroid: len(icvert) < 3")

        if nv == 3:  # Triangle
            return sum(x) / 3, sum(y) / 3

        xc, yc = 0.0, 0.0
        signedArea = 0.0
        for i in range(nv - 1):
            x0, y0, x1, y1 = x[i], y[i], x[i + 1], y[i + 1]
            a = x0 * y1 - x1 * y0
            signedArea += a
            xc += (x0 + x1) * a
            yc += (y0 + y1) * a

        x0, y0, x1, y1 = x1, y1, x[0], y[0]

        a = x0 * y1 - x1 * y0
        signedArea += a
        xc += (x0 + x1) * a
        yc += (y0 + y1) * a

        signedArea *= 0.5
        return xc / (6 * signedArea), yc / (6 * signedArea)

    @staticmethod
    def _get_rad_array(var, rep, rep2=None):
        if rep2 is None:
            rep2 = rep
        try:
            dim = len(var)
        except TypeError:
            dim, var = 1, [var]

        if dim != 1 and dim != rep and dim != rep2:
            msg = "disv_curvilinear_builder(var): var must be a scalar "
            msg += f"or have len(var)=={rep}"
            if rep2 != rep:
                msg += f"or have len(var)=={rep2}"
            raise IndexError(msg)

        if dim == 1:
            return np.full(rep, var[0], dtype=np.float64)
        else:
            return np.array(var, dtype=np.float64)

    def plot(self, show=True, title=""):
        import matplotlib.pyplot as plt

        nradial = self.nradial
        ncol = self.ncol
        single_center_cell = self.single_center_cell
        ncpl = self.disv_kw["ncpl"]
        vertices = self.disv_kw["vertices"]
        cell2d = self.disv_kw["cell2d"]

        x = []
        y = []
        for r in vertices:
            x.append(r[1])
            y.append(r[2])
        xx = []
        yy = []
        for r in cell2d[:ncpl]:
            xx.append(r[1])
            yy.append(r[2])

        fig = plt.figure(figsize=(6, 6), dpi=200)
        ax = fig.add_subplot()
        ax.set_aspect("equal", adjustable="box")
        if title != "":
            ax.set_title(title)

        ax.scatter(x, y, color="blue", s=7, marker="o", zorder=1)
        ax.scatter(xx, yy, color="red", s=3, marker="o", zorder=2)

        ic = -1
        for rad in range(nradial):
            for col in range(ncol):
                ic += 1
                ncon = cell2d[ic][3]
                conn = cell2d[ic][4:] + [cell2d[ic][4]]
                for i in range(ncon):
                    n1, n2 = conn[i], conn[i + 1]
                    px = [vertices[n1][1], vertices[n2][1]]
                    py = [vertices[n1][2], vertices[n2][2]]
                    ax.plot(px, py, color="black", zorder=0)
                if single_center_cell and rad == 0:
                    break
        if show:
            fig.show()
        return fig, ax


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    test1 = disv_curvilinear_builder(
        1,
        [0, 1, 2, 3],
        angle_start=0,
        angle_stop=90,
        angle_step=10,
        single_center_cell=False,
    )

    def tmp(**kwargs):
        assert "nlay" in kwargs

    tmp(**test1)
    del test1

    def check(nlay, radii, a_st, a_sp, a_stp, single, plot_time=0.0, title=""):
        ob = disv_curvilinear_builder(
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
                node = ob.get_node(rad, col)
                assert (rad, col) == ob.get_rad_col(node)

                if single and rad == 0:
                    break

        def ColorCycler():
            c = ("black", "green", "red", "grey", "magenta", "cyan", "yellow")
            i, d = -1, len(c)

            def dummy():
                nonlocal i, d
                i += 1
                if i == d:
                    i = 0
                return c[i]

            return dummy

        if title == "":
            title = (
                f"({ob.angle_step}°, {ob.angle_span}°) and "
                f"({ncol}, {nradial}) and SingleCenter={single}"
                f""
            )
        if plot_time < 0:
            fig, _ = ob.plot(False, title)
            plt.close(fig)
            return

        ncpl = ob.disv_kw["ncpl"]
        vertices = ob.disv_kw["vertices"]
        cell2d = ob.disv_kw["cell2d"]

        x = []
        y = []
        for r in vertices:
            x.append(r[1])
            y.append(r[2])
        xx = []
        yy = []
        for r in cell2d[:ncpl]:
            xx.append(r[1])
            yy.append(r[2])

        fig = plt.figure(figsize=(6, 6), dpi=200)
        ax = fig.add_subplot()
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title)

        ax.scatter(x, y, color="blue", s=7, marker="o")
        ax.scatter(xx, yy, color="red", s=3, marker="o")

        if plot_time > 0:
            plt.show(block=False)
            plt.pause(5 * plot_time)

        ic = -1
        ColorCycle = ColorCycler()
        for rad in range(nradial):
            for col in range(ncol):
                color = ColorCycle()
                ic += 1
                ncon = cell2d[ic][3]
                conn = cell2d[ic][4:] + [cell2d[ic][4]]
                for i in range(ncon):
                    n1, n2 = conn[i], conn[i + 1]
                    px = [vertices[n1][1], vertices[n2][1]]
                    py = [vertices[n1][2], vertices[n2][2]]
                    ax.plot(px, py, color=color)
                    if plot_time > 0:
                        plt.draw()
                        plt.pause(plot_time)
                if single and rad == 0:
                    break

        _ = ob.plot(False, title)
        plt.show()

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
