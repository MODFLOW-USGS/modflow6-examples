import copy
from itertools import cycle
from typing import List

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["DisvPropertyContainer"]


class DisvPropertyContainer:

    """
    Dataclass that stores MODFLOW 6 DISV grid information.

    This is a base class that stores DISV **kwargs information used
    by flopy for building a ``flopy.mf6.ModflowGwfdisv`` object.

    All indices are zero-based, but translated to one-base for the figures and
    by flopy for use with MODFLOW 6.

    If no arguments are provided then an empty object is returned.

    Parameters
    ----------
    nlay : int
        Number of layers.
    vertices : list[list[int, float, float]]
        List of vertices structured as
           ``[[iv, xv, vy], ...]``
        where
           ``iv`` is the vertex index,
           ``xv`` is the x-coordinate, and
           ``yv`` is the y-coordinate.
    cell2d : list[list[int, float, float, int, int...]]
        List of MODFLOW 6 cells structured as
           ```[[icell2d, xc, yc, ncvert, icvert], ...]```
        where
           ``icell2d`` is the cell index,
           ``xc`` is the x-coordinate for the cell center,
           ``yc`` is the y-coordinate for the cell center,
           ``ncvert`` is the number of vertices required to define the cell,
           ``icvert`` is a list of vertex indices that define the cell, and
              in clockwise order.
    top : np.ndarray
        Is the top elevation for each cell in the top model layer.
    botm : list[np.ndarray]
        List of bottom elevation by layer for all model cells.
    origin_x : float, default=0.0
        X-coordinate of the origin used as the reference point for other
        vertices. This is used for shift and rotate operations.
    origin_y : float, default=0.0
        X-coordinate of the origin used as the reference point for other
        vertices. This is used for shift and rotate operations.
    rotation : float, default=0.0
        Rotation angle in degrees for the model grid.
    shift_origin : bool, default=True
        If True and `origin_x` or `origin_y` is non-zero, then all vertices are
        shifted from an assumed (0.0, 0.0) origin to the (origin_x, origin_y)
        location.
    rotate_grid, default=True
        If True and `rotation` is non-zero, then all vertices are rotated by
         rotation degrees around (origin_x, origin_y).

    Attributes
    ----------
    nlay : int
        Number of layers.
    ncpl : int
        Number of cells per layer.
    nvert : int
        Number of vertices.
    vertices : list[list]
        List of vertices structured as ``[[iv, xv, vy], ...]``
    cell2d : list[list]
        List of 2D cells structured as ```[[icell2d, xc, yc, ncvert, icvert], ...]```
    top : np.ndarray
        Top elevation for each cell in the top model layer.
    botm : list[np.ndarray]
        List of bottom elevation by layer for all model cells.
    origin_x : float
        X-coordinate reference point used by grid.
    origin_y : float
        Y-coordinate reference point used by grid.
    rotation : float
        Rotation angle of grid about (origin_x, origin_y)

    Methods
    -------
    get_disv_kwargs()
        Get the keyword arguments for creating a MODFLOW-6 DISV package.
    plot_grid(...)
        Plot the model grid from `vertices` and `cell2d` attributes.
    change_origin(new_x_origin, new_y_origin)
        Change the origin of the grid.
    rotate_grid(rotation)
        Rotate the grid.
    get_centroid(icvert, vertices=None)
        Calculate the centroid of a cell given by list of vertices `icvert`.
    copy()
        Create and return a copy of the current object.
    """

    nlay: int
    ncpl: int
    nvert: int
    vertices: List[list]  # [[iv, xv, yv], ...]
    cell2d: List[list]  # [[ic, xc, yc, ncvert, icvert], ...]
    top: np.ndarray
    botm: List[np.ndarray]
    origin_x: float
    origin_y: float
    rotation: float

    def __init__(
        self,
        nlay=-1,
        vertices=None,
        cell2d=None,
        top=None,
        botm=None,
        origin_x=0.0,
        origin_y=0.0,
        rotation=0.0,
        shift_origin=True,
        rotate_grid=True,
    ):
        if nlay is None or nlay < 1:
            self._init_empty()
            return

        self.nlay = nlay
        self.ncpl = len(cell2d)
        self.nvert = len(vertices)

        self.vertices = [] if vertices is None else copy.deepcopy(vertices)
        self.cell2d = [] if cell2d is None else copy.deepcopy(cell2d)
        self.top = np.array([]) if top is None else copy.deepcopy(top)
        self.botm = [] if botm is None else copy.deepcopy(botm)

        self.origin_x, self.origin_y, self.rotation = 0.0, 0.0, 0.0

        if shift_origin:
            if abs(origin_x) > 1.0e-30 or abs(origin_y) > 1.0e-30:
                self.change_origin(origin_x, origin_y)
        elif not shift_origin:
            self.origin_x, self.origin_y = origin_x, origin_y

        if rotate_grid:
            self.rotate_grid(rotation)
        elif not shift_origin:
            self.rotation = rotation

    def get_disv_kwargs(self):
        """
        Get the dict of keyword arguments for creating a MODFLOW-6 DISV
        package using ``flopy.mf6.ModflowGwfdisv``.
        """
        return {
            "nlay": self.nlay,
            "ncpl": self.ncpl,
            "top": self.top,
            "botm": self.botm,
            "nvert": self.nvert,
            "vertices": self.vertices,
            "cell2d": self.cell2d,
        }

    def __repr__(self, cls="DisvPropertyContainer"):
        return (
            f"{cls}(\n\n"
            f"nlay={self.nlay}, ncpl={self.ncpl}, nvert={self.nvert}\n\n"
            f"origin_x={self.origin_x}, origin_y={self.origin_y}, "
            f"rotation={self.rotation}\n\n"
            f"vertices =\n{self._string_repr(self.vertices)}\n\n"
            f"cell2d =\n{self._string_repr(self.cell2d)}\n\n"
            f"top =\n{self.top}\n\n"
            f"botm =\n{self.botm}\n\n)"
        )

    def _init_empty(self):
        self.nlay = 0
        self.ncpl = 0
        self.nvert = 0
        self.vertices = []
        self.cell2d = []
        self.top = np.array([])
        self.botm = []
        self.origin_x = 0.0
        self.origin_y = 0.0
        self.rotation = 0.0

    def change_origin(self, new_x_origin, new_y_origin):
        shift_x_origin = new_x_origin - self.origin_x
        shift_y_origin = new_y_origin - self.origin_y

        self.shift_origin(shift_x_origin, shift_y_origin)

    def shift_origin(self, shift_x_origin, shift_y_origin):
        if abs(shift_x_origin) > 1.0e-30 or abs(shift_y_origin) > 1.0e-30:
            self.origin_x += shift_x_origin
            self.origin_y += shift_y_origin

            for vert in self.vertices:
                vert[1] += shift_x_origin
                vert[2] += shift_y_origin

            for cell in self.cell2d:
                cell[1] += shift_x_origin
                cell[2] += shift_y_origin

    def rotate_grid(self, rotation):
        """Rotate grid around origin_x, origin_y for given angle in degrees.

        References
        ----------
        [1] https://en.wikipedia.org/wiki/Transformation_matrix#Rotation

        """
        #
        if abs(rotation) > 1.0e-30:
            self.rotation += rotation

            sin, cos = np.sin, np.cos
            a = np.radians(rotation)
            x0, y0 = self.origin_x, self.origin_y
            # Steps to shift
            #   0) Get x, y coordinate to be shifted
            #   1) Shift coordinate's reference point to origin
            #   2) Rotate around origin
            #   3) Shift back to original reference point
            for vert in self.vertices:
                _, x, y = vert
                x, y = x - x0, y - y0
                x, y = x * cos(a) - y * sin(a), x * sin(a) + y * cos(a)
                vert[1] = x + x0
                vert[2] = y + y0

            for cell in self.cell2d:
                _, x, y, *_ = cell
                x, y = x - x0, y - y0
                x, y = x * cos(a) - y * sin(a), x * sin(a) + y * cos(a)
                cell[1] = x + x0
                cell[2] = y + y0

    @staticmethod
    def _string_repr(list_list, sep=",\n"):
        dim = len(list_list)
        s = []
        if dim == 0:
            return "[]"
        if dim < 7:
            for elm in list_list:
                s.append(repr(elm))
        else:
            for it in range(3):
                s.append(repr(list_list[it]))
            s.append("...")
            for it in range(-3, 0):
                s.append(repr(list_list[it]))
        return sep.join(s)

    def property_copy_to(self, DisvPropertyContainerType):
        if isinstance(DisvPropertyContainerType, DisvPropertyContainer):
            DisvPropertyContainerType.nlay = self.nlay
            DisvPropertyContainerType.ncpl = self.ncpl
            DisvPropertyContainerType.nvert = self.nvert
            DisvPropertyContainerType.vertices = self.vertices
            DisvPropertyContainerType.cell2d = self.cell2d
            DisvPropertyContainerType.top = self.top
            DisvPropertyContainerType.botm = self.botm
            DisvPropertyContainerType.origin_x = self.origin_x
            DisvPropertyContainerType.origin_y = self.origin_y
            DisvPropertyContainerType.rotation = self.rotation
        else:
            raise RuntimeError(
                "DisvPropertyContainer.property_copy_to "
                "can only copy to objects that inherit "
                "properties from DisvPropertyContainer"
            )

    def copy(self):
        cp = DisvPropertyContainer()
        self.property_copy_to(cp)
        return cp

    def keys(self):
        """
        Get the keys used by ``flopy.mf6.ModflowGwfdisv``.

        This method is only used to provide unpacking support for
        `DisvPropertyContainer` objects and subclasses.

        That is:
        ``flopy.mf6.ModflowGwfdisv(gwf, **DisvPropertyContainer)``

        Returns
        -------
        list
            List of keys used by ``flopy.mf6.ModflowGwfdisv``

        """
        return self.get_disv_kwargs().keys()

    def __getitem__(self, k):
        if hasattr(self, k):
            return getattr(self, k)
        raise KeyError(f"{k}")

    @staticmethod
    def _get_array(cls_name, var, rep, rep2=None):
        if rep2 is None:
            rep2 = rep

        try:
            dim = len(var)
        except TypeError:
            dim, var = 1, [var]

        try:  # Check of array needs to be flattened
            _ = len(var[0])
            tmp = []
            for row in var:
                tmp.extend(row)
            var = tmp
            dim = len(var)
        except TypeError:
            pass

        if dim != 1 and dim != rep and dim != rep2:
            msg = f"{cls_name}(var): var must be a scalar "
            msg += f"or have len(var)=={rep}"
            if rep2 != rep:
                msg += f"or have len(var)=={rep2}"
            raise IndexError(msg)

        if dim == 1:
            return np.full(rep, var[0], dtype=np.float64)
        else:
            return np.array(var, dtype=np.float64)

    def get_centroid(self, icvert, vertices=None):
        """
        Calculate the centroid of a cell for a given set of vertices.

        Parameters
        ----------
        icvert : list[int]
            List of vertex indices for the cell.
        vertices : list[list], optional
            List of vertices that `icvert` references to define the cell.
            If not present, then the `vertices` attribute is used.

        Returns
        -------
        tuple
            A tuple containing the X and Y coordinates of the centroid.

        References
        ----------
        [1] https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
        """
        if vertices is None:
            vertices = self.vertices

        nv = len(icvert)
        x = []
        y = []
        for iv in icvert:
            x.append(vertices[iv][1])
            y.append(vertices[iv][2])

        if nv < 3:
            raise RuntimeError("get_centroid: len(icvert) < 3")

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

    def plot_grid(
        self,
        title="",
        plot_time=0.0,
        show=True,
        figsize=(10, 10),
        dpi=None,
        xlabel="",
        ylabel="",
        cell2d_override=None,
        vertices_override=None,
        ax_override=None,
        cell_dot=True,
        cell_num=True,
        cell_dot_size=7.5,
        cell_dot_color="coral",
        vertex_dot=True,
        vertex_num=True,
        vertex_dot_size=6.5,
        vertex_dot_color="skyblue",
        grid_color="grey",
    ):
        """
        Plot the model grid with optional features.
        All inputs are optional.

        Parameters
        ----------
        title : str, default=""
            Title for the plot.
        plot_time : float, default=0.0
            Time interval for animation (if greater than 0).
        show : bool, default=True
            Whether to display the plot. If false, then plot_time is set to -1
            and function returns figure and axis objects.
        figsize : tuple, default=(10, 10)
            Figure size (width, height) in inches. Default is (10, 10).
        dpi : float, default=None
            Set dpi for Matplotlib figure.
            If set to None, then uses Matplotlib default dpi.
        xlabel : str, default=""
            X-axis label.
        ylabel : str, default=""
            Y-axis label.
        cell2d_override : list[list], optional
            List of ``cell2d`` cells to override the object's cell2d.
            Default is None.
        vertices_override : list[list], optional
            List of vertices to override the object's vertices.
            Default is None.
        ax_override : matplotlib.axes.Axes, optional
            Matplotlib axis object to use for generating plot instead of
            making a new figure and axis objects. If present, then show is
            set to False and plot_time to -1.
        cell_dot : bool, default=True
            Whether to add a filled circle at the cell center locations.
        cell_num : bool, default=True
            Whether to label cells with cell2d index numbers.
        cell_dot_size : bool, default=7.5
            The size, in points, of the filled circles and index numbers
            at the cell center.
        cell_dot_color : str, default="coral"
            The color of the filled circles at the cell center.
        vertex_num : bool, default=True
            Whether to label vertices with the vertex numbers.
        vertex_dot : bool, default=True
            Whether to add a filled circle at the vertex locations.
        vertex_dot_size : bool, default=6.5
            The size, in points, of the filled circles and index numbers
            at the vertex locations.
        vertex_dot_color : str, default="skyblue"
            The color of the filled circles at the vertex locations.
        grid_color : str or tuple[str], default="grey"
            The color of the grid lines.
            If plot_time > 0, then animation cycled through the colors
            for each cell outline.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axis.Axis) or None
            If `show` is False, returns the Figure and Axis objects;
            otherwise, returns None.

        Raises
        ------
        RuntimeError
            If either `cell2d_override` or `vertices_override` is provided
            without the other.

        Notes
        -----
        This method plots the grid using Matplotlib. It can label cells and
        vertices with numbers, show an animation of the plotting process, and
        customize the plot appearance.

        Note that figure size (`figsize`) is in inches and the total pixels is based on
        the `dpi` (dots per inch). For example, figsize=(3, 5) and
        dpi=110, results in a figure resolution of (330, 550).

        Changing `figsize` does not effect the size

        Elements (text, markers, lines) in Matplotlib use
        72 points per inch (ppi) as a basis for translating to dpi.
        Any changes to dpi result in a scaling effect. The default dpi of 100
        results in a line width 100/72 pixels wide.

        Similarly, a line width of 1 point with a dpi set to 72 results in
        a line that is 1 pixel wide. The following then occurs:

        2*72 dpi results in a line width 2 pixels;
        3*72 dpi results in a line width 3 pixels; and
        600 dpi results in a line width 600/72 pixels.

        Conversely, changing `figsize` increases the total pixel count, but
        elements maintain the same dpi. That is the figure wireframe will be
        larger, but the elements will have the same pixel widths. For example,
        a line width of 1 will have a width of 100/72 in for any figure size,
        as long as the dpi is set to 100.
        """

        if cell2d_override is not None and vertices_override is not None:
            cell2d = cell2d_override
            vertices = vertices_override
        elif cell2d_override is not None or vertices_override is not None:
            raise RuntimeError(
                "plot_vertex_grid: if you specify "
                "cell2d_override or vertices_override, "
                "then you must specify both."
            )
        else:
            cell2d = self.cell2d
            vertices = self.vertices

        if ax_override is None:
            fig = plt.figure(figsize=figsize, dpi=dpi)
            ax = fig.add_subplot()
        else:
            show = False
            ax = ax_override

        ax.set_aspect("equal", adjustable="box")

        if not show:
            plot_time = -1.0

        if not isinstance(grid_color, tuple):
            grid_color = (grid_color, )

        ColorCycler = grid_color
        if plot_time > 0.0 and grid_color == ("grey",):
            ColorCycler = ("green", "red", "grey", "magenta", "cyan", "yellow")
        ColorCycle = cycle(ColorCycler)

        xvert = []
        yvert = []
        for r in vertices:
            xvert.append(r[1])
            yvert.append(r[2])
        xcell = []
        ycell = []
        for r in cell2d:
            xcell.append(r[1])
            ycell.append(r[2])

        if title != "":
            ax.set_title(title)
        if xlabel != "":
            ax.set_xlabel(xlabel)
        if ylabel != "":
            ax.set_ylabel(ylabel)

        vert_size = vertex_dot_size
        cell_size = cell_dot_size

        if vertex_dot:
            ax.plot(
                xvert,
                yvert,
                linestyle="None",
                color=vertex_dot_color,
                markersize=vert_size,
                marker="o",
                markeredgewidth=0.0,
                zorder=2,
            )
        if cell_dot:
            ax.plot(
                xcell,
                ycell,
                linestyle="None",
                color=cell_dot_color,
                markersize=cell_size,
                marker="o",
                markeredgewidth=0.0,
                zorder=2,
            )

        if cell_num:
            for ic, xc, yc, *_ in cell2d:
                ax.text(
                    xc,
                    yc,
                    f"{ic + 1}",
                    fontsize=cell_size,
                    color="black",
                    fontfamily="Arial Narrow",
                    fontweight="black",
                    rasterized=False,
                    horizontalalignment="center",
                    verticalalignment="center",
                    zorder=3,
                )

        if vertex_num:
            for iv, xv, yv in vertices:
                ax.text(
                    xv,
                    yv,
                    f"{iv + 1}",
                    fontsize=vert_size,
                    fontweight="black",
                    rasterized=False,
                    color="black",
                    fontfamily="Arial Narrow",
                    horizontalalignment="center",
                    verticalalignment="center",
                    zorder=3,
                )

        if plot_time > 0:
            plt.show(block=False)
            plt.pause(5 * plot_time)

        for ic, xc, yc, ncon, *vert in cell2d:
            color = next(ColorCycle)

            conn = vert + [vert[0]]  # Extra node to complete polygon

            for i in range(ncon):
                n1, n2 = conn[i], conn[i + 1]
                px = [vertices[n1][1], vertices[n2][1]]
                py = [vertices[n1][2], vertices[n2][2]]
                ax.plot(px, py, color=color, zorder=1)

                if plot_time > 0:
                    plt.draw()
                    plt.pause(plot_time)

        if show:
            plt.show()
        elif ax_override is None:
            return fig, ax


if __name__ == "__main__":
    # 2 by 3 example
    vertices = [
        [0, 0.0, 0.0],
        [1, 10.0, 0.0],
        [2, 20.0, 0.0],
        [3, 30.0, 0.0],
        [4, 0.0, 10.0],
        [5, 10.0, 10.0],
        [6, 20.0, 10.0],
        [7, 30.0, 10.0],
        [8, 0.0, 20.0],
        [9, 10.0, 20.0],
        [10, 20.0, 20.0],
        [11, 30.0, 20.0],
    ]

    cell2d = [
        [0, 5.0, 5.0, 4, 4, 5, 1, 0],
        [1, 15.0, 5.0, 4, 5, 6, 2, 1],
        [2, 25.0, 5.0, 4, 6, 7, 3, 2],
        [3, 5.0, 15.0, 4, 8, 9, 5, 4],
        [4, 15.0, 15.0, 4, 9, 10, 6, 5],
        [5, 25.0, 15.0, 4, 10, 11, 7, 6],
    ]

    top = np.array([100.0, 100.0, 100.0, 100.0, 100.0, 100.0])
    botm = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    test = DisvPropertyContainer(
        nlay=1,
        vertices=vertices,
        cell2d=cell2d,
        top=top,
        botm=botm,
        origin_x=10.0,
        origin_y=10.0,
        rotation=30.0,
        shift_origin=True,
        rotate_grid=True,
    )

    fig, ax = test.plot_grid(show=False)

    test_empty = DisvPropertyContainer()
    test_copy = test.copy()

    pass
