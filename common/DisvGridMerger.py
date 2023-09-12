from __future__ import annotations

from math import sqrt
import copy

import matplotlib.path as mpltPath
import numpy as np

from DisvPropertyContainer import DisvPropertyContainer

__all__ = ["DisvGridMerger"]


class DisvGridMerger:
    """
    Class for merging, non-overlapping, MODFLOW 6 DISV grids. The merge is
    made by selecting a connection point and adjusting the (x,y) coordinates
    of one of the grids. The grid connection is made by starting with the first
    grid, called `__main__`, then adjusting the second grid to have __main__
    incorporate both grids. After that, subsequent grids are snapped to the
    __main__ grid to form the final merged grid.

    When a grid is shifted to snap to __main__ (`snap_vertices`), any vertices
    that are in proximity of the __main__ grid are merged (that is, the
    snapped grid drops the overlapping vertices and uses the existing __main__
    ones). Proximately is determined by having an x or y distance less than
    `connect_tolerance`.

    Vertices can also be forced to snap to the __main__ grid with `force_snap`.
    A force snap occurs after the second grid is shifted and snapped to
    __main__. The force snap drops the existing vertex and uses the forced one
    changing the shape of the cell. Note, if any existing vertices are located
    within the new shape of the cell, then they are added to the cell2d vertex
    list.

    Examples
    --------
    >>> # Example snaps two rectangular structured vertex grids.
    >>> # The first grid has 2 rows and 2 columns;
    >>> #    the second grid has 3 rows and 2 columns.
    >>> # DisvStructuredGridBuilder builds a DisvPropertyContainer object
    >>> #    that contains a structured vertex grid (rows and columns).
    >>> from DisvStructuredGridBuilder import DisvStructuredGridBuilder
    >>>
    >>> grid1 = DisvStructuredGridBuilder(nlay=1, nrow=2, ncol=2)
    >>> grid2 = DisvStructuredGridBuilder(nlay=1, nrow=3, ncol=2)
    >>>
    >>> # Optional step to see what vertex point to use
    >>> grid1.plot_grid()  # Plot and view vertex locations
    >>> grid2.plot_grid()  #    to identify connection points.
    >>>
    >>> # Steps to merge grid1 and grid2
    >>> mg = DisvGridMerger()  # init the object
    >>>
    >>> mg.add_grid("grid1", grid1)  # add grid1
    >>> mg.add_grid("grid2", grid2)  # add grid2
    >>>
    >>> # Snap grid1 upper right corner (vertex 3) to grid2 upper left
    >>> #    corner (vertex 1). Note the vertices must be zero-based indexed.
    >>> mg.set_vertex_connection("grid1", "grid2", 3 - 1, 1 - 1)
    >>>
    >>> # Grids do not require any force snapping because overlapping vertices
    >>> #    will be within the connect_tolerance. Otherwise, now would be
    >>> #    when to run the set_force_vertex_connection method.
    >>> # Merge the grids
    >>> mg.merge_grids()
    >>>
    >>> mg.merged.plot_grid()  # plot the merged grid

    Attributes
    ----------
    grids : dict
        A dictionary containing names of individual grids as keys and
        corresponding `DisvPropertyContainer` objects as values.
        The key `__main__` is used to refer to the final merged grid and is not
        allowed as a name for any `DisvPropertyContainer` object.
    merged : DisvPropertyContainer
        A `DisvPropertyContainer` object representing the merged grid.
    snap_vertices : dict
        A dictionary of vertex connections to be snapped during
        the merging process. This attribute is set with the
        ``set_vertex_connection`` method. The key is ``(name1, name2)`` and
        value is ``(vertex1, vertex2)``, where name1 and name2 correspond with
        keys from `grids` and ``vertex1`` and ``vertex2`` are the connection
        vertices for ``name1`` and ``name2``, respectively.
    connect_tolerance : dict
        A dictionary specifying the tolerance distance for vertex snapping.
        After a grid is snapped to __main__ via snap_vertices, any vertices
        that overlap within an x or y length of connect_tolerance are merged.
    snap_order : list
        A list of grid pairs indicating the order in which grids
        will be merged. This is variable is set after running the
        `merge_grids` method.
    force_snap : dict
        A dictionary of vertex connections that must be snapped,
        even if they don't satisfy the tolerance. The key is ``(name1, name2)``
        and value is ``[[v1, ...], [v2, ...]]``, where ``name1`` and ``name2``
        correspond with keys from `grids` and ``v1`` is a list of verties to
        snap from ``name1``, and ``v2`` is a list of vertices to snap to from
        ``name2``. The first ``v1``, corresponds with the first ``v2``,
        and so forth.
    force_snap_drop : dict
        A dictionary specifying which vertex to drop when force snapping. The
        key is ``(name1, name2)`` and value is ``[v_drop, ...]``, where
        ``v_drop`` is 1 to drop the vertex from ``name1``, and 2 to drop
        from ``name2``.
    force_snap_cellid : set
        A set that lists all the merged grid cellids that had one or more
        verties force snapped. This list is important for checking if the
        new vertex list and cell center are correct.
    vert2name : dict
        A dictionary mapping the merged grid's vertex numbers to the
        corresponding grid names and vertex indices. The key is the vertex
        from the new merged grid and the value is
        ``[[name, vertex_old], ...]``, where ``name`` is the original grid name
        and ``vertex_old`` is its correspondnig vertex from name.
    name2vert : dict
        A dictionary mapping grid names and vertex indices to the merged
        grid's vertices. The key is ``(name, vertex_old)``, where ``name`` is
        the original grid name and ``vertex_old`` is its correspondnig vertex
        from name. The value is the merged grid's vertex.
    cell2name : dict
        A dictionary mapping the merged grid's cellid's to the corresponding
        original grid names and cellid's. The key is the merged grid's cellid
        and value is ``(name, cellid_old)``, where ``name`` is the
        original grid name and ``cellid_old`` is its correspondnig cellid from
        name.
    name2cell : dict
        A dictionary mapping grid names and cellid's to the merged grid's
        cellid's. The key is ``(name, cellid_old)``, where ``name`` is the
        original grid name and ``cellid_old`` is its correspondnig cellid from
        name, and value is the merged grid's cellid.

    Notes
    -------
    The following is always true:

    ``cell2name[cell] ==  name2vert[cell2name[cell]]``

    ``name2vert[(name, vertex)] is in vert2name[name2vert[(name, vertex)]]``

    Methods
    -------
    get_disv_kwargs(name="__main__")
        Get the keyword arguments for creating a MODFLOW-6 DISV package for
        a specified grid.
    add_grid(name, grid)
        Add an individual grid to the merger.
    set_vertex_connection(name1, name2, vertex1, vertex2, autosnap_tolerance=1.0e-5)
        Set a vertex connection between two grids for snapping.
    set_force_vertex_connection(name1, name2, vertex1, vertex2, drop_vertex=2)
        Force a vertex connection between two grids.
    merge_grids()
        Merge the specified grids based on the defined vertex connections.
    plot_grid(name="__main__", ...)
        Selects the grid specified by ``name`` and passes the remaining
        kwargs to DisvPropertyContainer.plot_grid(...).
    """

    grids: dict[str, DisvPropertyContainer]
    merged: DisvPropertyContainer
    snap_vertices: dict
    connect_tolerance: dict
    snap_order: list
    force_snap: dict
    force_snap_drop: dict
    force_snap_cellid: set
    vert2name: dict
    name2vert: dict
    cell2name: dict
    name2cell: dict

    def __init__(self):
        self.grids = {}
        self.merged = DisvPropertyContainer()

        self.snap_vertices = {}
        self.connect_tolerance = {}
        self.snap_order = []
        self.force_snap = {}
        self.force_snap_drop = {}
        self.force_snap_cellid = set()

        self.vert2name = {}  # vertex: [[name, vertex], ...]
        self.name2vert = {}  # (name, vertex): vertex

        self.cell2name = {}  # cellid: (name, cellid)
        self.name2cell = {}  # (name, cellid): cellid

    def get_disv_kwargs(self, name="__main__"):
        return self.get_grid(name).get_disv_kwargs()

    def __repr__(self):
        names = ", ".join(self.grids.keys())
        return f"DisvGridMerger({names})"

    def property_copy_to(self, DisvGridMergerType):
        if isinstance(DisvGridMergerType, DisvGridMerger):
            DisvGridMergerType.merged = self.merged.copy()
            dcp = copy.deepcopy

            for name in self.grids:
                DisvGridMergerType.grids[name] = self.grids[name].copy()

            for name in self.snap_vertices:
                DisvGridMergerType.snap_vertices[name] = dcp(
                    self.snap_vertices[name]
                )

            for name in self.connect_tolerance:
                DisvGridMergerType.connect_tolerance[
                    name
                ] = self.connect_tolerance[name]

            for name in self.force_snap:
                DisvGridMergerType.force_snap[name] = dcp(
                    self.force_snap[name]
                )

            for name in self.force_snap_drop:
                DisvGridMergerType.force_snap_drop[name] = dcp(
                    self.force_snap_drop[name]
                )

            DisvGridMergerType.force_snap_cellid = dcp(self.force_snap_cellid)

            for name in self.vert2name:
                DisvGridMergerType.vert2name[name] = dcp(self.vert2name[name])

            for name in self.cell2name:
                DisvGridMergerType.cell2name[name] = dcp(self.cell2name[name])

            for name in self.name2vert:
                DisvGridMergerType.name2vert[name] = self.name2vert[name]

            for name in self.name2cell:
                DisvGridMergerType.name2cell[name] = self.name2cell[name]

            DisvGridMergerType.snap_order = dcp(self.snap_order)
        else:
            raise RuntimeError(
                "DisvGridMerger.property_copy_to "
                "can only copy to objects that inherit "
                "properties from DisvGridMerger"
            )

    def copy(self):
        cp = DisvGridMerger()
        self.property_copy_to(cp)
        return cp

    def get_merged_cell2d(self, name, cell2d_orig):
        return self.name2cell[(name, cell2d_orig)]

    def get_merged_vertex(self, name, vertex_orig):
        return self.name2vert[(name, vertex_orig)]

    def get_grid(self, name="__main__"):
        if name == "" or name == "__main__":
            return self.merged

        if name not in self.grids:
            raise KeyError(
                "DisvGridMerger.get_grid: requested grid, "
                f"{name} does not exist.\n"
                "Current grids stored are:\n"
                "\n".join(self.grids.keys())
            )

        return self.grids[name]

    def add_grid(self, name, grid):
        if name == "" or name == "__main__":
            raise RuntimeError(
                "\nDisvGridMerger.add_grid:\n"
                'name = "" or "__main__"\nis not allowed.'
            )
        if isinstance(grid, DisvPropertyContainer):
            grid = grid.copy()
        else:
            # grid = [nlay, vertices, cell2d, top, botm]
            grid = DisvPropertyContainer(*grid)

        self.grids[name] = grid

    def set_vertex_connection(
        self, name1, name2, vertex1, vertex2, autosnap_tolerance=1.0e-5
    ):
        if (name2, name1) in self.snap_vertices:
            name1, name2 = name2, name1
            vertex1, vertex2 = vertex2, vertex1

        key1 = (name1, name2)
        key2 = (name2, name1)

        self.snap_vertices[key1] = (vertex1, vertex2)
        self.connect_tolerance[key1] = autosnap_tolerance

        self.force_snap[key1] = [[], []]
        self.force_snap[key2] = [[], []]
        self.force_snap_drop[key1] = []
        self.force_snap_drop[key2] = []

    def set_force_vertex_connection(
        self, name1, name2, vertex1, vertex2, drop_vertex=2
    ):
        key1 = (name1, name2)
        key2 = (name2, name1)
        if key1 not in self.force_snap:
            self.force_snap[key1] = []
            self.force_snap[key2] = []
            self.force_snap_drop[key1] = []
            self.force_snap_drop[key2] = []

        drop_vertex_inv = 1 if drop_vertex == 2 else 2

        self.force_snap[key1][0].append(vertex1)
        self.force_snap[key1][1].append(vertex2)

        self.force_snap[key2][0].append(vertex2)
        self.force_snap[key2][1].append(vertex1)

        self.force_snap_drop[key1].append(drop_vertex)
        self.force_snap_drop[key2].append(drop_vertex_inv)

    def _get_vertex_xy(self, name, iv):
        vertices = self.get_grid(name).vertices
        for iv_orig, xv, yv in vertices:
            if iv == iv_orig:
                return xv, yv
        raise RuntimeError(
            "DisvGridMerger: " f"Failed to find vertex {iv} in grid {name}"
        )

    def _find_merged_vertex(self, xv, yv, tol):
        for iv, xv_chk, yv_chk in self.merged.vertices:
            if abs(xv - xv_chk) + abs(yv - yv_chk) < tol:
                return iv
        return None

    def _replace_vertex_xy(self, iv, xv, yv):
        for vert in self.merged.vertices:
            if iv == vert[0]:
                vert[1] = xv
                vert[2] = yv
                return
        raise RuntimeError(
            "DisvGridMerger: Unknown code error - "
            f"failed to locate vertex {iv}"
        )

    def _clear_attribute(self):
        self.snap_order.clear()
        self.force_snap_cellid.clear()

        self.merged.nlay = 0
        self.merged.nvert = 0
        self.merged.ncpl = 0
        self.merged.cell2d.clear()
        self.merged.vertices.clear()
        self.merged.top = np.array(())
        self.merged.botm.clear()

        self.cell2name.clear()
        self.name2cell.clear()
        self.vert2name.clear()
        self.name2vert.clear()

    def _grid_snap_order(self):
        # grids are snapped to one main grid using key = (name1, name2).
        # it is required that at least name1 or name2 already be defined
        # in the main grid. Function determines order to ensure this.
        #   snap_list -> List of (name1, name2) to parse
        snap_list = list(self.snap_vertices.keys())
        name_used = {snap_list[0][0]}  # First grid name always used first
        snap_order = []  # Final order to build main grid
        snap_append = []  # key's that need to be parsed again
        loop_limit = 50
        loop = 0
        error_msg = (
            "\nDisvGridMerger:\n"
            "Failed to find grid snap order.\n"
            "Snapping must occur in contiguous steps "
            "to the same main, merged grid.\n"
        )
        while len(snap_list) > 0:
            key = snap_list.pop(0)
            name1, name2 = key
            has_name1 = name1 in name_used
            has_name2 = name2 in name_used

            if has_name1 and has_name2:  # grid snapped to main twice
                raise RuntimeError(
                    error_msg + "Once a grid name has been snapped to "
                    "the main grid,\n"
                    "it cannot be snapped again.\n"
                    f"The current snap order determined is:\n\n"
                    f"{snap_order}\n\n"
                    "but the following two grids were already "
                    "snapped to the main grid:\n"
                    f"{name1}\nand\n{name2}\n"
                )

            if has_name1 or has_name2:  # have a name to snap too
                snap_order.append(key)
                name_used.add(name1)
                name_used.add(name2)
            else:  # neither name found, so save for later
                snap_append.append(key)

            if len(snap_list) == 0 and len(snap_append) > 0:
                snap_list.extend(snap_append)
                snap_append.clear()
                loop += 1
                if loop > loop_limit:
                    raise RuntimeError(
                        error_msg + "Determined snap order for the "
                        "main grid is:\n\n"
                        f"{snap_order}\n\n"
                        "but failed to snap the following "
                        "to the main grid:\n\n"
                        f"{snap_list}"
                    )
        return snap_order

    def merge_grids(self):
        self._clear_attribute()

        # First grid is unchanged, all other grids are changed as they
        #   snap to the final merged grid.
        name1 = next(iter(self.snap_vertices.keys()))[0]

        cell2d = self.merged.cell2d
        vertices = self.merged.vertices

        cell2d.extend(copy.deepcopy(self.grids[name1].cell2d))
        vertices.extend(copy.deepcopy(self.grids[name1].vertices))

        for ic, *_ in cell2d:
            self.cell2name[ic] = (name1, ic)
            self.name2cell[(name1, ic)] = ic

        for iv, *_ in vertices:
            self.vert2name[iv] = [(name1, iv)]
            self.name2vert[(name1, iv)] = iv

        ic_new = ic  # Last cell2d id from previous-previous loop
        iv_new = iv  # Last vertex id from previous loop
        snapped = {name1}
        force_snapped = set()
        for key in self._grid_snap_order():  # Loop through vertices to snap
            tol = self.connect_tolerance[key]
            v1, v2 = self.snap_vertices[key]
            name1, name2 = key
            if name2 in snapped:
                name1, name2 = name2, name1
                v1, v2 = v2, v1

            if name1 not in self.snap_order:
                self.snap_order.append(name1)
            if name2 not in self.snap_order:
                self.snap_order.append(name2)

            v1_orig = v1
            v1 = self.name2vert[(name1, v1_orig)]

            v1x, v1y = self._get_vertex_xy("__main__", v1)
            v2x, v2y = self._get_vertex_xy(name2, v2)

            difx = v1x - v2x
            dify = v1y - v2y

            for v2_orig, xv, yv in self.grids[name2].vertices:
                xv += difx
                yv += dify

                if (
                    key in self.force_snap
                    and v2_orig in self.force_snap[key][1]  # force snap vertex
                ):
                    ind = self.force_snap[key][1].index(v2_orig)
                    v1_orig_force = self.force_snap[key][0][ind]
                    iv = self.name2vert[(name1, v1_orig_force)]
                    if self.force_snap_drop[key] == 1:
                        # replace v1 vertex with the x,y from v2 to force snap
                        self._replace_vertex_xy(iv, xv, yv)
                        force_snapped.add((name1, v1_orig_force))
                    else:
                        force_snapped.add((name2, v2_orig))
                else:
                    iv = self._find_merged_vertex(xv, yv, tol)

                if iv is None:
                    iv_new += 1
                    vertices.append([iv_new, xv, yv])
                    self.vert2name[iv_new] = [(name2, v2_orig)]
                    self.name2vert[(name2, v2_orig)] = iv_new
                else:
                    self.vert2name[iv].append((name2, v2_orig))
                    self.name2vert[(name2, v2_orig)] = iv

            # Loop through cells and update center point and icvert
            for ic, xc, yc, ncvert, *icvert in self.grids[name2].cell2d:
                ic_new += 1
                xc += difx
                yc += dify

                self.cell2name[ic_new] = (name2, ic)
                self.name2cell[(name2, ic)] = ic_new

                item = [ic_new, xc, yc, ncvert]  # append new icvert's
                for iv_orig in icvert:
                    item.append(self.name2vert[(name2, iv_orig)])  # = iv_new
                cell2d.append(item)

        # Force snapped cells need to update cell center and check for
        # errors in the grid
        for name, iv_orig in force_snapped:
            for ic, _, _, _, *icvert in self.grids[name].cell2d:
                if iv_orig in icvert:  # cell was deformed by force snap
                    ic_new = self.name2cell[(name, ic)]
                    self.force_snap_cellid.add(ic_new)

        if len(self.force_snap_cellid) > 0:
            dist = lambda v1, v2: sqrt((v2[1]-v1[1])**2 + (v2[2]-v1[2])**2)
            mg = self.merged
            vert_xy = np.array([(x, y) for _, x, y in mg.vertices])
            for ic in self.force_snap_cellid:
                _, _, _, _, *vert = mg.cell2d[ic]
                if len(vert) != len(set(vert)):  # contains a duplicate vertex
                    seen = set()
                    seen_add = seen.add
                    vert = [v for v in vert if not (v in seen or seen_add(v))]
                    tmp = mg.cell2d[ic]
                    tmp[3] = len(vert)
                    mg.cell2d[ic] = tmp[:4] + vert
                # check if vertices are within cell.
                cell = [mg.vertices[iv][1:] for iv in vert]
                path = mpltPath.Path(cell)
                contain = np.where(path.contains_points(vert_xy))[0]
                for iv in contain:
                    # find closest polyline
                    if iv in vert:
                        continue
                    vert_check = mg.vertices[iv]  # vert_xy[iv, :]
                    d = np.inf
                    v_closest = -1
                    for v in vert:
                        d2 = dist(vert_check, mg.vertices[v])
                        if d > d2:
                            d = d2
                            v_closest = v
                    ind = vert.index(v_closest)
                    if ind == len(vert) - 1:  # Closest is at the end
                        d1 = dist(vert_check, mg.vertices[vert[0]])
                        d2 = dist(vert_check, mg.vertices[vert[-2]])
                        if d1 < d2:
                            ind = len(vert)
                    elif ind == 0:  # Closest is at the start check end members
                        d1 = dist(vert_check, mg.vertices[vert[-1]])
                        d2 = dist(vert_check, mg.vertices[vert[1]])
                        if d2 < d1:
                            ind = 1
                    else:
                        d1 = dist(vert_check, mg.vertices[vert[ind-1]])
                        d2 = dist(vert_check, mg.vertices[vert[ind+1]])
                        if d2 < d1:
                            ind += 1

                    # update cell2d for cell ic
                    vert.insert(ind, iv)
                    tmp = mg.cell2d[ic]
                    tmp[3] = len(vert)
                    # update cell center
                    tmp[1], tmp[2] = mg.get_centroid(vert)
                    mg.cell2d[ic] = tmp[:4] + vert

        self.merged.nvert = len(vertices)
        self.merged.ncpl = len(cell2d)

        nlay = 0
        for name in self.snap_order:
            if nlay < self.grids[name].nlay:
                nlay = self.grids[name].nlay
        self.merged.nlay = nlay

        top = []
        for name in self.snap_order:
            top.extend(self.grids[name].top)
        self.merged.top = np.array(top)

        for lay in range(nlay):
            bot = []
            for name in self.snap_order:
                if lay < self.grids[name].nlay:
                    bot.extend(self.grids[name].botm[lay])
                else:
                    bot.extend(self.grids[name].botm[-1])
            self.merged.botm.append(bot)

    def plot_grid(
        self,
        name="__main__",
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
        vertex_dot=True,
        vertex_num=True,
        vertex_dot_size=6.5,
    ):
        return self.get_grid(name).plot_grid(
            title,
            plot_time,
            show,
            figsize,
            dpi,
            xlabel,
            ylabel,
            cell2d_override,
            vertices_override,
            ax_override,
            cell_dot,
            cell_num,
            cell_dot_size,
            vertex_dot,
            vertex_num,
            vertex_dot_size,
        )


if __name__ == "__main__":
    from DisvCurvilinearBuilder import DisvCurvilinearBuilder
    from DisvStructuredGridBuilder import DisvStructuredGridBuilder

    grid1 = DisvStructuredGridBuilder(nlay=1, nrow=2, ncol=2)
    grid2 = DisvStructuredGridBuilder(nlay=1, nrow=3, ncol=2)
    # Optional step to see what vertex point to use
    grid1.plot_grid()  # Plot and view vertex locations
    grid2.plot_grid()  #    to identify connection points.
    # Steps to merge grid1 and grid2
    mg = DisvGridMerger()  # init the object
    mg.add_grid("grid1", grid1)  # add grid1
    mg.add_grid("grid2", grid2)  # add grid2
    # Snap grid1 upper right corner (vertex 3) to grid2 upper left
    #    corner (vertex 1). Note the vertices must be zero-based indexed.
    mg.set_vertex_connection("grid1", "grid2", 3 - 1, 1 - 1)
    # Force snap the bottom left corner of grid2 to grid one.
    #    This will reshape that cell to be a triangle.
    mg.set_force_vertex_connection("grid1", "grid2", 7-1, 10-1)
    # Merge the grids
    mg.merge_grids()
    mg.merged.plot_grid()  # plot the merged grid

    nlay = 1  # Number of layers
    nradial = 5  # Number of radial direction cells (radial bands)
    ncol = 5  # Number of columns in radial band (ncol)

    r_inner = 4  # Model inner radius ($ft$)
    r_outer = 9  # Model outer radius ($ft$)
    r_width = 1  # Model radial band width ($ft$)

    surface_elevation = 10.0  # Top of the model ($ft$)

    radii = np.arange(r_inner, r_outer + r_width, r_width)

    angle_start1 = 180
    angle_stop1 = 270
    angle_step1 = -ncol

    angle_start2 = 0
    angle_stop2 = 90
    angle_step2 = -ncol

    nrow = len(radii) - 1
    row_width = r_width
    col_width = r_width

    curv1 = DisvCurvilinearBuilder(
        nlay,
        radii,
        angle_start1,
        angle_stop1,
        angle_step1,
        surface_elevation=surface_elevation,
        layer_thickness=surface_elevation,
        single_center_cell=False,
    )

    curv2 = DisvCurvilinearBuilder(
        nlay,
        radii,
        angle_start2,
        angle_stop2,
        angle_step2,
        surface_elevation=surface_elevation,
        layer_thickness=surface_elevation,
        single_center_cell=False,
    )

    rect = DisvStructuredGridBuilder(
        nlay,
        nrow,
        ncol,
        row_width,
        col_width,
        surface_elevation,
        surface_elevation,
    )

    fig1, ax1 = curv1.plot_grid(show=False, dpi=150)
    fig2, ax2 = rect.plot_grid(show=False, dpi=150)
    fig3, ax3 = curv2.plot_grid(show=False, dpi=150)

    # fig1.savefig("fig1.png", dpi=600)
    # fig2.savefig("fig2.png", dpi=600)
    # fig3.savefig("fig3.png", dpi=600)

    gm = DisvGridMerger()

    gm.add_grid("curv1", curv1)
    gm.add_grid("rect", rect)
    gm.add_grid("curv2", curv2)

    gm.set_vertex_connection("curv1", "rect", 6 - 1, 1 - 1)
    gm.set_vertex_connection("rect", "curv2", 6 - 1, 36 - 1)

    gm.merge_grids()

    fig4, ax4 = gm.plot_grid(
        show=True,
        dpi=150,
    )


    # vertices = [[iv, xv, yv], ...]
    # cell2d = [[ic, xc, yc, ncvert, icvert], ...]
    
    # cell2d_orig
    # vertices_orig
    # cell2d
    # vertices
    # cell2d_index
    # vertices_index
    # connect_tolerance: distance to autosnap verticies
    # snap -> [(name1, name2)] = (vertex1, vertex2)
    # force_snap -> [(name1, name2)] = (vertex1, vertex2, drop_vertex)
    
    # vert2name -> vertex: [[name, vertex], ...]
    # name2vert -> (name, vertex): vertex
    # cell2name -> cellid: (name, cellid)
    # name2cell -> (name, cellid): cellid