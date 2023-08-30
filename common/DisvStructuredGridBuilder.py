import numpy as np
from DisvPropertyContainer import DisvPropertyContainer

__all__ = ["DisvStructuredGridBuilder"]


class DisvStructuredGridBuilder(DisvPropertyContainer):
    """
    A class for generating a structured MODFLOW 6 DISV grid.

    This class inherits from the `DisvPropertyContainer` class and provides
    methods to generate a rectangular, structured grid give the number rows
    (nrow), columns (ncol), row widths, and columns widths. Rows are
    discretized along the y-axis and columns along the x-axis. The row, column
    structure follows MODFLOW 6 structured grids. That is, the first row has
    the largest y-axis vertices and last row the lowest; and the first column
    has the lowest x-axis vertices and last column the highest.

    All indices are zero-based, but translated to one-base for the figures and
    by flopy for use with MODFLOW 6.

    The following shows the placement for (row, column) pairs
    in a nrow=3 and ncol=5 model:

    ``(0,0)  (0,1)  (0,2)  (0,3)  (0,4)``
    ``(1,0)  (1,1)  (1,2)  (1,3)  (1,4)``
    ``(2,0)  (2,1)  (2,2)  (2,3)  (2,4)``

    Array-like structures that are multidimensional (has rows and columns)
    are flatten by concatenating each row. Using the previous example,
    the following is the flattened representation:

    ``(0,0) (0,1) (0,2) (0,3) (0,4) (1,0) (1,1) (1,2) (1,3) (1,4) (2,0) (2,1) (2,2) (2,3) (2,4)``

    If no arguments are provided then an empty object is returned.

    Parameters
    ----------
    nlay : int
        Number of layers
    nrow : int
        Number of rows (y-direction cells).
    ncol : int
        Number of columns (x-direction cells).
    row_width : float or array_like
        Width of y-direction cells (each row). If a single value is provided,
        it will be used for all rows. Otherwise, it must be array_like
        of length ncol.
    col_width : float or array_like
        Width of x-direction cells (each column). If a single value is
        provided, it will be used for all columns. Otherwise, it must be
        array_like of length ncol.
    surface_elevation : float or array_like
        Surface elevation for the top layer. Can either be a single float
        for the entire `top`, or array_like of length `ncpl`.
        If it is a multidimensional array_like, then it is flattened to a
        single dimension along the rows (first dimension).
    layer_thickness : float or array_like
        Thickness of each layer. Can either be a single float
        for model cells, or array_like of length `nlay`, or
        array_like of length `nlay`*`ncpl`.
    origin_x : float, default=0.0
        X-coordinate reference point the lower-left corner of the model grid.
        That is, the outermost corner of ``row=nrow-1` and `col=0`.
        Rotations are performed around this point.
    origin_y : float, default=0.0
        Y-coordinate reference point the lower-left corner of the model grid.
        That is, the outermost corner of ``row=nrow-1` and `col=0`.
        Rotations are performed around this point.
    rotation : float, default=0.0
        Rotation angle in degrees for the model grid around (origin_x, origin_y).

    Attributes
    ----------
    nrow : int
        Number of rows in the grid.
    ncol : int
        Number of columns in the grid.
    row_width : np.ndarray
        Width of y-direction cells (each row).
    col_width : np.ndarray
        Width of x-direction cells (each column).

    Methods
    -------
    get_disv_kwargs()
        Get the keyword arguments for creating a MODFLOW-6 DISV package.
    plot_grid(...)
        Plot the model grid from `vertices` and `cell2d` attributes.
    get_cellid(row, col, col_check=True)
        Get the cellid given the row and column indices.
    get_row_col(cellid)
        Get the row and column indices given the cellid.
    get_vertices(row, col)
        Get the vertex indices for a cell given the row and column indices.
    iter_row_col()
        Iterate over row and column indices for each cell.
    iter_row_cellid(row)
        Iterate over cellid's in a specific row.
    iter_column_cellid(col)
        Iterate over cellid's in a specific column.
    """

    nrow: int
    ncol: int
    row_width: np.ndarray
    col_width: np.ndarray

    def __init__(
        self,
        nlay=-1,
        nrow=-1,  # number of Y direction cells
        ncol=-1,  # number of X direction cells
        row_width=10.0,  # width of Y direction cells (each row)
        col_width=10.0,  # width of X direction cells (each column)
        surface_elevation=100.0,
        layer_thickness=100.0,
        origin_x=0.0,
        origin_y=0.0,
        rotation=0.0,
    ):
        if nlay is None or nlay < 1:
            self._init_empty()
            return

        ncpl = ncol * nrow  # Nodes per layer

        self.nrow = nrow
        self.ncol = ncol
        ncell = ncpl * nlay

        # Check if layer_thickness needs to be flattened
        cls_name = "DisvStructuredGridBuilder"
        top = self._get_array(cls_name, surface_elevation, ncpl)
        thick = self._get_array(cls_name, layer_thickness, nlay, ncell)
        self.row_width = self._get_array(cls_name, row_width, ncol)
        self.col_width = self._get_array(cls_name, col_width, nrow)

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

        # Build the grid

        # Setup vertices
        vertices = []

        # Get row 1 top:
        yv_model_top = self.col_width.sum()

        # Assemble vertices along x-axis and model top
        iv = 0
        xv, yv = 0.0, yv_model_top
        vertices.append([iv, xv, yv])
        for c in range(ncol):
            iv += 1
            xv += self.row_width[c]
            vertices.append([iv, xv, yv])

        # Finish the rest of the grid a row at a time
        for r in range(nrow):
            iv += 1
            yv -= self.col_width[r]
            xv = 0.0
            vertices.append([iv, xv, yv])
            for c in range(ncol):
                iv += 1
                xv += self.row_width[c]
                vertices.append([iv, xv, yv])

        # cell2d: [icell2d, xc, yc, ncvert, icvert]
        cell2d = []
        ic = -1
        # Finish the rest of the grid a row at a time
        for r in range(nrow):
            for c in range(ncol):
                ic += 1
                icvert = self.get_vertices(r, c)
                xc, yc = self.get_centroid(icvert, vertices)
                cell2d.append([ic, xc, yc, 4, *icvert])

        super().__init__(
            nlay, vertices, cell2d, top, bot, origin_x, origin_y, rotation
        )

    def __repr__(self):
        return super().__repr__("DisvStructuredGridBuilder")

    def _init_empty(self):
        super()._init_empty()
        nul = np.array([])
        self.nrow = 0
        self.ncol = 0
        self.row_width = nul
        self.col_width = nul

    def property_copy_to(self, DisvStructuredGridBuilderType):
        if isinstance(
            DisvStructuredGridBuilderType, DisvStructuredGridBuilder
        ):
            super().property_copy_to(DisvStructuredGridBuilderType)
            DisvStructuredGridBuilderType.nrow = self.nrow
            DisvStructuredGridBuilderType.ncol = self.ncol
            DisvStructuredGridBuilderType.row_width = self.row_width.copy()
            DisvStructuredGridBuilderType.col_width = self.col_width.copy()
        else:
            raise RuntimeError(
                "DisvStructuredGridBuilder.property_copy_to "
                "can only copy to objects that inherit "
                "properties from DisvStructuredGridBuilder"
            )

    def copy(self):
        cp = DisvStructuredGridBuilder()
        self.property_copy_to(cp)
        return cp

    def get_cellid(self, row, col):
        """
        Get the cellid given the row and column indices.

        Parameters
        ----------
        row : int
            Row index.
        col : int
            Column index.

        Returns
        -------
        int
            cellid index
        """
        return row * self.ncol + col

    def get_row_col(self, cellid):
        """
        Get the row and column indices given the cellid.

        Parameters
        ----------
        cellid : int
            cellid index

        Returns
        -------
        (int, int)
            Row index, Column index
        """
        row = cellid // self.ncol
        col = cellid - row * self.ncol
        return row, col

    def get_vertices(self, row, col):
        """
        Get the vertex indices for a cell given the row and column indices.

        Parameters
        ----------
        row : int
            Row index.
        col : int
            Column index.

        Returns
        -------
        list[int]
            List of vertex indices that define the cell at (row, col).
        """
        nver = self.ncol + 1
        return [
            row * nver + col,
            row * nver + col + 1,
            (row + 1) * nver + col + 1,
            (row + 1) * nver + col,
        ]

    def iter_row_col(self):
        """Generator that iterates through each rows' columns.

        Yields
        -------
        (int, int)
            Row index, column index
        """
        for cellid in range(self.ncpl):
            yield self.get_row_col(cellid)

    def iter_row_cellid(self, row):
        """Generator that iterates through the cellid within a row.
        That is, the cellid for all columns within the specified row.

        Parameters
        ----------
        row : int
            Row index.

        Yields
        -------
        int
            cellid index
        """
        for col in range(self.ncol):
            yield self.get_cellid(row, col)

    def iter_column_cellid(self, col):
        """Generator that iterates through the cellid within a column.
        That is, the cellid for all rows within the specified column.

        Parameters
        ----------
        col : int
            Column index.

        Yields
        -------
        int
            cellid index
        """
        for row in range(self.nrow):
            yield self.get_cellid(row, col)


if __name__ == "__main__":
    simple2by2 = DisvStructuredGridBuilder(1, 2, 2)
    simple2by3 = DisvStructuredGridBuilder(1, 2, 3)
    simple2by2.plot_grid()
    simple2by3.plot_grid()

    test1 = DisvStructuredGridBuilder(1, 2, 3, 10.0, 10.0, 100.0, 100)
    test1.plot_grid(plot_time=0.2)

    test2 = DisvStructuredGridBuilder(
        4, 5, 7, 50.0, 50.0, 100.0, 25.0, 10, 10, 30.0
    )
    test2.plot_grid(plot_time=0.2)
    pass
