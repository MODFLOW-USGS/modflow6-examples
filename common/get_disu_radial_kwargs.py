import numpy as np


def get_disu_radial_kwargs(
    nlay,
    nradial,
    radius_outer,
    surface_elevation,
    layer_thickness,
    get_vertex=False,
):
    """
    Simple utility for creating radial unstructured elements
    with the disu package.

    Input assumes that each layer contains the same radial band,
    but their thickness can be different.

    Parameters
    ----------
    nlay: number of layers (int)
    nradial: number of radial bands to construct (int)
    radius_outer: Outer radius of each radial band (array-like float with nradial length)
    surface_elevation: Top elevation of layer 1 as either a float or nradial array-like float values.
                       If given as float, then value is replicated for each radial band.
    layer_thickness: Thickness of each layer as either a float or nlay array-like float values.
                     If given as float, then value is replicated for each layer.
    """
    pi = 3.141592653589793

    def get_nn(lay, rad):
        return nradial * lay + rad

    def get_rad_array(var, rep):
        try:
            dim = len(var)
        except:
            dim, var = 1, [var]

        if dim != 1 and dim != rep:
            raise IndexError(
                f"get_rad_array(var): var must be a scalar or have len(var)=={rep}"
            )

        if dim == 1:
            return np.full(rep, var[0], dtype=np.float64)
        else:
            return np.array(var, dtype=np.float64)

    nodes = nlay * nradial
    surf = get_rad_array(surface_elevation, nradial)
    thick = get_rad_array(layer_thickness, nlay)

    iac = np.zeros(nodes, dtype=int)
    ja = []
    ihc = []
    cl12 = []
    hwva = []

    area = np.zeros(nodes, dtype=float)
    top = np.zeros(nodes, dtype=float)
    bot = np.zeros(nodes, dtype=float)

    for lay in range(nlay):
        st = nradial * lay
        sp = nradial * (lay + 1)
        top[st:sp] = surf - thick[:lay].sum()
        bot[st:sp] = surf - thick[: lay + 1].sum()

    for lay in range(nlay):
        for rad in range(nradial):
            # diagonal/self
            n = get_nn(lay, rad)
            ja.append(n)
            iac[n] += 1
            if rad > 0:
                area[n] = pi * (
                    radius_outer[rad] ** 2 - radius_outer[rad - 1] ** 2
                )
            else:
                area[n] = pi * radius_outer[rad] ** 2
            ihc.append(n + 1)
            cl12.append(n + 1)
            hwva.append(n + 1)
            # up
            if lay > 0:
                ja.append(n - nradial)
                iac[n] += 1
                ihc.append(0)
                cl12.append(0.5 * (top[n] - bot[n]))
                hwva.append(area[n])
            # to center
            if rad > 0:
                ja.append(n - 1)
                iac[n] += 1
                ihc.append(1)
                cl12.append(0.5 * (radius_outer[rad] - radius_outer[rad - 1]))
                hwva.append(2.0 * pi * radius_outer[rad - 1])

            # to outer
            if rad < nradial - 1:
                ja.append(n + 1)
                iac[n] += 1
                ihc.append(1)
                hwva.append(2.0 * pi * radius_outer[rad])
                if rad > 0:
                    cl12.append(
                        0.5 * (radius_outer[rad] - radius_outer[rad - 1])
                    )
                else:
                    cl12.append(radius_outer[rad])
            # bottom
            if lay < nlay - 1:
                ja.append(n + nradial)
                iac[n] += 1
                ihc.append(0)
                cl12.append(0.5 * (top[n] - bot[n]))
                hwva.append(area[n])

    # Build rectangular equivalent of radial coordinates (unwrap radial bands)
    if get_vertex:
        perimeter_outer = np.fromiter(
            (2.0 * pi * rad for rad in radius_outer),
            dtype=float,
            count=nradial,
        )
        xc = 0.5 * radius_outer[0]
        yc = 0.5 * perimeter_outer[-1]
        # all cells have same y-axis cell center; yc is costant
        #
        # cell2d: [icell2d, xc, yc, ncvert, icvert]; first node: cell2d = [[0, xc, yc, [2, 1, 0]]]
        cell2d = []
        for lay in range(nlay):
            n = get_nn(lay, 0)
            cell2d.append([n, xc, yc, 3, 2, 1, 0])
        #
        xv = radius_outer[0]
        # half perimeter is equal to the y shift for vertices
        sh = 0.5 * perimeter_outer[0]
        vertices = [
            [0, 0.0, yc],
            [1, xv, yc - sh],
            [2, xv, yc + sh],
        ]  # vertices: [iv, xv, yv]
        iv = 3
        for r in range(1, nradial):
            # radius_outer[r-1] + 0.5*(radius_outer[r] - radius_outer[r-1])
            xc = 0.5 * (radius_outer[r - 1] + radius_outer[r])
            for lay in range(nlay):
                n = get_nn(lay, r)
                # cell2d: [icell2d, xc, yc, ncvert, icvert]
                cell2d.append([n, xc, yc, 4, iv - 2, iv - 1, iv + 1, iv])

            xv = radius_outer[r]
            # half perimeter is equal to the y shift for vertices
            sh = 0.5 * perimeter_outer[r]
            vertices.append([iv, xv, yc - sh])  # vertices: [iv, xv, yv]
            iv += 1
            vertices.append([iv, xv, yc + sh])  # vertices: [iv, xv, yv]
            iv += 1
        cell2d.sort(key=lambda row: row[0])  # sort by node number

    ja = np.array(ja, dtype=np.int32)
    nja = ja.shape[0]
    hwva = np.array(hwva, dtype=np.float64)
    kw = {}
    kw["nodes"] = nodes
    kw["nja"] = nja
    kw["nvert"] = None
    kw["top"] = top
    kw["bot"] = bot
    kw["area"] = area
    kw["iac"] = iac
    kw["ja"] = ja
    kw["ihc"] = ihc
    kw["cl12"] = cl12
    kw["hwva"] = hwva

    if get_vertex:
        kw["nvert"] = len(vertices)  # = 2*nradial + 1
        kw["vertices"] = vertices
        kw["cell2d"] = cell2d
        kw["angldegx"] = np.zeros(nja, dtype=float)
    else:
        kw["nvert"] = 0

    return kw


if __name__ == "__main__":
    # Test case

    nlay = 25  # Number of layers
    nradial = 22  # number of radial bands (innermost band is a full circle)

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
    ]  # outer radius of each radial band

    layer_thickness = 2.0  # thickness of each radial layer
    surface_elevation = 50.0  # surface elevation model

    radial_kwargs = get_disu_radial_kwargs(
        nlay,
        nradial,
        radius_outer,
        surface_elevation,
        layer_thickness,
        get_vertex=True,
    )

    print(radial_kwargs)
