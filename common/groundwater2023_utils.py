import numpy as np
import shapely

geometries = {
    "sv_boundary": """0.0 0.0
    0.0 20000.0
    12500.0 20000.0
    12500.0 0.0""",
    "sv_river": """4250.0 8750.0 
    4250.0 0.0""",
    "sv_river_box": """3500.0 0.0
    3500.0 9500.0
    5000.0 9500.0
    5000.0 0.0""",
    "sv_wells": """7250. 17250.
    7750. 2750.
    2750 3750.""",
    "sv_lake": """1500. 18500.
    3500. 18500.
    3500. 15500.
    4000. 15500.
    4000. 14500.
    4500. 14500.
    4500. 12000.
    2500. 12000.
    2500. 12500.
    2000. 12500.
    2000. 14000.
    1500. 14000.
    1500. 15000.
    1000. 15000.
    1000. 18000.
    1500. 18000.""",
}


def string2geom(geostring, conversion=None):
    if conversion is None:
        multiplier = 1.0
    else:
        multiplier = float(conversion)
    res = []
    for line in geostring.split("\n"):
        line = line.strip()
        line = line.split(" ")
        x = float(line[0]) * multiplier
        y = float(line[1]) * multiplier
        res.append((x, y))
    return res


def densify_geometry(line, step, keep_internal_nodes=True):
    xy = []  # list of tuple of coordinates
    lines_strings = []
    if keep_internal_nodes:
        for idx in range(1, len(line)):
            lines_strings.append(shapely.geometry.LineString(line[idx - 1 : idx + 1]))
    else:
        lines_strings = [shapely.geometry.LineString(line)]

    for line_string in lines_strings:
        length_m = line_string.length  # get the length
        for distance in np.arange(0, length_m + step, step):
            point = line_string.interpolate(distance)
            xy_tuple = (point.x, point.y)
            if xy_tuple not in xy:
                xy.append(xy_tuple)
        # make sure the end point is in xy
        if keep_internal_nodes:
            xy_tuple = line_string.coords[-1]
            if xy_tuple not in xy:
                xy.append(xy_tuple)

    return xy


def circle_function(center=(0, 0), radius=1.0, dtheta=10.0):
    angles = np.arange(0.0, 360.0, dtheta) * np.pi / 180.0
    xpts = center[0] + np.cos(angles) * radius
    ypts = center[1] + np.sin(angles) * radius
    return np.array([(x, y) for x, y in zip(xpts, ypts)])
