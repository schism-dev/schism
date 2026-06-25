import argparse
import copy
from time import time
import multiprocessing as mp

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid
from geopandas import GeoDataFrame

import utils


# Operational parameters
EXPECTED_ELEVATION_FILL_VALUE = -99999.0
DRY_NODE_TOLERANCE = 1.e-6
LAND_DISTURBANCE_MASK_THRESHOLD = -5
TRIANGULATION_MASK_THRESHOLD = -90000

VERTICAL_DATUM = 'xGEOID20b'
UNITS = 'meters'
OUTPUT_DRIVER = 'GPKG'
OUTPUT_LAYER = 'disturbance'
OUTPUT_PREFIX = 'stofs_3d_atl.t12z.disturbance'
NOWCAST_OUTPUT_COUNT = 24
FORECAST_OUTPUT_COUNT = 96

CONTOUR_LEVEL_START = -0.5
CONTOUR_LEVEL_STOP = 2.1
CONTOUR_LEVEL_STEP = 0.1
CONTOUR_LEVEL_MIN = -5
CONTOUR_LEVEL_MAX = 20
CRS_EPSG = 4326
DEFAULT_FACE_COLOR = [0, 0, 0, 1]


def get_disturbance(elevation, depth, levels, fill_value, out_filename):

    # disturbance
    disturbance = copy.deepcopy(elevation)
    idxs_land_node = depth < 0
    disturbance[idxs_land_node] = np.maximum(
        0, elevation[idxs_land_node] + depth[idxs_land_node]
    )

    # set mask on dry nodes
    idxs_dry = np.where(elevation + depth <= DRY_NODE_TOLERANCE)
    disturbance[idxs_dry] = fill_value

    # set mask on nodes with small max disturbance (<-5 m) on land
    idxs_small = disturbance < LAND_DISTURBANCE_MASK_THRESHOLD
    idxs_mask_maxdist = idxs_small * idxs_land_node
    disturbance[idxs_mask_maxdist] = fill_value

    # get mask for triangulation
    imask = disturbance < TRIANGULATION_MASK_THRESHOLD
    mask = np.any(
        np.where(imask[triangulation.triangles], True, False), axis=1
    )
    triangulation.set_mask(mask)
    gdf = contour_to_gdf(disturbance, levels, triangulation)

    # gdf.to_file(out_filename, driver="GeoJSON")
    gdf.to_file(out_filename, driver=OUTPUT_DRIVER, layer=OUTPUT_LAYER)
    # gdf.to_file(out_filename)


def contour_to_gdf(disturbance, levels, triangulation):

    MinVal = levels[0]
    MaxVal = levels[-1]

    # MinMax = [(-99999, levels[0])]
    MinMax = [(disturbance.min(), levels[0])]
    for i in np.arange(len(levels) - 1):
        MinMax.append((levels[i], levels[i + 1]))
    # MinMax.append((levels[-1], np.max(disturbance)))

    fig = plt.figure()
    ax = fig.add_subplot()

    my_cmap = plt.cm.jet
    contour = ax.tricontourf(
        triangulation, disturbance, vmin=MinVal, vmax=MaxVal,
        levels=levels, cmap=my_cmap, extend='min')

    polygons, colors = [], []
    data = []

    # Get paths and facecolors directly from the ContourSet
    paths = contour.get_paths()
    facecolors = contour.get_facecolor()

    # Each path corresponds to a specific contour level/interval
    for i, path in enumerate(paths):
        mpoly = []
        try:
            path.should_simplify = False
            poly_list = path.to_polygons()

            if len(poly_list) > 0:
                # Sort by area descending to ensure the largest polygon is
                # the exterior ring.
                # (Matplotlib occasionally returns holes before exteriors)
                poly_list = sorted(
                    poly_list,
                    key=lambda x: Polygon(x).area if len(x) > 3 else 0,
                    reverse=True,
                )

                exterior = poly_list[0]
                holes = [h for h in poly_list[1:] if len(h) > 3]

                if len(exterior) > 3:
                    mpoly.append(make_valid(Polygon(exterior, holes)))
        except Exception as e:
            print(f'Warning: Geometry error when making polygon #{i}: {e}')

        # Skip adding to data if no valid geometry was generated
        if not mpoly:
            continue

        # Combine into MultiPolygon if multiple discrete shapes exist for this
        # level.
        geom = MultiPolygon(mpoly) if len(mpoly) > 1 else mpoly[0]

        polygons.append(geom)

        # Safely extract facecolor for this level
        color = (
            facecolors[i].tolist()
            if i < len(facecolors) else DEFAULT_FACE_COLOR
        )
        colors.append(color)

        data.append({
            'id': i + 1,
            'minWaterLevel': MinMax[i][0],
            'maxWaterLevel': MinMax[i][1],
            'verticalDatum': VERTICAL_DATUM,
            'units': UNITS,
            'geometry': geom
        })

    plt.close('all')

    gdf = GeoDataFrame(data)

    # Get color in Hex
    colors_elev = []
    my_cmap = plt.cm.jet

    for i in range(len(gdf)):
        color = my_cmap(i / len(gdf))
        colors_elev.append(mpl.colors.to_hex(color))

    gdf['rgba'] = colors_elev

    # set crs
    gdf = gdf.set_crs(CRS_EPSG)

    return gdf


def parse_args():
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--input_filename", required=True, help="file name in SCHISM format"
    )
    return argparser.parse_args()


def read_grid(ds):
    x = ds['SCHISM_hgrid_node_x'][:]
    y = ds['SCHISM_hgrid_node_y'][:]
    depth = ds['depth'][:]
    elements = ds['SCHISM_hgrid_face_nodes'][:, :]

    return x, y, depth, elements


def build_triangulation(x, y, elements):
    t0 = time()
    tris = utils.split_quads(elements)
    print(f'Spliting quads took {time() - t0} seconds')

    return utils.triangulation(x, y, tris)


def validate_elevation_fill_value(ds):
    elevation_fill_value = ds["elevation"].missing_value
    if not np.isclose(elevation_fill_value, EXPECTED_ELEVATION_FILL_VALUE):
        raise ValueError(
            "Expected elevation:missing_value to be "
            f"{EXPECTED_ELEVATION_FILL_VALUE}, "
            f"but got {elevation_fill_value}. "
            "The masking logic assumes this fill value."
        )

    return elevation_fill_value


def get_pond_mask(ds, elev_shape):
    if "isolatedPondNode" in ds.variables:
        print(
            "Found variable 'isolatedPondNode' in dataset, applying pond mask"
        )
        return ds["isolatedPondNode"][:, :].astype(bool)

    print(
        "Variable 'isolatedPondNode' not found in dataset, "
        "proceeding without pond mask"
    )
    return np.zeros(elev_shape, dtype=bool)


def read_masked_elevation(ds):
    elevation_fill_value = validate_elevation_fill_value(ds)

    elev = ds['elevation'][:, :]
    dryFlagNode = ds['dryFlagNode'][:, :]
    pond_mask = get_pond_mask(ds, elev.shape)

    elev[dryFlagNode[:, :] == 1] = elevation_fill_value
    elev[np.isnan(elev)] = elevation_fill_value
    elev[abs(elev) >= abs(elevation_fill_value)] = elevation_fill_value
    elev[pond_mask] = elevation_fill_value

    return elev, elevation_fill_value


def get_contour_levels():
    levels = [
        round(i, 1)
        for i in np.arange(
            CONTOUR_LEVEL_START, CONTOUR_LEVEL_STOP, CONTOUR_LEVEL_STEP
        )
    ]
    levels.append(CONTOUR_LEVEL_MAX)
    levels.insert(0, CONTOUR_LEVEL_MIN)
    # print(levels)
    # print(f'Calculating and masking disturbance took {time() - t0} seconds')

    return levels


def get_output_filenames():
    filenames = [
        f'{OUTPUT_PREFIX}.n{i - 1:03d}.gpkg'
        for i in range(NOWCAST_OUTPUT_COUNT, 0, -1)
    ]
    for i in range(FORECAST_OUTPUT_COUNT):
        filenames.append(f'{OUTPUT_PREFIX}.f{i + 1:03d}.gpkg')
    # print(filenames)

    return filenames


def write_hourly_disturbances(
    elev, depth, times, levels, fill_value, filenames
):
    npool = len(times) if len(times) < mp.cpu_count() else mp.cpu_count() - 1
    # print(npool)

    t0 = time()

    pool = mp.Pool(npool)
    pool.starmap(
        get_disturbance,
        [
            (
                np.squeeze(elev[i, :]),
                depth,
                levels,
                fill_value,
                filenames[i],
            )
            for i in range(len(times))
        ],
    )

    pool.close()
    del pool

    print(
        'Calculating and masking disturbance for all times took '
        f'{time() - t0} seconds'
    )


def main():
    global triangulation

    args = parse_args()
    input_filename = args.input_filename

    ds = Dataset(input_filename)
    x, y, depth, elements = read_grid(ds)
    triangulation = build_triangulation(x, y, elements)

    times = ds['time'][:]

    elev, elevation_fill_value = read_masked_elevation(ds)
    levels = get_contour_levels()
    filenames = get_output_filenames()

    write_hourly_disturbances(
        elev, depth, times, levels, elevation_fill_value, filenames
    )


if __name__ == "__main__":
    main()
