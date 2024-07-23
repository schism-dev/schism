"""This script creates a DEM source map for the DEM-loaded hgrid,
the output is a shapefile with each polygon representing a DEM source area."""


import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
from tqdm import tqdm
from pylib_experimental.schism_file import cread_schism_hgrid


def create_polygons(points: np.ndarray, elements: np.ndarray):
    """
    Function to create polygons from unstructured grid data.
    Each polygon is a contiguous area of the same nodal z value.

    Inputs:
        points: np.ndarray, shape (n, 3), where n is the number of points and the
                columns are (x, y, value)
        elements: np.ndarray, shape (m, 4), where m is the number of elements and the
                columns are (node1, node2, node3, node4)

    Outputs:
        polygons: list of tuples, each tuple contains a Polygon and the value of the points that make up the polygon
    """

    polygons = []
    values = []
    print('Creating polygons...')
    for element in tqdm(elements):
        nodes_indices = element[:4]
        nodes_indices = nodes_indices[nodes_indices >= 0]
        value = points[nodes_indices[0], 2]  # take the first value
        vertices = points[nodes_indices]
        polygon = Polygon(vertices)
        polygons.append(polygon)
        values.append(value)

    gdf = gpd.GeoDataFrame({'geometry': polygons, 'value': values})
    dissolved_gdf= gdf.dissolve(by='value')
    dissolved_gdf['dem_id'] = dissolved_gdf.index  # add dem_id column, same as index (values)

    return dissolved_gdf


def make_dem_dict(working_dir: str):
    """Make a dictionary of shortend dem_name and shp_value,
    where shp_value is a unique value for each dem source area."""

    # dem_id at each node
    hgrid_dem_id_file = f'{working_dir}/hgrid.ll.dem_loaded.gr3_dem_id'
    dem_id_nodal = np.loadtxt(hgrid_dem_id_file).astype(int)

    # dem_name for each dem_id
    dem_names_file = f'{working_dir}/hgrid.ll.dem_loaded.gr3_dem_name'
    dem_name_dict = {}
    with open(dem_names_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip().replace(':', ' ').split()
            dem_name_dict[int(line[0])] = line[1]

    # condense the dem_name_dict by treating similar names as the same
    dem_name_dict_condensed = {}
    for dem_id, dem_name in dem_name_dict.items():
        dem_name_condensed = dem_name.replace('.', '_').split('_')[0]  # the first part of the name
        if dem_name_condensed not in dem_name_dict_condensed.keys():
            dem_name_dict_condensed[dem_id] = dem_name_condensed
    # assign a unique value to each shortened dem_name
    dem2num_condensed = {}
    num2dem_condensed = {}
    shp_value = 0
    for _, value in dem_name_dict_condensed.items():
        if value not in dem2num_condensed.keys():
            dem2num_condensed[value] = shp_value
            num2dem_condensed[shp_value] = value
            shp_value += 1

    # map each node dem id to the condensed dem name then to the shp_value
    dem_name_nodal = [dem_name_dict[id] for id in dem_id_nodal]
    dem_short_name_nodal = [dem_name.replace('.', '_').split('_')[0] for dem_name in dem_name_nodal]
    dem_num_nodal = [dem2num_condensed[dem_name] for dem_name in dem_short_name_nodal]

    return dem_num_nodal, num2dem_condensed


def make_dem_source_map():
    """Master function to create the DEM source map."""

    WDIR = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12/Bathy_edit/DEM_loading/Original/'

    hgrid = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12/hgrid.ll')

    dem_id_nodal, dem_dict = make_dem_dict(working_dir=WDIR)
    gdf = create_polygons(np.c_[hgrid.x, hgrid.y, dem_id_nodal], hgrid.elnode)

    # add name attribute to each polygon
    for i, poly in gdf.iterrows():
        gdf.loc[i, 'dem_name'] = dem_dict[poly['dem_id']]

    gdf.to_file(f'{WDIR}/dem_source_map.shp')
    print('DEM source map created')


if __name__ == '__main__':
    make_dem_source_map()
