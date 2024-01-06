#!/usr/bin/env python3

'''
This script generates a polygon shapefile covering the spatial domain of hgrid.*,
with each polygon showing the corresponding source DEM.
'''
import numpy as np
import pandas as pd
from pylib import schism_grid
import numpy as np
import json
from pylib_essentials.schism_file import read_schism_hgrid_cached
from shapely.geometry import Polygon
import geopandas as gpd
from pathlib import Path


def x_in_y(x, y):
    '''
    Find the index of x in y, where x and y are numpy 1D arrays.
    '''
    x = np.array(x)
    y = np.array(y)

    index = np.argsort(x)
    sorted_x = x[index]
    sorted_index = np.searchsorted(sorted_x, y)

    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y

    result = np.ma.array(yindex, mask=mask)
    return result

if __name__ == "__main__":
    wdir = './'  # this is the working directory of pload*.py, which should contain
                 # hgrid.ll, hgrid.ll.new_dem_id, hgrid.ll.new_dem_name, and DEM_info.json
    dem_info_dict = json.load(open(fr'{wdir}/DEM_info.json'))  # existing json containing the DEM info
    dem_id_shp_fname = Path(fr'{wdir}/hgrid.dem_id.shp')  # this is the end product of this script
    attributes = ['Note', 'Annotation', 'File']

    if not dem_id_shp_fname.exists():  # this is the end product of this script

        # outputs from pload*.py
        dem_tile_name_df = pd.read_csv(fr'{wdir}/hgrid.ll.new_dem_name', sep=':', header=None)
        hg_dem_id = np.loadtxt(fr'{wdir}/hgrid.ll.new_dem_id')
        hg = schism_grid(fr'{wdir}/hgrid.ll')

        # make a dict of dem tile ids and the tile names
        tile_names_dict = {dem_tile_name_df[0][i]: dem_tile_name_df[1][i] for i in range(len(dem_tile_name_df))}
        tile_name2dem_id = {}
        tile_id2dem_id = {}
        for tile_id, tile_name in tile_names_dict.items():  # i.e., tile id as arranged by pload*py, and file names of each tile
            for dem_id, dem_ident_string in enumerate(list(dem_info_dict.keys())):  # i.e., the identifying string of each DEM
                if dem_ident_string.lower() in tile_name.lower():
                    tile_name2dem_id[tile_name] = dem_id
                    tile_id2dem_id[tile_id] = dem_id

        # write the index of the unique prefix to *.2dm
        hg.dp = hg_dem_id.astype(int)
        hg.dp = np.array([tile_id2dem_id[x] for x in hg.dp])
        # hg.grd2sms(fr'{wdir}/hgrid.dem_id.2dm')

        nodes = np.c_[hg.x, hg.y, hg.dp]
        nodes = np.r_[nodes, np.c_[np.nan, np.nan, np.nan]]  # add a dummy node to accomodate for -1 in elnode
        # replace -2 with -1 to accomodate for the dummy node id of triangles in elnode
        hg.elnode[hg.elnode == -2] = -1
        # map the nodal values to elements
        elnode_dp = nodes[:, 2][hg.elnode]

        # define interface elements, this is needed because the DEM interpolation is node based
        # option 1: set interface elements to -999
        # for n in [3, 4]:  # triangle and quad
        #     idx = np.argwhere(hg.i34 == n).flatten()
        #     valid = np.all(elnode_dp[idx, 1:n] == elnode_dp[idx, :n-1], axis=1)  # check if all nodes have the same z
        #     elnode_dp[idx[~valid], :] = -999
        # dpe = elnode_dp[:, 0].astype(int)
        # option 2: set interface elements to the max of the nodal values
        dpe = np.nanmax(elnode_dp, axis=1).astype(int)

        # assemble triangles and quads into polygons dataframe
        polygons = []
        for i, [i34, elnode] in enumerate(zip(hg.i34, hg.elnode)):
            elnode = elnode[:i34]
            element_coords = nodes[elnode, :2]
            polygons.append(Polygon(element_coords))
        gdf = gpd.GeoDataFrame({'geometry': polygons, 'z': dpe})

        # merge polygons with the same elemental z
        gdf = gdf.dissolve(by='z')

        # add an attribute for the dem_prefix (i.e., dem_ident_string in DEM_info.json) based on z
        gdf['dem_prefix'] = [list(dem_info_dict.keys())[x] for x in gdf.index]

    else:
        # read the final shapefile and make changes
        gdf = gpd.read_file(dem_id_shp_fname)
        # delete columns except 'z', 'dem_prefix', and 'geometry'
        gdf.drop(columns=gdf.columns.difference(['z', 'dem_prefix', 'geometry']), inplace=True)

    # add an attribute of the links to the dem notes
    # match the dem_prefix with the dem_info_dict
    for key, value in dem_info_dict.items():  # loop through the DEM identifying string in the dem_info_dict
        # match the dem_prefix if the key is a substring of the dem_prefix case insensitive
        for prefix in gdf['dem_prefix'].unique():
            if key.lower() in prefix.lower():  # found a match in the dem_prefix
                for attr in attributes:  # push the needed attributes from the DEM info to the gdf
                    gdf.loc[gdf['dem_prefix'] == prefix, attr] = value[attr]
    # if no match, then use a dummy link
    gdf.loc[gdf['Annotation'].isna(), 'Note'] = 'https://vims0-my.sharepoint.com/:o:/r/personal/feiye_vims_edu/_layouts/15/doc2.aspx?sourcedoc=%7Bd6aa1259-6a28-4b12-9caf-58f05600b3b5%7D&action=edit&wd=target(DEM.one%7CAE181E84-4618-4FAF-8AF2-32FE3A6C3D29%2F)onenote%3Ahttps%3A%2F%2Fvims0-my.sharepoint.com%2Fpersonal%2Ffeiye_vims_edu%2FDocuments%2FNotebooks%2FNOAA%20TWL%20project%E2%80%8B%2FDEM.one'

    gdf.to_file(dem_id_shp_fname)
    pass
