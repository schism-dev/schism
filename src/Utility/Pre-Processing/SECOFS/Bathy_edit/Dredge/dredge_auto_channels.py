#!/usr/bin/env python3
import os
from pathlib import Path

#set the following environment variable to use pygeos
os.environ['USE_PYGEOS'] = '0'

import numpy as np
import geopandas as gpd

from pylib_essentials.schism_file import read_schism_hgrid_cached, grd2sms, schism_grid

def dredge_auto_channels(hgrid_obj:schism_grid, dredge_polygon_file, dredge_depth):
    '''
    This function is used to dredge the channel nodes automatically.
    :param hg_file: Path to the hgrid.gr3 file
    :param dredge_polygon_file: Path to the shapefile containing the polygons of the channel
    :param dredge_depth: dredge depth in meters
    :return:
    '''
    # load river polygons
    river_polys = gpd.read_file(dredge_polygon_file)  # epsg:4326, this is from clip_autoarcs.py
    # shrink the polygons by 1 m to exclude bank nodes
    river_polys.geometry = river_polys.geometry.to_crs('esri:102008').buffer(-1).to_crs('epsg:4326')

    # determine in-channel nodes
    hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(hgrid_obj.x, hgrid_obj.y), crs='epsg:4326')
    joined_gdf = gpd.sjoin(hg_points, river_polys, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()

    in_channel = np.zeros_like(hgrid_obj.dp, dtype=bool)
    in_channel[idx] = True
    hgrid_obj.plot(value=in_channel.astype(int), fmt=1)

    # dredge the in-channel nodes
    hgrid_obj.dp[idx] += dredge_depth

    return hgrid_obj

if __name__ == '__main__':
    # ----------- inputs -------------
    hg_file = Path('./hgrid_dem_levee_loaded_NCF_loaded_xGEOID20b_chart_loaded.gr3')
    dredge_depth = 2
    # ---------------------------------

    hgrid_obj= schism_grid(str(hg_file))  # epsg:4326

    hgrid_obj = dredge_auto_channels(hgrid_obj=hgrid_obj, dredge_polygon_file=Path('./total_river_polys_clipped_test.shp'), dredge_depth=dredge_depth)
    grd2sms(hgrid_obj, (f'{hg_file.parent}/{hg_file.stem}_dredged_{dredge_depth}m.2dm'))
    hgrid_obj.save(f'{hg_file.parent}/{hg_file.stem}_dredged_{dredge_depth}m.gr3', fmt=1)

    pass
