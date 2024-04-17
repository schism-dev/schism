#!/usr/bin/env python

from pylib_essentials.schism_file import read_schism_hgrid_cached
import geopandas as gpd
import numpy as np
from copy import deepcopy
from shapely.geometry import Point

gd1 = read_schism_hgrid_cached('Original/hgrid.ll.dem_loaded.gr3')
gd2 = read_schism_hgrid_cached('BlueTopo/hgrid.ll.dem_loaded.gr3')

region_gdf = gpd.read_file('v18_s2_v1_polys_dissolved.shp')
# set projection in place
region_gdf.set_crs('esri:102008', inplace=True)

# points_list = [Point(xy) for xy in np.c_[gd1.x, gd1.y]]
# points_gdf = gpd.GeoDataFrame(geometry=points_list, crs='epsg:4326')
points_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x=gd1.x, y=gd1.y), crs='epsg:4326')

points_in_region = gpd.sjoin(points_gdf, region_gdf.to_crs('epsg:4326'), how='inner', op='within').index.values
print(f'{len(points_in_region)} nodes in region.')

# get the max depth in the region
gd = deepcopy(gd1)
gd.dp[points_in_region] = np.maximum(gd1.dp[points_in_region], gd2.dp[points_in_region])
gd.save('hgrid_max_dp.gr3')

pass
