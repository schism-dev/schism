#!/usr/bin/env python3
'''
This script clips the auto arcs to the watershed boundary,
which is the part of domain that is not ocean, levee, or manually refined coastal polygons.
A 50-m buffer is applied to the levee and refined coastal polygons to ensure that the auto arcs
does not directly intersect with the levee and refined coastal polygons.
However, no buffer is applied to the interface between watershed and ocean (i.e., shoreline)
, allowing intersections between the auto arcs and shoreline to provide channels where
there is no manual refinement.
'''

import os
import geopandas as gpd

# ------------------------- inputs ---------------------------
WDIR = '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v31/Clip/'
CRS = 'esri:102008'

# manual polygons defined in the coastal coverage,
# may need to delete some polygons so that they are treated as watershed,
# check coastal_remove coverage for polygons to be removed
coastal_shpfile = f'{WDIR}/inputs/coastal.shp'  # esri:102008

# land boundary + coastal, i.e., watershed region between the coastline and the land boundary
# as well as manual polygons in the coastal coverage, remove the land between CB and DB
lbnd_coastal_shpfile = f'{WDIR}/inputs/lbnd_coastal.shp'  # esri:102008

# use linestring as the levee shapefile
levee_shpfile = f'{WDIR}/inputs/levee.shp'  # esri:102008
LEVEE_BUF_DISTANCE = 5

# used to select un-refined polygons, which are allowed to be intersected by auto arcs
select_nearshore_shpfile = f'{WDIR}/inputs/select_nearshore_v4.shp'  # esri:102008

auto_arcs_file = f'{WDIR}/inputs/total_arcs.shp'  # lonlat
auto_polys_file = f'{WDIR}/inputs/total_river_polys.shp'  # lonlat
# ------------------------- end inputs ---------------------------

# make outputs directory
output_dir = f'{WDIR}/outputs'
os.makedirs(output_dir, exist_ok=True)
# save a copy of this script to the working directory
os.system(f'cp {__file__} {WDIR}')

coastal = gpd.read_file(coastal_shpfile)
coastal.set_crs(CRS, inplace=True)

lbnd_coastal = gpd.read_file(lbnd_coastal_shpfile)
lbnd_coastal.set_crs(CRS, inplace=True)

# raw_watershed is a rough estimate of watershed region,
# which excludes manual polygons but may include levee
# This includes most of the watershed region EXCEPT for the manual polygons in the coastal coverage
raw_watershed = lbnd_coastal.overlay(coastal, how='difference').dissolve()

# Extract the refined coastal polygons, i.e.,
# coastal minus the un-refined polyogns (those intersecting the select_nearshore),
# Reined polygons are manual polygons supposedly better configured than the auto arcs.
# It also implies that refined coastal polygons accomodate the connectivity to watershed
# rivers, e.g., a lake shoreline should have denser nodes where it connects to a river.
# However, this is not always the case because some manualy polygons are not accurate,
# No intersections between refined polygons and auto arcs are allowed, hence the buffer.
# Unrefined polygons are allowed to intersect with auto arcs to ensure channel connectivity
select_nearshore = gpd.read_file(select_nearshore_shpfile)
nearshore_coastal = gpd.sjoin(coastal, select_nearshore, how="inner", predicate='intersects')
refined_coastal = coastal.overlay(nearshore_coastal, how='difference').dissolve()
refined_coastal_buf = gpd.GeoDataFrame(geometry=refined_coastal.buffer(50))
refined_coastal_buf.to_file(f'{output_dir}/coastal_refined.shp')
print(f'coastal_refined.shp saved to {output_dir}')

# get final watershed
print('subtracting refined coastal from watershed')
watershed = raw_watershed.overlay(refined_coastal_buf, how='difference').dissolve()

print('subtracting levee from watershed')
levee_buf = gpd.GeoDataFrame(geometry=gpd.read_file(levee_shpfile).buffer(LEVEE_BUF_DISTANCE))
levee_buf.to_file(f'{output_dir}/levee_buf.shp')
watershed = watershed.overlay(levee_buf, how='difference').dissolve()
watershed.to_file(f'{output_dir}/watershed.shp')

# add additional watershed region (manually specified) to the watershed
print('break here if using qgis for clipping, which is faster')
watershed = gpd.read_file(f'{output_dir}/watershed.shp')

# clip the auto arcs to the watershed boundary
total_arcs = gpd.read_file(auto_arcs_file).to_crs(CRS)
total_arcs_clipped = total_arcs.clip(watershed)
total_arcs_clipped.to_file(f'{output_dir}/total_arcs_clipped.shp', crs=CRS)

# clip the polygons formed by auto arcs to the watershed boundary;
# this step may fail, in that case do it manually in qgis
total_arcs_polygons = gpd.read_file(auto_polys_file).to_crs(CRS)
total_arcs_polygons_clipped = total_arcs_polygons.buffer(0).clip(watershed)
# put total_arcs_polygons_clipped back to gdf and dissolve
total_arcs_polygons_clipped = gpd.GeoDataFrame(geometry=total_arcs_polygons_clipped).dissolve()

total_arcs_polygons_clipped.to_file(f'{output_dir}/total_river_polys_clipped_test.shp', crs=CRS)
print('done')
