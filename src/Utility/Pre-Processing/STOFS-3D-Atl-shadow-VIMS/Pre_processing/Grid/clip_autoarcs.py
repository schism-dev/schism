'''
This script clips the auto arcs to the watershed boundary,
which is the part of domain that is not ocean, levee, or manually refined coastal polygons.
A 50-m buffer is applied to the levee and refined coastal polygons to ensure that the auto arcs
does not directly intersect with the levee and refined coastal polygons.
However, no buffer is applied to the interface between watershed and ocean (i.e., shoreline)
, allowing intersections between the auto arcs and shoreline to provide channels where
there is no manual refinement.
'''

import geopandas as gpd

wdir = '/sciclone/schism10/Hgrid_projects/STOFS3D-v7/new19/auto_arc_process/'
crs = 'esri:102008'

coastal = gpd.read_file(f'{wdir}/coastal.shp', crs=crs)
lbnd_coastal = gpd.read_file(f'{wdir}/lbnd_coast.shp', crs=crs)

raw_watershed = lbnd_coastal.overlay(coastal, how='difference').dissolve()

coastal_refined_buf = gpd.GeoDataFrame(geometry=gpd.read_file(f'{wdir}/coastal_refined.shp', crs=crs).buffer(50))
levee_buf = gpd.GeoDataFrame(geometry=gpd.read_file(f'{wdir}/levee.shp', crs=crs).buffer(50))

watershed = raw_watershed.overlay(coastal_refined_buf, how='difference').dissolve()
watershed = watershed.overlay(levee_buf, how='difference').dissolve()

watershed.to_file(f'{wdir}/watershed.shp', crs=crs)
watershed = gpd.read_file(f'{wdir}/watershed.shp', crs=crs)

total_arcs = gpd.read_file(f'{wdir}/total_arcs.shp').to_crs(crs)
total_arcs_clipped = total_arcs.clip(watershed)
total_arcs_clipped.to_file(f'{wdir}/total_arcs_clipped.shp', crs=crs)

pass