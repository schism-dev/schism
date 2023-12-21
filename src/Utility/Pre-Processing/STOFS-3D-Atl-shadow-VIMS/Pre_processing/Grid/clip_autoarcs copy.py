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
from shapely.ops import polygonize

wdir = '/sciclone/schism10/feiye/Test/RUN02b_JZ/Dredge/'
crs = 'esri:102008'
'''
coastal = gpd.read_file(f'{wdir}/coastal.shp')
coastal.set_crs(crs, inplace=True)

lbnd_coastal = gpd.read_file(f'{wdir}/lbnd_coastal.shp')
lbnd_coastal.set_crs(crs, inplace=True)

raw_watershed = lbnd_coastal.overlay(coastal, how='difference').dissolve()

# make the refined coastal polygons, i.e., coastal minus those intersecting the select_nearshore polygons
select_nearshore = gpd.read_file(f'/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Grid/select_nearshore_v2.shp')
nearshore_coastal = gpd.sjoin(coastal, select_nearshore, how="inner", predicate='intersects')
refined_coastal = coastal.overlay(nearshore_coastal, how='difference').dissolve()
refined_coastal_buf = gpd.GeoDataFrame(geometry=refined_coastal.buffer(50))
refined_coastal_buf.to_file(f'{wdir}/coastal_refined.shp', crs=crs)

levee_buf = gpd.GeoDataFrame(geometry=gpd.read_file(f'{wdir}/levee.shp').buffer(50))

watershed = raw_watershed.overlay(refined_coastal_buf, how='difference').dissolve()
watershed = watershed.overlay(levee_buf, how='difference').dissolve()

watershed.to_file(f'{wdir}/watershed.shp', crs=crs)
'''

watershed = gpd.read_file(f'{wdir}/watershed.shp')

# total_arcs = gpd.read_file(f'{wdir}/total_arcs.shp').to_crs(crs)
# total_arcs_clipped = total_arcs.clip(watershed)
# total_arcs_clipped.to_file(f'{wdir}/total_arcs_clipped.shp', crs=crs)

total_arcs_polygons = gpd.read_file(f'{wdir}/local_river_polys.shp').to_crs(crs)
total_arcs_polygons_clipped = total_arcs_polygons.clip(watershed)
total_arcs_polygons_clipped.to_file(f'{wdir}/total_river_polys_clipped.shp', crs=crs)

pass