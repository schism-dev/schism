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

# ------------------------- inputs ---------------------------
wdir = '/sciclone/schism10/Hgrid_projects/SECOFS/new20_JZ/'
crs = 'esri:102008'
coastal_shpfile = f'{wdir}/coastal.shp'  # esri:102008
lbnd_coastal_shpfile = f'{wdir}/lbnd_coastal.shp'  # esri:102008
levee_shpfile = '/sciclone/schism10/feiye/SECOFS/Inputs/I02e/SMS/levee.shp'
select_nearshore_shpfile = '/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Grid/select_nearshore_v2.shp'  
# ------------------------- end inputs ---------------------------

coastal = gpd.read_file(coastal_shpfile)
coastal.set_crs(crs, inplace=True)

lbnd_coastal = gpd.read_file(lbnd_coastal_shpfile)
lbnd_coastal.set_crs(crs, inplace=True)
# raw_watershed is a rough estimate of watershed region, which excludes manual polygons but may include levee
raw_watershed = lbnd_coastal.overlay(coastal, how='difference').dissolve()

# extract the refined coastal polygons, i.e., coastal minus those intersecting the select_nearshore polygons, where no intersections with auto arcs are allowed, hence the buffer
select_nearshore = gpd.read_file(select_nearshore_shpfile)
nearshore_coastal = gpd.sjoin(coastal, select_nearshore, how="inner", predicate='intersects')
refined_coastal = coastal.overlay(nearshore_coastal, how='difference').dissolve()
refined_coastal_buf = gpd.GeoDataFrame(geometry=refined_coastal.buffer(50))
refined_coastal_buf.to_file(f'{wdir}/coastal_refined.shp', crs=crs)
# auto-arcs cannot intersect with the levee, hence the buffer
levee_buf = gpd.GeoDataFrame(geometry=gpd.read_file(levee_shpfile).buffer(50))
# get final watershed by subtracting polygons that cannot be intersected by auto arcs from the raw watershed
watershed = raw_watershed.overlay(refined_coastal_buf, how='difference').dissolve()
watershed = watershed.overlay(levee_buf, how='difference').dissolve()
watershed.to_file(f'{wdir}/watershed.shp', crs=crs)
watershed = gpd.read_file(f'{wdir}/watershed.shp')

# clip the auto arcs to the watershed boundary
total_arcs = gpd.read_file(f'{wdir}/total_arcs.shp').to_crs(crs)
total_arcs_clipped = total_arcs.clip(watershed)
total_arcs_clipped.to_file(f'{wdir}/total_arcs_clipped.shp', crs=crs)

# clip the polygons formed by auto arcs to the watershed boundary; this step may fail, in that case do it manually in qgis
total_arcs_polygons = gpd.read_file(f'{wdir}/total_river_polys.shp').to_crs(crs)
total_arcs_polygons_clipped = total_arcs_polygons.buffer(0).clip(watershed)
# put total_arcs_polygons_clipped back to gdf and dissolve
total_arcs_polygons_clipped = gpd.GeoDataFrame(geometry=total_arcs_polygons_clipped).dissolve()

total_arcs_polygons_clipped.to_file(f'{wdir}/total_river_polys_clipped_test.shp', crs=crs)

pass