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
wdir = '/sciclone/schism10/Hgrid_projects/GoM/'
original = gpd.read_file(f'{wdir}/original.shp')
full_domain = gpd.read_file(f'{wdir}/full.shp')
islands = gpd.read_file(f'{wdir}/islands.shp')
buffer_distance = 50e-5   # set to 0 if no buffer is needed
# ------------------------- end inputs ---------------------------

# buffer the original domain so that the clipped arcs don't intersect with the original domain
original_buf = gpd.GeoDataFrame(geometry=original.buffer(buffer_distance))
# the extended region is the "watershed" region where the auto arcs are placed
extended_region = full_domain.overlay(original_buf, how='difference').dissolve()

# add any islands back
islands_buf = gpd.GeoDataFrame(geometry=islands.buffer(-buffer_distance))  # "minus" because we want the river arcs contained in the islands
extended_region = extended_region.append(islands_buf)

# clip the auto arcs to the extended region
total_arcs = gpd.read_file(f'{wdir}/total_arcs.shp')
total_arcs_clipped = total_arcs.clip(extended_region)
total_arcs_clipped.to_file(f'{wdir}/total_arcs_clipped.shp')

pass