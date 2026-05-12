import os
import numpy as np
from RiverMapper.SMS import SMS_MAP
import pickle
from pathlib import Path
import geopandas as gpd


arcs_map_file = Path('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/Outputs/CUDEM_merged_thalwegs_1e6_single_fix_simple_sms_cleaned_32cores/total_inner_arcs.map')

if os.path.exists(f'{arcs_map_file}.pkl'):
    with open(f'{arcs_map_file}.pkl', 'rb') as file:
        arcs_map = pickle.load(file)
else:
    arcs_map = SMS_MAP(arcs_map_file)
    # cache the inner arc map
    with open(f'{arcs_map_file}.pkl', 'wb') as file:
        pickle.dump(arcs_map, file)

arc_idx = 0
outer_arcs = []; inner_arcs = []; all_arcs = []
while arc_idx < arcs_map.n_arcs:
    n_row = int(arcs_map.arcs[arc_idx].points[0, 2])  # third column's integer part is the number of arcs in this group
    outer_arcs += [arcs_map.arcs[arc_idx], arcs_map.arcs[arc_idx+n_row-1]]  # first and last arc in this group
    inner_arcs += [arc for arc in arcs_map.arcs[arc_idx+1:arc_idx+n_row-1]]  # inner arcs in this group
    all_arcs.append([inner_arc for inner_arc in arcs_map.arcs[arc_idx:arc_idx+n_row]])  # all arcs in this group
    arc_idx += n_row

inner_arcs_map = SMS_MAP(arcs=inner_arcs)
inner_arcs_gdf = inner_arcs_map.to_GeoDataFrame()
# expand the lines to polygons by applying a buffer of 10 m
inner_arcs_gdf.geometry = inner_arcs_gdf.geometry.buffer(1e-4)
inner_arcs_gdf.to_file(f'{arcs_map_file.parent}/inner_arcs.shp', crs='epsg:4326')  

outer_arcs_map = SMS_MAP(arcs=outer_arcs)
# extract points from the outer arcs to nx3 np array
outer_arc_points = np.vstack([arc.points for arc in outer_arcs_map.arcs])
outer_arc_points_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(outer_arc_points[:, 0], outer_arc_points[:, 1]), crs='epsg:4326')
outer_arc_points_gdf.to_file(f'{arcs_map_file.parent}/outer_arcs_points.shp', crs='epsg:4326')

pass