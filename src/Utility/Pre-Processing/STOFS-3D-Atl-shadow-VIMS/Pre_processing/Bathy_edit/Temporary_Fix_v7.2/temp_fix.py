from pylib_experimental.schism_file import cread_schism_hgrid
import copy
import geopandas as gpd
import numpy as np

gd_ll = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/I15d_v7/hgrid.ll')

gd_meters = copy.deepcopy(gd_ll)
gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')
hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y), crs='esri:102008')

# deepen rivers
river_fix_dict = {
    'BayouLafourche': {
        'fname': 'BayouLafourche_lonlat.shp',
        'depth': 0.99,  # depth in meters, ensure 1 layer (< 1 m)
        'buffer': -3  # from banks
    },
}

for river_name, river_info in river_fix_dict.items():
    river_gdf = gpd.read_file(river_info['fname'])
    river_gdf.set_crs('epsg:4326', inplace=True)

    river_gdf = river_gdf.to_crs('esri:102008')
    river_gdf.geometry = river_gdf.buffer(river_info['buffer'])

    joined_gdf = gpd.sjoin(hg_points, river_gdf, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

    gd_ll.dp[idx] = np.maximum(gd_ll.dp[idx], river_info['depth'])
    print(f'Forcing minimum depth to {river_info["depth"]} meters for {river_name}.')

# raise levees
levee_fix_dict = {
    'MTG_IWW_intersection': {
        'fname': 'MTG_IWW_intersection.shp',
        'depth': -3.012,  # depth in meters
    }
}
for levee_name, levee_info in levee_fix_dict.items():
    levee_gdf = gpd.read_file(levee_info['fname'])
    levee_gdf.set_crs('epsg:4326', inplace=True)

    levee_gdf = levee_gdf.to_crs('esri:102008')

    joined_gdf = gpd.sjoin(hg_points, levee_gdf, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

    gd_ll.dp[idx] = np.minimum(gd_ll.dp[idx], levee_info['depth'])
    print(f'Forcing maximum depth to {levee_info["depth"]} meters for {levee_name}.')

gd_ll.grd2sms('temporary_fix_river.2dm')
gd_ll.save('temporary_fix_river.gr3', fmt=1)
pass
