from pylib_experimental.schism_file import cread_schism_hgrid
import copy
import geopandas as gpd
import numpy as np


script_dir = '/sciclone/home/feiye/stofs3d-atl/Pre_processing/Bathy_edit/Temporary_Fix_v7.2/'
reference_hgrid_file = "/sciclone/schism10/feiye/STOFS3D-v8/I15g_v7/Bathy_edit/hgrid_dem_edit.ll"
output_dir = '/sciclone/schism10/feiye/STOFS3D-v8/I15d4_v7/'  # Bathy_edit/Temporary_Fix_v7.2/'

gd_ll = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/I15e_v7/hgrid.gr3')

gd_meters = copy.deepcopy(gd_ll)
gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')
hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y), crs='esri:102008')

# deepen rivers
river_fix_dict = {
    'BayouLafourche': {
        'fname': f'{script_dir}/BayouLafourche_lonlat.shp',
        'depth': 2,
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

    print(f'Forcing minimum depth to {river_info["depth"]} meters for {river_name}.')
    gd_ll.dp[idx] = np.maximum(gd_ll.dp[idx], river_info['depth'])

# raise levees
levee_fix_dict = {
    'MTG_IWW_intersection': {
        'fname': f'{script_dir}/MTG_IWW_intersection.shp',  # lonlat
        'depth': -3.012,  # depth in meters
    }
}
for levee_name, levee_info in levee_fix_dict.items():
    levee_gdf = gpd.read_file(levee_info['fname'])
    levee_gdf.set_crs('epsg:4326', inplace=True)

    levee_gdf = levee_gdf.to_crs('esri:102008')

    joined_gdf = gpd.sjoin(hg_points, levee_gdf, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

    print(f'Forcing maximum depth to {levee_info["depth"]} meters for {levee_name}.')
    gd_ll.dp[idx] = np.minimum(gd_ll.dp[idx], levee_info['depth'])


# raise land surface
land_fix_dict = {
    'Lafourche_Parish': {
        'fname': f'{script_dir}/Lafourche_Parish.shp',  # lonlat
        'reference_hgrid': reference_hgrid_file
    }
}
for land_name, land_info in land_fix_dict.items():
    land_gdf = gpd.read_file(land_info['fname'])
    land_gdf = land_gdf.to_crs('esri:102008')

    joined_gdf = gpd.sjoin(hg_points, land_gdf, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

    original_dp = gd_ll.dp[idx]
    print(f'{sum(original_dp >= 1.0)} points with large depth in {land_name}.')
    reference_hgrid = cread_schism_hgrid(land_info['reference_hgrid'])
    print(f'Forcing minimum depth (maximum ground elevation) for {land_name}.')
    gd_ll.dp[idx] = np.minimum(gd_ll.dp[idx], reference_hgrid.dp[idx])


# Borrow depth from adjacent nodes
print('Borrowing depth from adjacent nodes.')
polygon_shp = gpd.read_file(f'{script_dir}/Lafourche_Parish_tweaks.shp')  # lonlat
polygon_shp = polygon_shp.to_crs('esri:102008')
for i, polygon in polygon_shp.iterrows():
    print('Processing polygon', i)
    polygon = gpd.GeoDataFrame(geometry=[polygon.geometry]).set_crs('esri:102008')
    idx = gpd.sjoin(hg_points, polygon, predicate='within', how='inner').index
    if len(idx) > 0:
        gd_ll.dp[idx] = np.min(gd_ll.dp[idx])  # highest point in the polygon
    else:
        raise ValueError(f'No points found in polygon {i}.')

# Save the modified grid
gd_ll.grd2sms(f'{output_dir}/temporary_fix_river.2dm')
gd_ll.save(f'{output_dir}/temporary_fix_river.gr3', fmt=1)
print('Done')
