from pylib_essentials.schism_file import cread_schism_hgrid
import copy
import geopandas as gpd
import numpy as np
from stofs3d_setup.utils.projection import project_geodataframe, project_grid

gd_ll = cread_schism_hgrid('./14a.gr3')

gd_meters = copy.deepcopy(gd_ll)
project_grid(gd_meters, 'epsg:4326', 'esri:102008')
hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y), crs='esri:102008')

river_fix_dict = {
    'HoumaCanal': {
        'fname': 'HoumaCanal_lonlat.shp',
        'depth': 5,  # depth in meters
        'buffer': -3  # from banks
    },
    'BayouLafourche': {
        'fname': 'BayouLafourche_lonlat.shp',
        'depth': 1,  # depth in meters
        'buffer': -3  # from banks
    },
    'NewOrleansOutfallCanals': {
        'fname': 'NewOrleansOutfallCanals_lonlat.shp',
        'depth': 1,  # depth in meters
        'buffer': -5  # from levee centerline
    },
}

for river_name, river_info in river_fix_dict.items():
    river_gdf = gpd.read_file(river_info['fname'])
    river_gdf.set_crs('epsg:4326', inplace=True)

    river_gdf = project_geodataframe(
        river_gdf, 'esri:102008', f"{river_name} river polygon projection"
    )
    river_gdf.geometry = river_gdf.buffer(river_info['buffer'])

    joined_gdf = gpd.sjoin(hg_points, river_gdf, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

    gd_ll.dp[idx] = np.maximum(gd_ll.dp[idx], river_info['depth'])
    print(f'Forcing depth to {river_info["depth"]} meters for {river_name}.')

gd_ll.grd2sms('temporary_fix_river.2dm')
gd_ll.save('temporary_fix_river.gr3', fmt=1)
pass
