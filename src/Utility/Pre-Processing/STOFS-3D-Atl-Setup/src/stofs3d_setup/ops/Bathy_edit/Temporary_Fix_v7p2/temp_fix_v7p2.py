import copy
import geopandas as gpd
import numpy as np
from pylib import read
from stofs3d_setup.utils.projection import project_geodataframe, project_grid


def temp_fix_v7p2(gd_ll, wdir, reference_hgrid_file):
    '''
    A temporary fix for v7.2:
    1) Deepen Bayou Lafourche channel depths,
    2) Raise land surface to prevent leakage into Lafourche Parish based on a reference hgrid.
    3) Raise levees at the intersection of MTG and the IWW near Bayou Lafourche, effectively placing a flood gate

    gd_ll: hgrid object in lonlat, will be modified in place and returned
    '''

    gd_meters = copy.deepcopy(gd_ll)
    project_grid(gd_meters, 'epsg:4326', 'esri:102008')
    hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y), crs='esri:102008')

    # deepen rivers
    river_fix_dict = {
        'BayouLafourche': {
            'fname': f'{wdir}/BayouLafourche_lonlat.shp',
            'depth': 2,
            'buffer': -3  # from banks
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

        print(f'Forcing minimum depth to {river_info["depth"]} meters for {river_name}.')
        gd_ll.dp[idx] = np.maximum(gd_ll.dp[idx], river_info['depth'])

    # raise levees
    levee_fix_dict = {
        'MTG_IWW_intersection': {
            'fname': f'{wdir}/MTG_IWW_intersection.shp',  # lonlat
            'depth': -3.012,  # depth in meters
        }
    }
    for levee_name, levee_info in levee_fix_dict.items():
        levee_gdf = gpd.read_file(levee_info['fname'])
        levee_gdf.set_crs('epsg:4326', inplace=True)

        levee_gdf = project_geodataframe(
            levee_gdf, 'esri:102008', f"{levee_name} levee polygon projection"
        )

        joined_gdf = gpd.sjoin(hg_points, levee_gdf, how="inner", predicate='within')
        idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

        print(f'Forcing maximum depth to {levee_info["depth"]} meters for {levee_name}.')
        gd_ll.dp[idx] = np.minimum(gd_ll.dp[idx], levee_info['depth'])

    # raise land surface
    land_fix_dict = {
        'Lafourche_Parish': {
            'fname': f'{wdir}/Lafourche_Parish.shp',  # lonlat
            'reference_hgrid': reference_hgrid_file
        }
    }
    for land_name, land_info in land_fix_dict.items():
        land_gdf = gpd.read_file(land_info['fname'])
        land_gdf = project_geodataframe(
            land_gdf, 'esri:102008', f"{land_name} land polygon projection"
        )

        joined_gdf = gpd.sjoin(hg_points, land_gdf, how="inner", predicate='within')
        idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

        original_dp = gd_ll.dp[idx]
        print(f'{sum(original_dp >= 1.0)} points with large depth in {land_name}.')
        reference_hgrid = read(land_info['reference_hgrid'])
        print(f'Forcing minimum depth (maximum ground elevation) for {land_name}.')
        gd_ll.dp[idx] = np.minimum(gd_ll.dp[idx], reference_hgrid.dp[idx])

    # Borrow depth from adjacent nodes
    print('Borrowing depth from adjacent nodes.')
    polygon_shp = gpd.read_file(f'{wdir}/Lafourche_Parish_tweaks.shp')  # lonlat
    polygon_shp = project_geodataframe(
        polygon_shp, 'esri:102008', "temporary-fix tweak-polygon projection"
    )
    for i, polygon in polygon_shp.iterrows():
        print('Processing polygon', i)
        polygon = gpd.GeoDataFrame(geometry=[polygon.geometry]).set_crs('esri:102008')
        idx = gpd.sjoin(hg_points, polygon, predicate='within', how='inner').index
        if len(idx) > 0:
            gd_ll.dp[idx] = np.min(gd_ll.dp[idx])  # highest point in the polygon
        else:
            raise ValueError(f'No points found in polygon {i}.')

    # Save the modified grid
    return gd_ll


if __name__ == '__main__':
    wdir = '/sciclone/home/feiye/stofs3d-atl/Pre_processing/Bathy_edit/Temporary_Fix_v7.2/'
    reference_hgrid_file = "/sciclone/schism10/feiye/STOFS3D-v8/I15g_v7/Bathy_edit/hgrid_dem_edit.ll"
    gd_ll = read('/sciclone/schism10/feiye/STOFS3D-v8/I15e_v7/hgrid.gr3')
    gd_ll = temp_fix_v7p2(gd_ll=None, wdir=wdir, reference_hgrid_file=reference_hgrid_file)
