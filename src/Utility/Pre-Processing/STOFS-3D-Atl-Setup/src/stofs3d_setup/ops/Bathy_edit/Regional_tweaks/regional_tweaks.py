"""
Apply temporary regional tweaks to hgrid depth,
only in regions where the channel representation is not accurate.
"""

import numpy as np
import copy
from pylib import inside_polygon, read_schism_bpfile, grd2sms
try:
    # Try to import from pylib_experimental to speed up the grid reading
    from pylib_experimental.schism_file import cread_schism_hgrid as read_hgrid
    print("Using 'cread_schism_hgrid' from 'pylib_experimental.schism_file'")
except ImportError:
    # If that fails, fall back to pylib
    from pylib import read as read_hgrid
    print("Using 'read' from 'pylib'")


DEFAULT_REGIONAL_TWEAKS = {
    'min_5m_ll_noPR': 5,
    'SabinePass': 7,
    'BergenPoint': 5,
    'Washington_3': 15,
    'Elk_river': 2,
    'Hudson_river': 16,
    'James_river': 14,
    'NorthEast_river': 5,
    'Rappahannock_river': 6,
    'Susquehanna_river': 10,
    'York_river': 10,
    'Androscoggin_Kennebec_rivers': 3,
    'Merrimack_river': 3,
    'Patuxent_river': 5,
    'Penobscot_river': 5,
    'Saco_river': 3,
    'StCroix_river': 5,
    'Oyster_landing': 1,
    'st_lawrence1': 10,
    'st_lawrence2': 10,
    'st_lawrence3': 10,
}

DEFAULT_REGIONAL_TWEAKS_v7p2 = {  # incoorporating SECOFS updates, for SMS v27 and after
    'min_5m_ll_noPR': 5,
    'SabinePass': 7,
    'BergenPoint': 5,
     # 'Washington_3': 15,  # deleted
    'Elk_river': 2,
    'Hudson_river': 16,
    'James_river': 2,  # changed
    'NorthEast_river': 5,
    'Rappahannock_river': 6,
    'Susquehanna_river': 10,
    'York_river': 4,  # changed
    'Androscoggin_Kennebec_rivers': 3,
    'Merrimack_river': 3,
    'Patuxent_river': 5,
    'Penobscot_river': 5,
    'Saco_river': 3,
    'StCroix_river': 5,
    'Oyster_landing': 1,
    'st_lawrence1': 10,
    'st_lawrence2': 10,
    'st_lawrence3': 10,
}
DEFAULT_REGION_DIR = '/sciclone/schism10/Hgrid_projects/DEMs/regions/'


def tweak_hgrid_depth(hgrid, regional_tweaks=None, regions_dir=DEFAULT_REGION_DIR, min_max='min'):
    '''
    Set the minimum/maximum depth in the regions specified in regional_tweaks.
    '''
    if regional_tweaks is None:
        regional_tweaks = DEFAULT_REGIONAL_TWEAKS_v7p2

    for region, depth in regional_tweaks.items():
        region_file = f'{regions_dir}/{region}.reg'
        bp = read_schism_bpfile(region_file, fmt=1)
        in_region = inside_polygon(np.c_[hgrid.x, hgrid.y], bp.x, bp.y).astype(bool)

        if min_max == 'min':
            tweak_idx = in_region * (hgrid.dp < depth)
        elif min_max == 'max':
            tweak_idx = in_region * (hgrid.dp > depth)
        else:
            raise ValueError(f'Invalid min_max value: {min_max}')

        hgrid.dp[tweak_idx] = depth
        print(f'Applied min depth {depth} in region {region}')
    return hgrid


def shape_tweak(hgrid, gpkg_file, target_prj='esri:102008'):
    """
    Apply regional tweaks to hgrid depth based on a GeoPackage file
    containing the regions and their corresponding depth values.

    The GeoPackage file should have a layer with a geometry column and a 'value' attribute column.
    Value:
        999999.9 (deepen to max depth inside that polygon);
        smaller positive value (deepen to max(value, current depth) inside that polygon);
        -999999.9 (raise to min depth inside that polygon);
        larger negative value (raise ground to min(value, current depth) inside that polygon).

    Parameters:
    - hgrid: The hgrid object to be tweaked.
    - gpkg_file: Path to the GeoPackage file containing the regions and depth values.

    Returns:
    - The tweaked hgrid object.
    """
    import geopandas as gpd

    # reproject to meters
    hg = copy.deepcopy(hgrid)
    hgrid_xy_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(hg.x, hg.y), crs='EPSG:4326')

    # Read the GeoPackage file
    gdf = gpd.read_file(gpkg_file)

    hgrid_tweaked_idx = copy.deepcopy(hg)
    hgrid_tweaked_idx.dp[:] = 0

    for _, polygon in gdf.iterrows():
        value = polygon['value']  # Assuming there's a 'value' column for depth tweak
        idx = gpd.tools.sjoin(
            hgrid_xy_gdf, gpd.GeoDataFrame(geometry=[polygon.geometry], crs='EPSG:4326'),
            how='inner', predicate='within').index

        if len(idx) == 0:
            raise ValueError(
                f'No points found inside the polygon with value {value}. '
                'Check the geometry and coordinate reference system of the GeoPackage.')

        # Apply the depth tweak
        if value > 999999:  # Deepen to max depth
            hg.dp[idx] = np.max(hg.dp[idx])
            print(f'Deepened to max depth {np.max(hg.dp[idx])} in region {idx}\n')
        elif value < -999999:  # Raise to min depth
            hg.dp[idx] = np.min(hg.dp[idx])
            print(f'Raised to min depth {np.min(hg.dp[idx])} in region {idx}\n')
        elif value > 0:  # Deepen to max(value, current depth)
            hg.dp[idx] = np.maximum(hg.dp[idx], value)
            print(f'Deepened to max({value}, current depth) in region {idx}\n')
        elif value < 0:  # Raise to min(value, current depth)
            hg.dp[idx] = np.minimum(hg.dp[idx], value)
            print(f'Raised to min({value}, current depth) in region {idx}\n')
        else:
            raise ValueError(f'Invalid value in GeoPackage: {value}')

        # sanity check for raising deep points to shallow depth
        if value < 0 and any(hgrid.dp[idx] >= 1.0):
            print(f'Warning: Raised deep points to shallow depth in region {idx}, which may result in multiple vertical layers on shallow points.\n')
            original_deepest_idx = np.argmax(hgrid.dp[idx])
            original_shallowest_idx = np.argmin(hgrid.dp[idx])
            print(
                f'Original deepest point in this region: Point {idx[original_deepest_idx]} '
                f'with depth {hgrid.dp[idx[original_deepest_idx]]}'
                f' is raised to match the shallow depth of Point {idx[original_shallowest_idx]} '
                f'with depth {hgrid.dp[idx[original_shallowest_idx]]} '
            )

        hgrid_tweaked_idx.dp[idx] = 1  # Mark these points as tweaked

    # diagnostic outputs
    print(f'Total points tweaked: {hgrid_tweaked_idx.dp.sum()}')
    return hg, hgrid_tweaked_idx


def sample():
    '''
    Sample usage of the tweak_hgrid_depth function.
    '''
    wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I09b/Bathy_edit/DEM_loading/'

    hgrid = read_hgrid(f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3')
    hgrid = tweak_hgrid_depth(hgrid)
    hgrid.write(f'{wdir}/hgrid.tweaked.gr3')


def sample_shape_tweak():
    '''
    Sample usage of the shape_tweak function.
    '''
    wdir = '/sciclone/schism10/feiye/STOFS3D-v8/v7.2.1_static_inputs_2026_03_25/'
    gpkg_file = f'{wdir}/v7.2.1_fix.gpkg'

    hgrid = read_hgrid(f'{wdir}/v7.2.gr3')
    hgrid, hgrid_tweaked_idx = shape_tweak(hgrid, gpkg_file)
    grd2sms(hgrid_tweaked_idx, f'{wdir}/hgrid.tweaked_idx.2dm')
    grd2sms(hgrid, f'{wdir}/hgrid.shape_tweaked.2dm')
    hgrid.save(f'{wdir}/hgrid.shape_tweaked.gr3', fmt=1)


if __name__ == '__main__':
    sample_shape_tweak()
    print('Done')
