#!/usr/bin/env python3

"""
Script to set levee profiles in SCHISM hgrid files.
This script reads levee heights from the National Levee Database and additional
levee information, and sets the levee profiles in the hgrid object.

The levee top nodes are selected by buffering the levee centerline shapefile
by a specified distance.

The levee heights are set by matching the levee top nodes to the nearest levee height points
within 0.01 degree (roughly 1000 m).

Temporary measures:
* Temporary fixes are applied to include additional nodes as levee top nodes
  to prevent leaking of water through the levees, see temporary_levee_fix_102008.shp in
  centerline_shp_dict.
* Some levees in National Levee Database are not trustworthy, so they are set to
  the nearest valid levee heights
* Some levees are not in the National Levee Database, so they are set based on local levee shapefiles,
  see set_local_levee_profile for more details.
  Set dp to be the smallest of all local levee dp, i.e., highest levee, from MTG 2022 and 2024 levee heights.
* Some levees are not in NLD or local levee shapefiles, so they are set based on estimated additional levee information,
  see set_additional_dp for more details.
"""

import os
import tarfile
import copy

import numpy as np
import geopandas as gpd
from scipy import spatial

from spp_core.Download.download_nld import nld2map, nld2map2
try:
    from pylib_experimental.schism_file import cread_schism_hgrid as schism_read
except ImportError:
    from pylib import schism_grid as schism_read
from pylib import grd2sms


def set_levee_profile(gd=None, wdir='./', centerline_shp_dict=None, min_levee_height=6.56168):
    '''
    set levee profile based on the
    National Levee Database, and the centerline shapefile

    all top nodes of the levee will be attached to hgrid_obj as hgrid_obj.ilevee

    min_levee_height: 
      if levee height < 1 ft, it may not be true, so clip it to 6.56168 ft = 2 m
      note that this is just a rough estimate, the resulting min levee height will be 1 ft,
      because only lower levees are adjusted
    '''
    if centerline_shp_dict is None:
        raise ValueError('centerline_shp_dict is not provided')

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = "./Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()

    levee_names = ['FEMA_region_levees']
    levee_name_str = "_".join(levee_names)
    levee_xyz = np.zeros((0, 3), dtype=float)

    for _ in levee_names:
        # read levee heights as xyz
        # _, xyz = nld2map(nld_fname=f'{wdir}/System.geojson')
        _, xyz = nld2map2(nld_fname=f'{wdir}/SystemRoute.geojson')
        # assert np.allclose(xyz, xyz2), "nld2map and nld2map2 should return the same xyz"
        levee_xyz = np.r_[levee_xyz, xyz]
    levee_x = levee_xyz[:, 0]
    levee_y = levee_xyz[:, 1]
    levee_height = levee_xyz[:, 2]

    # deal with small levee heights
    invalid_mask = levee_height < 1  # some levee heights are very close to 0, e.g., 0.000130833403 ft
    valid_levee_height = levee_height[~invalid_mask]
    _, idx = spatial.cKDTree(np.c_[levee_x[~invalid_mask], levee_y[~invalid_mask]]).query(
        np.c_[levee_x[invalid_mask], levee_y[invalid_mask]], k=1)
    # set the invalid levee heights to the nearest valid levee height
    levee_height[invalid_mask] = valid_levee_height[idx]
    # and clip at 6.56168 ft = 2 m, a minimum levee height for invalid levee heights
    levee_height[invalid_mask] = np.clip(levee_height[invalid_mask], min_levee_height, None)  # clip at 2 m
    np.savetxt(f'{wdir}/{levee_name_str}_clipped_{min_levee_height}ft.xyz', np.c_[levee_x, levee_y, levee_height], fmt='%.10f')

    # convert levee heights to meters
    levee_height *= 0.3048

    # plt.plot(np.sort(levee_height))
    # plt.show()

    if gd is None:
        gd = schism_read(f'{wdir}/hgrid.ll')  # ; gd.save(f'{wdir}/hgrid.pkl')

    gd.lon = gd.x
    gd.lat = gd.y
    gd.proj(prj0='epsg:4326', prj1='esri:102008')  # this overwrites gd.x, gd.y

    # find levee center line points in hgrid, use projected cooridnates to avoid truncation error
    # initialize: all false (no levee)
    if not hasattr(gd, 'ilevee'):  # record the top nodes of the levee in the hgrid object
        gd.ilevee = np.zeros_like(gd.dp, dtype=bool)
    for _, levee_info in centerline_shp_dict.items():
        levee_centerline_gdf = gpd.read_file(levee_info['fname'])
        levee_centerline_gdf.geometry = levee_centerline_gdf.geometry.buffer(levee_info['buffer'])
        # gd.x, gd.y should already be esri:102008
        hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd.x, gd.y), crs='esri:102008')
        # determine which points are inside the levee buffer polygons
        joined_gdf = gpd.sjoin(hg_points, levee_centerline_gdf, how="inner", predicate='within')
        idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons
        gd.ilevee[idx] = True
    gd.save(f'{wdir}/{levee_name_str}.gr3', value=gd.ilevee.astype(int))

    dist, idx = spatial.cKDTree(np.c_[levee_x, levee_y]).query(np.c_[gd.lon[gd.ilevee], gd.lat[gd.ilevee]])
    short_dist = dist < 0.01  # degree, roughly 1000 m
    idx_levee_in_range = np.argwhere(gd.ilevee)[:, 0][short_dist]
    gd.dp[idx_levee_in_range] = - levee_height.astype(float)[idx][short_dist]

    gd.x = gd.lon
    gd.y = gd.lat

    return gd  # levee loaded hgrid.ll


def set_additional_dp(gd_ll=None, additional_levee_info=None):
    '''additional levee info should specify shapefile name and forced dp for each levee'''

    if additional_levee_info is None:
        raise ValueError('additional_levee_info is not provided')

    # Check levee info existence
    levee_info_dir = './Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = "./Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall('./')
        my_tar.close()

    if gd_ll is None:
        raise ValueError('gd_ll is not provided')
    if additional_levee_info == {}:
        print('additional_levee_info is empty, no additional levee dp is set and original hgrid is returned')
        return gd_ll

    gd_meters = copy.deepcopy(gd_ll)
    gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')

    for levee_name, levee_info in additional_levee_info.items():
        print(f"forcing dp <= {levee_info['dp']} for {levee_name}")
        levee_centerline_gdf = gpd.read_file(levee_info['shapefile'])
        # buffer 3 m
        levee_centerline_gdf.geometry = levee_centerline_gdf.geometry.buffer(3)
        hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y), crs='esri:102008')
        # determine which points are inside the levee buffer polygons
        joined_gdf = gpd.sjoin(hg_points, levee_centerline_gdf, how="inner", predicate='within')
        i_inpoly = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

        gd_ll.dp[i_inpoly] = np.minimum(gd_ll.dp[i_inpoly], levee_info['dp'])

    return gd_ll


def set_local_levee_profile(gd_ll=None, local_levee_info=None, i_levee_top_node=None, max_levee=False):
    '''
    This function is to set the levee profile for the local levee, which is not in the National Levee Database,
    Note that the meshing is mostly based on centerlines from National Levee Database, but the local levee shapefile
    may have slight differences from the centerlines.
    To simplify the searching process, local levee nodes are selected from found levee top nodes that are
    within a buffer (depending on local_levee_info) of the local levee centerlines.

    local_levee_info is a dictionary with the following keys:
    shapefile: the shapefile of the local levee centerline
    height_points: an nx3 array, where n is the number of height points, and the columns are lon, lat, z
    buffer: the buffer distance in meters to select the local levee nodes
    For example:
    local_levee_info = {
        'MTG': {
            'shapefile': f'Levee_info/Polygons/MTG_Centerline_102008.shp',
            'height_points': local_levee_heights,
            'buffer': 2000
        }
    }

    i_levee_top_node is a boolean array indicating the top nodes of the levees,
    which can be obtained from gd_ll.ilevee after setting levees from National Levee Database
    by set_levee_profile.

    max_levee: if True, set dp to be the smallest of all local levee dp, i.e., highest levee
                if False, overwrite the dp with the local levee dp
    '''

    if gd_ll is None:
        raise ValueError('gd_ll is not provided')
    if local_levee_info is None:
        print('local_levee_info is empty, no local levee dp is set and original hgrid is returned')
        return gd_ll

    gd_meters = copy.deepcopy(gd_ll)
    gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')

    for i, [levee_name, levee_info] in enumerate(local_levee_info.items()):
        print(f"forcing local levee dp for {levee_name}")
        levee_centerline_gdf = gpd.read_file(levee_info['shapefile'])
        levee_centerline_gdf.geometry = levee_centerline_gdf.geometry.buffer(levee_info['buffer'])
        hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y), crs='esri:102008')
        joined_gdf = gpd.sjoin(hg_points, levee_centerline_gdf, how="inner", predicate='within')
        in_poly_idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons
        i_inpoly = np.zeros_like(gd_ll.dp, dtype=bool)  # initialize: all false (no local levee)
        i_inpoly[in_poly_idx] = True
        # only take the top node of the levee
        i_top_local_levee = i_inpoly & i_levee_top_node

        # add local levee top nodes to gd_ll.ilevee
        if not hasattr(gd_ll, 'ilevee'):
            gd_ll.ilevee = i_top_local_levee
        else:
            gd_ll.ilevee = np.logical_or(gd_ll.ilevee, i_top_local_levee)

        # find nearest height point for each top node of the local levee
        _, node2heightpoints = spatial.cKDTree(levee_info['height_points'][:, :2]).query(
            np.c_[gd_ll.lon[i_top_local_levee], gd_ll.lat[i_top_local_levee]])
        local_levee_dp = - levee_info['height_points'][node2heightpoints, 2]  # invert the z value to get dp

        if max_levee:
            # set dp to be the smallest of all local levee dp, i.e., highest levee
            if i == 0:
                gd_ll.dp[i_top_local_levee] = local_levee_dp
            else:
                gd_ll.dp[i_top_local_levee] = np.minimum(gd_ll.dp[i_top_local_levee], local_levee_dp)
        else:  # overwrite local levee dp in the order of local_levee_info
            gd_ll.dp[i_top_local_levee] = local_levee_dp
            
        # gd_ll.save(f'./hgrid_local_levee_loaded_ll.gr3')
        # gd_ll.save(f'./hgrid_local_levee_{levee_name}.gr3', value=i_top_local_levee.astype(int))

    return gd_ll


def set_levees(hgrid_obj, wdir):
    """
    A temporary levee setting script for the MTG levees.
    The levee heights are in NAVD88, but the hgrid is in xGEOID20b,
    so the conversion difference (output from the xGEOID20b conversion step)
    is utilized to save time.

    Levee heights are loaded in 3 steps:
    1) National Levee Database
    2) Additional manually tweaked levee heights
    3) Local levee heights (not in the National Levee Database)
        This part may need to be adjusted based on the local levee shapefile
    """

    # levee nodes that take values from National Levee Database
    centerline_shp_dict = {
        # buffer 7 m, because the foot of the levee is about 15 m wide,
        # so levee tops are selected but not the levee base
        'main': {'fname': 'Levee_info/Polygons/levee_v3p2_centerline_102008.shp', 'buffer': 7},
        # buffer 100 m, to select tweaked levee nodes (by improve_hgrid.py), which are mostly at intersections
        'temporary': {'fname': 'Levee_info/Polygons/temporary_levee_fix_102008.shp', 'buffer': 100},
    }

    # additional tweaks on levee heights
    additional_levee_info = {
        # 'Bonnet Carre Spill Way': {
        #     'dp': -9,
        #     'shapefile': 'Levee_info/Polygons/v3_BonnetCarre_centerline_esri102008.shp',
        # },
        'Herbert Hoover Dam': {
            'dp': -9,
            'shapefile': 'Levee_info/Polygons/v3_HH_Dam_centerline_esri102008.shp',
        },
        'Additional Mississippi River Levee': {
            'dp': -9,
            'shapefile': 'Levee_info/Polygons/v3_additional_levee_centerline_esri102008.shp',
        },
        'Upstream Mississippi River Levee': {
            'dp': -18.2,  # 60 ft, this is the levee height at the most downstream point of the shapefile
            'shapefile': 'Levee_info/Polygons/v3_la_upstream_missi_centerline_esri102008.shp',
        }
    }

    # local levees
    # MTG 2022
    local_levee_heights_gdf = gpd.read_file(f'{wdir}/Levee_info/Polygons/2022_MTG_CL_Survey_NAVD88_lonlat_feet.shp')
    elev = local_levee_heights_gdf['Field4'].to_numpy() * 0.3048  # convert feet to meters
    local_levee_heights_2022 = np.c_[local_levee_heights_gdf.geometry.x, local_levee_heights_gdf.geometry.y, elev]

    # MTG 2024, which is a subset of the 2022 levee but higher
    local_levee_heights_gdf = gpd.read_file(
        f'{wdir}/Levee_info/Polygons/2024-05-22_MTG Levee Elevations_NAVD88_lonlat_feet.shp')
    elev = local_levee_heights_gdf['field_4'].to_numpy() * 0.3048  # convert feet to meters
    local_levee_heights_2024 = np.c_[local_levee_heights_gdf.geometry.x, local_levee_heights_gdf.geometry.y, elev]

    # MTG 2025
    local_levee_heights_gdf = gpd.read_file(
        f'{wdir}/Levee_info/Polygons/2025_CL_Profile_Data_DCC_25-2_lonlat_ft.shp')  # NAVD88
    elev = local_levee_heights_gdf['Elev (FT)'].to_numpy() * 0.3048  # convert feet to meters
    local_levee_heights_2025 = np.c_[local_levee_heights_gdf.geometry.x, local_levee_heights_gdf.geometry.y, elev]

    local_levee_info = {
        'MTG_2022': {
            'shapefile': 'Levee_info/Polygons/2022_MTG_Centerline_102008.shp',
            'height_points': local_levee_heights_2022,
            'buffer': 2000  # buffer 2000 m, a wide buffer is okay because only the levee top nodes will be selected
        },
        'MTG_2024': {
            'shapefile': 'Levee_info/Polygons/2024_MTG_Centerline_102008.shp',
            'height_points': local_levee_heights_2024,
            'buffer': 2000  # buffer 2000 m, a wide buffer is okay because only the levee top nodes will be selected
        },
        'MTG_2025': {
            'shapefile': 'Levee_info/Polygons/2025_MTG_Centerline_102008.shp',
            'height_points': local_levee_heights_2025,
            'buffer': 250  # buffer 250 m, a wide buffer is okay because only the levee top nodes will be selected
        }
    }
    # ---------------------------------------------

    print('loading levee heights from National Levee Database')
    # all top nodes of the levee will be attached to hgrid_obj as hgrid_obj.ilevee
    # this also writes ilevee (levee top nodes) in {levee_name}.gr3 for diagnostic purposes
    hgrid_obj = set_levee_profile(gd=hgrid_obj, wdir=wdir, centerline_shp_dict=centerline_shp_dict)
    grd2sms(hgrid_obj, f'{wdir}/hgrid_with_NLD_levees.2dm')

    # print('force minimum dp to be above -7 m for all levee top points')
    # hgrid_obj.dp[hgrid_obj.ilevee] = np.minimum(hgrid_obj.dp[hgrid_obj.ilevee], -7)

    print('loading additional tweaks on levee heights')
    hgrid_obj = set_additional_dp(gd_ll=hgrid_obj, additional_levee_info=additional_levee_info)

    print('loading local levee heights')
    hgrid_obj = set_local_levee_profile(
        gd_ll=hgrid_obj, local_levee_info=local_levee_info, i_levee_top_node=hgrid_obj.ilevee)

    grd2sms(hgrid_obj, f'{wdir}/hgrid_with_all_levees.2dm')
    return hgrid_obj


def v7p3():
    """
    Load 2025 NLD and MTG levee heights into hgrid,
    which updates v7.2 to v7.3
    """

    # prepare base hgrid to load new levee heights
    WDIR = '/sciclone/schism10/feiye/STOFS3D-v8/I15y_v7/Set_Levee/'
    os.chdir(WDIR)  # to set the directory
    hg = schism_read(f'{WDIR}/v7.2_2025_04_12.gr3')
    hg_original = copy.deepcopy(hg)
    
    # set levee heights
    hg.dp[:] = -9999  # initialize dp to -9999
    hg.save(f'{WDIR}/v7.2_2025_04_12_dp=-9999.gr3')

    set_levees(hgrid_obj=hg, wdir=WDIR)
    hg.save(f'{WDIR}/hgrid_with_levees.gr3')

    valid_value_mask = hg.dp > -9000
    print(f"valid values after setting levees: {np.sum(valid_value_mask)} of {hg.np}")

    # datum shift only for valid values
    xgeoid_minus_navd = schism_read('/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/Bathy_edit/xGEOID/hgrid_xGEOID20b_dp_minus_NAVD_dp.gr3')
    hg.dp[valid_value_mask] += xgeoid_minus_navd.dp[valid_value_mask]

    # plug in original dp for non-levee nodes
    hg.dp[~valid_value_mask] = hg_original.dp[~valid_value_mask]
    hg.save(f'{WDIR}/hgrid_with_levees_xGEOID.gr3', fmt=1)


if __name__ == '__main__':
    v7p3()
    print('Done')
