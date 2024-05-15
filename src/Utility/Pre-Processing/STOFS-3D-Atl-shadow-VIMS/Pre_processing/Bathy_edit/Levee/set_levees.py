#!/usr/bin/env python3

import numpy as np
from schism_py_pre_post import Datafiles
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid, grd2sms
from schism_py_pre_post.Download.download_nld import nld2map
from pathlib import Path
import geopandas as gpd
import copy
import os
from scipy import spatial
import tarfile


def set_levee_profile(gd=None, wdir='./', centerline_shp_list=[]):

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()

    levee_names = ['FEMA_region_levees']
    levee_name_str = "_".join(levee_names)
    levee_xyz = np.zeros((0, 3), dtype=float)

    for levee_name in levee_names:
        # read levee heights as xyz
        _, xyz = nld2map(nld_fname=f'/sciclone/schism10/Hgrid_projects/Levees/Levee_v3/FEMA_regions/FEMA_region_levees/System.geojson')
        levee_xyz = np.r_[levee_xyz, xyz]
    levee_x = levee_xyz[:, 0]
    levee_y = levee_xyz[:, 1]
    levee_height = levee_xyz[:, 2]
    levee_height[levee_height < 1] = 27  # raise very low levees to 9 meters, because their heights are not trustworthy
    levee_height *= 0.3048  # convert to meters
    # plt.plot(np.sort(levee_height))
    # plt.show()

    if gd is None:
        gd = schism_grid(f'{wdir}/hgrid.ll')  # ; gd.save(f'{wdir}/hgrid.pkl')

    gd.lon = gd.x
    gd.lat = gd.y
    gd.proj(prj0='epsg:4326', prj1='esri:102008')  # this overwrites gd.x, gd.y

    # find levee center line points in hgrid, use UTM to avoid truncation error
    shapefile_names = centerline_shp_list
    # initialize: all false (no levee)
    ilevee = np.zeros_like(gd.dp, dtype=bool)
    for shapefile_name in shapefile_names:
        levee_centerline_gdf = gpd.read_file(shapefile_name)
        # buffer 3 m
        levee_centerline_gdf.geometry = levee_centerline_gdf.geometry.buffer(3)
        hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd.x, gd.y), crs='esri:102008')  # gd.x, gd.y are already converted to esri:102008
        # determine which points are inside the levee buffer polygons
        joined_gdf = gpd.sjoin(hg_points, levee_centerline_gdf, how="inner", predicate='within')
        idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons
        ilevee[idx] = True

    gd.save(f'{wdir}/{levee_name_str}.gr3', value=ilevee)

    II = spatial.cKDTree(np.c_[levee_x, levee_y]).query(np.c_[gd.lon[ilevee], gd.lat[ilevee]])[1]
    dist = np.sqrt((gd.lon[ilevee] - levee_x[II])**2 + (gd.lat[ilevee] - levee_y[II])**2)
    short_dist = dist < 0.01  # degree, roughly 1000 m

    # gd.dp[:] = 0
    idx_levee_in_range = np.argwhere(ilevee)[:, 0][short_dist]
    gd.dp[idx_levee_in_range] = - levee_height.astype(float)[II][short_dist]

    gd.x = gd.lon
    gd.y = gd.lat

    # os.system(f"cp {wdir}/hgrid_{levee_name_str}_loaded_ll.gr3 {wdir}/hgrid.ll")
    # proj(
    #     f'{wdir}/hgrid.ll', 0, 'epsg:4326',
    #     f'{wdir}/hgrid.utm.gr3', 0, 'epsg:26918',
    # )

    return gd  # levee loaded hgrid.ll

def set_additional_dp(gd_ll=None, additional_levee_info:dict={}):
    '''additional levee info should specify shapefile name and forced dp for each levee'''

    # Check levee info existence
    levee_info_dir = f'.//Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
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
        print(f"forcing dp = {levee_info['dp']} for {levee_name}")
        levee_centerline_gdf = gpd.read_file(levee_info['shapefile'])
        # buffer 3 m
        levee_centerline_gdf.geometry = levee_centerline_gdf.geometry.buffer(3)
        hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y), crs='esri:102008')
        # determine which points are inside the levee buffer polygons
        joined_gdf = gpd.sjoin(hg_points, levee_centerline_gdf, how="inner", predicate='within')
        i_inpoly = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons

        gd_ll.dp[i_inpoly] = np.minimum(gd_ll.dp[i_inpoly], levee_info['dp'])

    return gd_ll

def set_levees(hgrid_obj:schism_grid, wdir):

    # additional tweaks on levee heights
    additional_levee_info = {
        'Bonnet Carre Spill Way': {
            'dp': -9,
            'shapefile': f'Levee_info/Polygons/v3_BonnetCarre_centerline_esri102008.shp',
        },
        'Herbert Hoover Dam': {
            'dp': -9,
            'shapefile': f'Levee_info/Polygons/v3_HH_Dam_centerline_esri102008.shp',
        },
        'Additional Mississippi River Levee': {
            'dp': -9,
            'shapefile': f'Levee_info/Polygons/v3_additional_levee_centerline_esri102008.shp',
        },
        'Upstream Mississippi River Levee': {
            'dp': -25,
            'shapefile': f'Levee_info/Polygons/v3_la_upstream_missi_centerline_esri102008.shp',
        }
    }
    # ---------------------------------------------

    print('loading levee heights from National Levee Database')
    centerline_shp_list = ['Levee_info/Polygons/levee_v3p2_centerline_102008.shp']
    hgrid_obj = set_levee_profile(gd=hgrid_obj, wdir=wdir, centerline_shp_list=centerline_shp_list)

    print('loading additional tweaks on levee heights')  # some heights needs to be reverted
    hgrid_obj = set_additional_dp(gd_ll=hgrid_obj, additional_levee_info=additional_levee_info)

    return hgrid_obj
