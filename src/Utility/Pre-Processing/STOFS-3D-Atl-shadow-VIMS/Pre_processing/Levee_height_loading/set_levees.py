import numpy as np
from schism_py_pre_post import Datafiles
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid, inside_polygon
from schism_py_pre_post.Download.download_nld import nld2map
import copy
import os
from scipy import spatial
import tarfile
import shapefile
import shutil


def set_constant_levee_height(gd=None, wdir='./'):
    # set constant levee heights for levees not in NLD.

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()

    if gd is None:
        gd = schism_grid(f'{wdir}/hgrid.ll')

    gd_meters = copy.deepcopy(gd)
    gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')

    for set_depth, shapefile in zip(
        [-9, -9], [f'{levee_info_dir}/Polygons/additional_levee_poly_102008.shp', f"{levee_info_dir}/Polygons/HH_Dam_buffer_13m_102008.shp"]
    ):
        i_inpoly = find_node_in_shpfiles(shapefile_names=[shapefile], gd=gd_meters)
        gd.dp[i_inpoly] = np.minimum(gd.dp[i_inpoly], set_depth)

    return gd

def set_levee_profile(gd=None, wdir='./'):

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()

    levee_names = ['LA_levees', 'FL_levees']
    levee_name_str = "_".join(levee_names)
    levee_xyz = np.zeros((0, 3), dtype=float)

    for levee_name in levee_names:
        # read levee heights as xyz
        _, xyz = nld2map(nld_fname=f'{levee_info_dir}/{levee_name}/System.geojson')
        levee_xyz = np.r_[levee_xyz, xyz]
    levee_x = levee_xyz[:, 0]
    levee_y = levee_xyz[:, 1]
    levee_height = levee_xyz[:, 2]
    levee_height[levee_height < 1] = 27  # raise very low levees to 9 meters
    levee_height *= 0.3048  # convert to meters
    # plt.plot(np.sort(levee_height))
    # plt.show()

    if gd is None:
        gd = schism_grid(f'{wdir}/hgrid.ll')  # ; gd.save(f'{wdir}/hgrid.pkl')

    gd.lon = gd.x
    gd.lat = gd.y
    gd.proj(prj0='epsg:4326', prj1='esri:102008')  # this overwrites gd.x, gd.y

    # find levee center line points in hgrid, use UTM to avoid truncation error
    shapefile_names = [
        f"{levee_info_dir}/Polygons/la_levee_center_line_buffer_102008.shp",
        f"{levee_info_dir}/Polygons/fl_levees_buffer_10m_102008.shp",
    ]
    ilevee = np.zeros(gd.dp.shape)
    for shapefile_name in shapefile_names:
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()
        for i, shp in enumerate(shapes):
            print(f'shp {i+1} of {len(shapes)}')
            poly_xy = np.array(shp.points).T
            ilevee += inside_polygon(np.c_[gd.x, gd.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
    ilevee = ilevee.astype('bool')

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


def set_additional_dp_v11_91(gd_ll=None, gd_dem=None, wdir='./'):
    # v11.91, set upstream mississippi river levee height to 25 m, to prevent overtopping near the source
    #         set Bonnet Carre Spill Way gate to 8.5 m
    #         revert levee height near Ostrica, LA to DEM values

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()


    if gd_ll is None:
        gd_ll = schism_grid(f'{wdir}/hgrid.ll')
    if gd_dem is None:
        gd_dem = schism_grid(f'{wdir}/hgrid.DEM_loaded.ll')

    gd_meters = copy.deepcopy(gd_ll)
    gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')

    for set_depth, shapefile in zip(
        [-9, -20],
        [f'{levee_info_dir}/Additional_Polygons/BonnetCarre_102008.shp',
         f'{levee_info_dir}/Additional_Polygons/la_levee_center_line_upstream_missi_13m_buffer_102008.shp']
    ):
        i_inpoly = find_node_in_shpfiles(shapefile_names=[shapefile], gd=gd_meters)
        gd_ll.dp[i_inpoly] = np.minimum(gd_ll.dp[i_inpoly], set_depth)

    i_inpoly = find_node_in_shpfiles(shapefile_names=[f'{levee_info_dir}/Additional_Polygons/la_levee_center_line_Ostrica_buffer_102008.shp'], gd=gd_meters)
    gd_ll.dp[i_inpoly] = gd_dem.dp[i_inpoly]


    # os.system(f"cp {wdir}/hgrid.additional_dp.ll {wdir}/hgrid.ll")
    # proj(
    #     f'{wdir}/hgrid.ll', 0, 'epsg:4326',
    #     f'{wdir}/hgrid.utm.gr3', 0, 'epsg:26918',
    # )
    # os.system(f"cp {wdir}/hgrid.utm.gr3 {wdir}/hgrid.utm.26918")

    return gd_ll

def tweak_depths(hgrid_name='', gd:schism_grid=None, original_hgrid_name='hgrid.ll'):
    shutil.move(hgrid_name, hgrid_name+'.gr3')
    hgrid_name += '.gr3'
    if gd is None:
        gd = schism_grid(hgrid_name)
    else:
        hgrid_name = gd.source_file

    dirname = os.path.dirname(hgrid_name)

    gd_ll_original = schism_grid(f'{dirname}/{original_hgrid_name}')  # before loading DEM (x, y may be slightly changed after DEM loading)


    os.system(f'mv {hgrid_name} {dirname}/hgrid.ll')
    gd = schism_grid(f'{dirname}/hgrid.ll')
    gd.x, gd.y = gd_ll_original.x, gd_ll_original.y  # force original x, y
    gd_DEM_loaded = copy.deepcopy(gd)

    print('set default levee heights (-9 m)')
    gd = set_constant_levee_height(gd=gd, wdir=dirname)

    print('loading levee heights from National Levee Database')
    gd = set_levee_profile(gd=gd, wdir=dirname)

    print('loading additional tweaks on levee heights')
    gd = set_additional_dp_v11_91(gd_ll=gd, gd_dem=gd_DEM_loaded, wdir=dirname)

    print('outputing hgrid.ll')
    gd.save(f'{dirname}/hgrid.ll')

    # set feeder dp
    # gd = set_feeder_dp(
    #     feeder_info_dir='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/',
    #     new_grid_dir=dirname
    # )
    # os.system(f'mv {dirname}/hgrid.ll {dirname}/hgrid.ll_before_feeder_dp')
    # gd.save(f'{dirname}/hgrid.ll')

    return gd


if __name__ == "__main__":
    # Sample usage
    wdir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/v14.42_post_proc2/'
    gd = tweak_depths(f'{wdir}/hgrid.ll.new')
