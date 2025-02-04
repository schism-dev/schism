import numpy as np
import os
from pathlib import Path

# pip install git+https://github.com/feiye-vims/schism_py_pre_post.git
from schism_py_pre_post.Shared_modules.hotstart_proc import Hotstart
from schism_py_pre_post.Download.download_usgs_with_api import get_usgs_obs_for_stofs3d
from schism_py_pre_post.Download.download_cbp_with_api import get_cbp_obs_for_stofs3d
from schism_py_pre_post.Shared_modules.gen_subregion_ic2 import gen_subregion_ic_stofs3d
from schism_py_pre_post.Grid.Grid_geometry import find_points_in_polyshp
import schism_py_pre_post  # this is needed for the __file__ attribute

# pip install git+https://github.com/wzhengui/pylibs.git
from pylib import schism_grid
from replace_eta2_aviso import interp_to_points_2d, transform_ll_to_cpp


def gen_elev_ic(hgrid=None, h0=0.1, city_shape_fnames=None, aviso_file=None):
    '''
    set initial elevation: 0 in the ocean, h0 below ground on higher grounds and in cities
    city_shape_fname: the shapefile must be in the same projection as the hgrid
    in this case lon/lat (because aviso data is in lon/lat)
    '''
    if city_shape_fnames is None:
        in_city = np.ones(hgrid.dp.shape, dtype=bool)
    else:
        in_city = find_points_in_polyshp(pt_xy=np.c_[hgrid.x, hgrid.y], shapefile_names=city_shape_fnames)

    high_ground = hgrid.dp < 0

    elev_ic = np.zeros(hgrid.dp.shape, dtype=float)

    # set coastal points
    land = np.logical_or(high_ground, in_city)
    elev_ic[land] = - hgrid.dp[land] - h0

    # set ocean points
    import xarray as xr
    ocean = np.logical_not(land)
    x1, y1 = transform_ll_to_cpp(hgrid.x, hgrid.y)
    # read aviso data
    ds = xr.open_dataset(aviso_file)
    lon2 = ds.longitude.values
    lat2 = ds.latitude.values
    x2, y2 = transform_ll_to_cpp(lon2, lat2)
    ssh = np.squeeze(ds.adt.values[0, :, :])
    ds.close()
    # interpolate aviso data onto model grid
    ssh_int = interp_to_points_2d(y2, x2, np.c_[y1, x1], ssh)

    elev_ic[ocean] = ssh_int[ocean] - 0.42  # nEw, uniform shift for STOFS-3D-Atl v7

    return elev_ic


def tweak_stofs3d_hotstart(
    wdir='./', hotstart_date_str='2021-05-01', city_shapefile_names=[], aviso_file='aviso.nc',
    hycom_TS_file='TS_1.nc', hycom_hot_file='hotstart.nc.hycom',
):
    '''
    The "wdir" should contain the following files:
        hgrid.gr3 (should be the same as hgrid.ll for STOFS3D), vgrid.in
        TS_1.nc: Temperature and Salinity from HYCOM, the first time in the nc file should match model start time.
                 This file should already be prepared in the previous step of hotstart.nc generation
        hotstart.nc.hycom: hotstart.nc based on HYCOM only, renamed with the '.hycom' suffix.
        aviso.nc: sea surface height above geoid, in lon/lat

    You may select pre-defined city polygons to force the water level to be below ground;
    only specify the file basename, the script will copy these files to your wdir.

    Use scripts in ../AVISO/ to download aviso.nc
    '''
    griddir = wdir
    output_obs_dir = f'{wdir}/Obs/'  # observations will be downloaded to this directory
    my_hot_file = f'{wdir}/hotstart.nc'  # this is the final hotstart.nc file
    hycom_TS_file = f'{wdir}/TS_1.nc'  # this file should already be prepared in the previous step of hotstart.nc generation from HYCOM
    hycom_hot_file = f'{wdir}/hotstart.nc.hycom'  # this file should already be prepared in the previous step of hotstart.nc generation from HYCOM, renamed with the '.hycom' suffix.

    # copy datafiles
    mydir = os.path.dirname(schism_py_pre_post.__file__)
    for shp in city_shapefile_names:
        shp_basename = Path(shp).stem
        if not os.path.exists(f'{wdir}/{shp_basename}.shp'):
            os.system(f'cp {mydir}/Datafiles/{shp_basename}.* {wdir}')

    print(f'downloading obs from USGS to {output_obs_dir}')
    get_usgs_obs_for_stofs3d(outdir=output_obs_dir, start_date_str=hotstart_date_str)

    print('downloading coastal obs from CBP; only salinity, no temperature (cbp samples are sparse in time)')
    get_cbp_obs_for_stofs3d(outdir=output_obs_dir, sample_time=hotstart_date_str, varname=['sal'])

    print('interpolating obs onto model grid')
    gen_subregion_ic_stofs3d(wdir=wdir, obsdir=output_obs_dir, hycom_TS_file=hycom_TS_file, date_str=hotstart_date_str)

    print(f'copying {hycom_TS_file} to {my_hot_file}, which will be modified')
    if os.path.exists(my_hot_file):
        os.system(f"rm {my_hot_file}")
    os.system(f"cp {hycom_hot_file} {my_hot_file}")

    my_hot = Hotstart(
        grid_info=griddir,
        hot_file=my_hot_file
    )

    # increase T by 1 oC
    # my_hot.tr_nd.val[:, :, 0] += 1.0
    # my_hot.tr_nd0.val[:, :, 0] += 1.0
    # my_hot.tr_el.val[:, :, 0] += 1.0

    print('tweaking coastal values based on obs')
    for i, var in enumerate(['tem', 'sal']):
        hg = schism_grid(f'{wdir}/ecgc_coastal_{var}.gr3')  # this file is from get*_obs_for_stofs3d
        idx = hg.dp > -9998
        for k in range(my_hot.grid.vgrid.nvrt):
            my_hot.tr_nd.val[idx, k, i] = hg.dp[idx]

    print('setting salinity to 0 on higher grounds')
    rat = np.maximum(np.minimum(1.0, (my_hot.grid.hgrid.dp + 3.0) / 3.0), 0.0)  # linearly varying from 0 to 3 m
    my_hot.tr_nd.val[:, :, 1] *= np.transpose(np.tile(rat, (my_hot.grid.vgrid.nvrt, 1)))
    my_hot.trnd_propogate()  # propogate trnd values to trnd0 and tr_el

    print('setting initial elevation: aviso values in the ocean, just below ground on higher grounds and in cities')
    my_hot.eta2.val[:] = gen_elev_ic(
        hgrid=my_hot.grid.hgrid, h0=0.1,
        city_shape_fnames=[f'{wdir}/{x}' for x in city_shapefile_names],
        aviso_file=f'{wdir}/{aviso_file}'  # no shift based on R11, R24, R11a
    )

    print('writing the final hotstart.nc')
    my_hot.writer(fname=my_hot_file)

if __name__ == "__main__":
    # Modify hycom-based hotstart.nc.hycom with coastal observation values of T and S, and a better initial elevation
    tweak_stofs3d_hotstart(
        wdir='/sciclone/schism10/feiye/STOFS3D-v7/Inputs/v7_2024_reforecast/',
        hotstart_date_str='2023-12-01',
        city_shapefile_names = ["LA_urban_polys_lonlat.shp"],  # polygon shapefile specifying cities
        aviso_file='aviso.nc'
    )

