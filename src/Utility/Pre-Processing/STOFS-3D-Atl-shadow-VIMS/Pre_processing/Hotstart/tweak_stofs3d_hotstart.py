import numpy as np
import os
import STOFS3D_scripts
from schism_py_pre_post.Shared_modules.hotstart_proc import Hotstart
from schism_py_pre_post.Download.download_usgs_with_api import get_usgs_obs_for_stofs3d
from schism_py_pre_post.Download.download_cbp_with_api import get_cbp_obs_for_stofs3d
from schism_py_pre_post.Shared_modules.gen_subregion_ic2 import gen_subregion_ic_stofs3d
from schism_py_pre_post.Grid.Grid_geometry import find_points_in_polyshp
from pylib import schism_grid
from pathlib import Path


def gen_elev_ic(hgrid=None, h0=0.1, city_shape_fnames=None):
    '''
    set initial elevation: 0 in the ocean, h0 below ground on higher grounds and in cities
    city_shape_fname: the shapefile must be in the same projection as the hgrid
    '''
    if city_shape_fnames is None:
        in_city = np.ones(hgrid.dp.shape, dtype=bool)
    else:
        in_city = find_points_in_polyshp(pt_xy=np.c_[hgrid.x, hgrid.y], shapefile_names=city_shape_fnames)

    above_NAVD88_0 = hgrid.dp < 0

    elev_ic = np.zeros(hgrid.dp.shape, dtype=float)

    ind = np.logical_or(above_NAVD88_0, in_city)
    elev_ic[ind] = - hgrid.dp[ind] - h0

    return elev_ic
    
def tweak_stofs3d_hotstart(wdir ='./', hotstart_date_str='2021-05-01', city_shapefile_names=[]):
    '''
    The "wdir" should contain the following files:
        hgrid.gr3 (should be the same as hgrid.ll for STOFS3D), vgrid.in
        TS_1.nc: Temperature and Salinity from HYCOM, the first time in the nc file should match model start time.
                 This file should already be prepared in the previous step of hotstart.nc generation 
        hotstart.nc.hycom: hotstart.nc based on HYCOM only, renamed with the '.hycom' suffix.

    You may select pre-defined city polygons to force the water level to be below ground;
    only specify the file basename, the script will copy these files to your wdir.
    '''

    # input section
    griddir = wdir
    output_obs_dir = f'{wdir}/Obs/'
    hycom_TS_file = f'{wdir}/TS_1.nc'
    hycom_hot_file = f'{wdir}/hotstart.nc.hycom'
    my_hot_file = f'{wdir}/hotstart.nc'
    # end input section

    # copy datafiles
    mydir = os.path.dirname(STOFS3D_scripts.__file__)
    for shp in city_shapefile_names:
        shp_basename = Path(shp).stem
        os.system(f'cp {mydir}/Datafiles/Tweak_hotstart/{shp_basename}.* {wdir}')

    # download coastal obs from usgs
    get_usgs_obs_for_stofs3d(outdir=output_obs_dir, start_date_str=hotstart_date_str)

    # download coastal obs from CBP
    get_cbp_obs_for_stofs3d(outdir=output_obs_dir, sample_time=hotstart_date_str)

    # interpolate obs onto model grid
    gen_subregion_ic_stofs3d(wdir=wdir, obsdir=output_obs_dir, hycom_TS_file=hycom_TS_file, date_str=hotstart_date_str)

    # make a copy of the hycom-based hotstart.nc
    if os.path.exists(my_hot_file):
        os.system(f"rm {my_hot_file}")
    os.system(f"cp {hycom_hot_file} {my_hot_file}")

    # tweak coastal values based on obs
    my_hot = Hotstart(
        grid_info=griddir,
        hot_file=my_hot_file
    )

    for i, var in enumerate(['tem', 'sal']):
        hg = schism_grid(f'{wdir}/ecgc_coastal_{var}.gr3')  # this file is from get*_obs_for_stofs3d
        idx = hg.dp > -9998
        for k in range(my_hot.grid.vgrid.nvrt):
            my_hot.tr_nd.val[idx, k, i] = hg.dp[idx]

    # set salinity to 0 on higher grounds
    rat = np.maximum(np.minimum(1.0, (my_hot.grid.hgrid.dp + 3.0) / 3.0), 0.0)  # linearly varying from 0 to 3 m
    my_hot.tr_nd.val[:, :, 1] *= np.transpose(np.tile(rat, (my_hot.grid.vgrid.nvrt, 1)))
    my_hot.trnd_propogate()  # propogate trnd values to trnd0 and tr_el

    # set initial elevation: 0 in the ocean, just below ground on higher grounds and in cities
    my_hot.eta2.val[:] = gen_elev_ic(hgrid=my_hot.grid.hgrid, h0=0.1, city_shape_fnames=[f'{wdir}/{x}' for x in city_shapefile_names])

    # write
    my_hot.writer(fname=my_hot_file)

if __name__ == "__main__":
    '''
    Modify hycom-based hotstart.nc.hycom with coastal observation values
    See instructions in "def tweak_stofs3d_hotstart()"
    '''
    tweak_stofs3d_hotstart(
        wdir='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/I24/Hot_test/',
        hotstart_date_str='2021-05-01',
        city_shapefile_names = ["city_polys_from_v10_lonlat.shp"]
    )
    