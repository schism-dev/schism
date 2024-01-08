import os
import copy

import xarray as xr
#import scipy as sp
import numpy as np
from netCDF4 import Dataset
from pyschism.forcing.hycom.hycom2schism import interp_to_points_2d

from pylib import read_schism_hgrid, loadz, zdata, WriteNC
from schism_py_pre_post.Shared_modules.hotstart_proc import Hotstart

def transform_ll_to_cpp(lon, lat, lonc=-77.07, latc=24.0):
    longitude = lon/180*np.pi
    latitude = lat/180*np.pi
    radius = 6378206.4
    loncc = lonc/180*np.pi
    latcc = latc/180*np.pi
    lon_new = [radius*(longitude[i]-loncc)*np.cos(latcc) for i in np.arange(len(longitude))]
    lat_new = [radius*latitude[i] for i in np.arange(len(latitude))]
    return np.array(lon_new), np.array(lat_new)


if __name__ == '__main__':

    wdir = './'
    griddir = wdir
    hycom_hot_file = f'{wdir}/hotstart.nc.tweaked'
    my_hot_file = f'{wdir}/hotstart.nc'

    if os.path.exists('grid.npz'):
        gd = loadz('grid.npz').hgrid
    else:
        gd = read_schism_hgrid('hgrid.gr3')
        gd.save('grid.npz')

    lon1 = gd.x
    lat1 = gd.y
    x1,y1 = transform_ll_to_cpp(lon1, lat1)
    bxy = np.c_[y1, x1]

    #time period from 2017-12-01T00 to 2019-01-01T00
    ds = xr.open_dataset('aviso.nc')
    lon2 = ds.longitude.values
    lat2 = ds.latitude.values
    x2, y2=transform_ll_to_cpp(lon2, lat2)

    #sea surface hegith above geoid
    ssh = np.squeeze(ds.adt.values[0, :, :])
    ssh_int = interp_to_points_2d(y2, x2, bxy, ssh)

    ds.close()

    ##read old hotstart.nc
    #ds = Dataset('hotstart.nc')
    #eta2 = ds['eta2'][:]
    #cumsum_eta = ds['cumsum_eta'][:]

    # make a copy of the hycom-based hotstart.nc
    if os.path.exists(my_hot_file):
        os.system(f"rm {my_hot_file}")
    os.system(f"cp {hycom_hot_file} {my_hot_file}")

    my_hot = Hotstart(grid_info = griddir, hot_file = my_hot_file)
    eta2_new = copy.deepcopy(my_hot.eta2.val)
    idxs = eta2_new == 0
    eta2_new[idxs] = ssh_int[idxs]  - 0.42 #adjust elev

    my_hot.eta2.val[:] = eta2_new

    my_hot.writer(fname=my_hot_file)
