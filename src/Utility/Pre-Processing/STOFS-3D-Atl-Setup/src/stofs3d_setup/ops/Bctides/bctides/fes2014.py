import logging
import os
import pathlib
import glob

import appdirs
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate.fitpack2 import RectBivariateSpline

from pyschism.forcing.bctides.base import TidalDataProvider


logger = logging.getLogger(__name__)

FES2014_TIDES_EXTRA = 'ocean_tide_extrapolated'
FES2014_EASTWARD_VEL = 'eastward_velocity'
FES2014_NORTHWARD_VEL = 'northward_velocity'


def raise_missing_file(fpath, fname):
    raise FileNotFoundError('\n'.join([
        f'No FES2014 file found at "{fpath}".',
        'New users will need to register and download a copy of '
        f'the FES2014 NetCDF file (specifically `{fname}`) '
        'from https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html',
        'Once you obtain netCDF files, you can follow one of the '
        'following options: ',
        f'1) copy or symlink the file to "{fpath}"',
        f'2) set the environment variable `{fname}` to point'
        ' to the file',
    ]))


class FES2014(TidalDataProvider):

    def __init__(self, resource=None):
        self.resource = resource

    def get_elevation(self, constituent, vertices):
        logger.info('Querying FES2014 for elevation constituent '
                    f'{constituent}.')
        amp = self._get_interpolation(
            'elevation', 'amplitude', constituent, vertices)
        phase = self._get_interpolation(
            'elevation', 'phase', constituent, vertices)
        return amp, phase

    def get_velocity(self, constituent, vertices):
        logger.info('Querying FES2014 for velocity constituent '
                    f'{constituent}.')
        uamp = self._get_interpolation(
            'eastward_vel', 'Ua', constituent, vertices)
        uphase = self._get_interpolation(
            'eastward_vel', 'Ug', constituent, vertices)
        vamp = self._get_interpolation(
            'northward_vel', 'Va', constituent, vertices)
        vphase = self._get_interpolation(
            'northward_vel', 'Vg', constituent, vertices)
        return uamp, uphase, vamp, vphase

    @property
    def constituents(self):
        return ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1',
                'Q1', 'Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1']
        if not hasattr(self, '_constituents'):
            files = glob.glob(appdirs.user_data_dir('fes2014') + '/ocean_tide_extrapolated/*.nc')
            self._constituents = [f.split('/')[-1].split('.')[0] for f in files]
        return self._constituents

    def _get_resource(self, variable, constituent) -> Dataset:
        resource = self._resource[variable][constituent]
        if resource is not None:
            return Dataset(resource)
        datapath = pathlib.Path(appdirs.user_data_dir('fes2014'))
        if variable == 'elevation':
            fname = datapath / f'{FES2014_TIDES_EXTRA}/{constituent.lower()}.nc'
        if variable == 'eastward_vel':
            fname = datapath / f'{FES2014_EASTWARD_VEL}/{constituent.lower()}.nc'
        if variable == 'northward_vel':
            fname = datapath / f'{FES2014_NORTHWARD_VEL}/{constituent.lower()}.nc'
        return Dataset(fname)

    def _get_interpolation(self, phys_var, ncvar, constituent, vertices):
        ds = self._get_resource(phys_var, constituent) 
        lon = ds['lon'][:]
        lat = ds['lat'][:]
        dxs = np.unique(np.diff(lon))
        dys = np.unique(np.diff(lat))
        if len(dxs) !=1 or len(dys) != 1:
           raise ValueError(f'{phys_var}: lon, lat of {constituent}.nc not uniform! ')  
        dx = dxs[0]
        dy = dys[0]

        #get interp index
        xi = np.asarray(
            [x + 360. if x < 0. else x for x in vertices[:, 0]]).flatten()
        yi = vertices[:, 1].flatten()
        idx = np.floor((xi - lon[0]) / dx).astype('int')
        mask = np.nonzero((lon[idx] - xi) > 0)[0]
        idx[mask] = idx[mask] - 1
        idy = np.floor((yi - lat[0]) / dy).astype('int')
        mask = np.nonzero((lat[idy] - yi) > 0)[0]
        idy[mask] = idy[mask] - 1
        xrat = (xi - lon[idx]) / (lon[idx+1] - lon[idx])
        yrat = (yi - lat[idy]) / (lat[idy+1] - lat[idy])
        if np.sum((xrat > 1) | (xrat < 0) | (yrat > 1) | (yrat < 0)) != 0:
            raise ValueError(f'xrat or yrat > 1 or < 0')
       
        zi = ds[ncvar][:,:]
        # vm = 100 junk for amplitude, vm = 370 junk for phase
        if ncvar == 'amplitude' or ncvar == 'Ua' or ncvar == 'Va':
            vm = 100
            zi = zi / 100
        if ncvar == 'phase' or ncvar == 'Ug' or ncvar == 'Vg':
            vm = 370
            #Put phases into [0, 360)
            mask = zi<0
            zi[mask] = zi[mask] + 360
        v0 = np.c_[zi[idy, idx], zi[idy, idx+1], zi[idy+1, idx], zi[idy+1, idx+1]].T
        vmax = v0.max(axis=0)
        vmin = v0.min(axis=0)
        #handle phase jump
        if ncvar == 'phase' or ncvar == 'Ug' or ncvar == 'Vg':
            for kk in np.nonzero(abs(vmax - vmin) > 180)[0]:
                mask = abs(v0[:, kk] - vmax[kk]) > 180
                v0[mask, kk] = v0[mask, kk] + 360
        v1 = v0[0] * (1 - xrat) + v0[1] * xrat
        v2 = v0[2] * (1 - xrat) + v0[3] * xrat
        values = v1 * (1 - yrat) + v2 * yrat
        mask = np.nonzero((vmax > vm) * (vmin <= vm) * (vmin >= 0))[0]
        values[mask] = vmin[mask]
        if sum((vmax > vm) * ((vmin > vm) | (vmin < 0))) != 0:
            raise ValueError('All junk values for {phys_var} {constituent}')
        
        ds.close()
        return values

    @property
    def resource(self):
        return self._resource

    @resource.setter
    def resource(self, resource):
        if resource is None:
            resource = {'elevation': {}, 'eastward_vel': {}, 'northward_vel': {}}
            for constituent in self.constituents:
                for key in resource.keys():
                    resource[key].update({constituent: None})
        else:
            raise NotImplementedError('Check that static files exist.')
            for file in pathlib.Path(resource).glob('*.nc'):
                print(file)
        self._resource = resource
