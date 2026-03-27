import logging
import pathlib

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

from pyschism.forcing.bctides.base import TidalDataProvider

# https://icdc.cen.uni-hamburg.de/en/hamtide.html
base_url = 'https://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/hamtide/'

logger = logging.getLogger(__name__)


class HAMTIDE(TidalDataProvider):
    ''' Wrapper for querying HAMTIDE model harmonic constituents.

    The source data can be either files from the static storage of on-the-fly
    querying from OpenDAP server (default).

    Reference:
        Taguchi, E., D. Stammer and W. Zahel (2010), Estimation of deep ocean
        tidal energy dissipation based on the high-resolution data-assimilative
        HAMTIDE model (to be submitted to J. Geophys. Res.).
    '''

    def __init__(self, resource=None):
        self.resource = resource

    def get_elevation(self, constituent, vertices):
        logger.info('Querying HAMTIDE for elevation constituent '
                    f'{constituent}.')
        amp = 0.01 * self._get_interpolation(
            'elevation', 'AMPL', constituent, vertices)
        phase = self._get_interpolation(
            'elevation', 'PHAS', constituent, vertices)
        return amp, phase

    def get_velocity(self, constituent, vertices):
        logger.info('Querying HAMTIDE for velocity constituent '
                    f'{constituent}.')
        uamp = 0.01 * self._get_interpolation(
            'velocity', 'UAMP', constituent, vertices)
        uphase = self._get_interpolation(
            'velocity', 'UPHA', constituent, vertices)
        vamp = 0.01 * self._get_interpolation(
            'velocity', 'VAMP', constituent, vertices)
        vphase = self._get_interpolation(
            'velocity', 'VPHA', constituent, vertices)
        return uamp, uphase, vamp, vphase

    @property
    def constituents(self):
        return ['S2', 'Q1', 'P1', 'O1', 'N2', 'M2', 'K2', 'K1']

    @property
    def x(self):
        if not hasattr(self, '_x'):
            self._x = Dataset(base_url + 'k2.hamtide11a.nc')['LON'][:].data
        return self._x

    @property
    def y(self):
        if not hasattr(self, '_y'):
            self._y = Dataset(base_url + 'k2.hamtide11a.nc')['LAT'][:].data
        return self._y

    def _get_resource(self, variable, constituent) -> Dataset:
        resource = self._resource[variable][constituent]
        if resource is not None:
            return Dataset(resource)
        if variable == 'elevation':
            fname = f'{constituent.lower()}.hamtide11a.nc'
        if variable == 'velocity':
            fname = f'HAMcurrent11a_{constituent.lower()}.nc'
        return Dataset(base_url + fname)

    def _get_interpolation(self, phys_var, ncvar, constituent, vertices):
        xq = np.asarray(
            [x + 360. if x < 0. else x for x in vertices[:, 0]]).flatten()
        yq = vertices[:, 1].flatten()
        dx = (self.x[-1] - self.x[0]) / len(self.x)
        xidx = np.logical_and(
            self.x >= np.min(xq) - 2.*dx,
            self.x <= np.max(xq) + 2.*dx
        )
        dy = (self.y[-1] - self.y[0]) / len(self.y)
        yidx = np.logical_and(
            self.y >= np.min(yq) - 2.*dy,
            self.y <= np.max(yq) + 2.*dy
        )
        xi, yi = np.meshgrid(self.x[xidx], self.y[yidx])
        xi = xi.flatten()
        yi = yi.flatten()
        zi = self._get_resource(
                phys_var, constituent)[ncvar][yidx, xidx].flatten()
        values = griddata(
                (xi[~zi.mask], yi[~zi.mask]),
                zi[~zi.mask],
                (xq, yq),
                method='linear',
                fill_value=np.nan
            )
        nan_idxs = np.where(np.isnan(values))
        values[nan_idxs] = griddata(
                (xi[~zi.mask], yi[~zi.mask]),
                zi[~zi.mask],
                (xq[nan_idxs], yq[nan_idxs]),
                method='nearest',
        )
        return values

    @property
    def resource(self):
        return self._resource

    @resource.setter
    def resource(self, resource):
        if resource is None:
            resource = {'elevation': {}, 'velocity': {}}
            for constituent in self.constituents:
                for key in resource.keys():
                    resource[key].update({constituent: None})
        else:
            raise NotImplementedError('Check that static files exist.')
            for file in pathlib.Path(resource).glob('*.nc'):
                print(file)
        self._resource = resource
