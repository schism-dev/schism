import logging
import os
import pathlib

import appdirs
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate.fitpack2 import RectBivariateSpline

from pyschism.forcing.bctides.base import TidalDataProvider


logger = logging.getLogger(__name__)

TPXO_ELEVATION = 'h_tpxo9.v1.nc'
TPXO_VELOCITY = 'u_tpxo9.v1.nc'


def raise_missing_file(fpath, fname):
    raise FileNotFoundError('\n'.join([
        f'No TPXO file found at "{fpath}".',
        'New users will need to register and request a copy of '
        f'the TPXO9 NetCDF file (specifically `{fname}`) '
        'from the authors at https://www.tpxo.net.',
        'Once you obtain `h_tpxo9.v1.nc`, you can follow one of the '
        'following options: ',
        f'1) copy or symlink the file to "{fpath}"',
        f'2) set the environment variable `{fname}` to point'
        ' to the file',
    ]))


class TPXO(TidalDataProvider):

    def __init__(self, h_file=None, u_file=None):
        self._h_file = h_file
        self._u_file = u_file

    def get_elevation(self, constituent, vertices):
        logger.info('Querying TPXO for elevation constituent '
                    f'{constituent}.')
        amp = self._get_interpolation(
            'elevation', 'ha', constituent, vertices)
        phase = self._get_interpolation(
            'elevation', 'hp', constituent, vertices)
        return amp, phase

    def get_velocity(self, constituent, vertices):
        logger.info('Querying TPXO for velocity constituent '
                    f'{constituent}.')
        uamp = self._get_interpolation(
            'velocity', 'ua', constituent, vertices) / 100.
        uphase = self._get_interpolation(
            'velocity', 'up', constituent, vertices)
        vamp = self._get_interpolation(
            'velocity', 'va', constituent, vertices) / 100.
        vphase = self._get_interpolation(
            'velocity', 'vp', constituent, vertices)
        return uamp, uphase, vamp, vphase

    @property
    def constituents(self):
        return ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1',
                'Q1', 'Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1']
        if not hasattr(self, '_constituents'):
            self._constituents = [
                c.capitalize() for c in self.h['con'][:].astype(
                    '|S1').tobytes().decode('utf-8').split()]
        return self._constituents

    @property
    def lon_z(self) -> np.ndarray:
        return self.h['lon_z'][:, 0].data

    @property
    def lat_z(self) -> np.ndarray:
        return self.h['lat_z'][0, :].data

    @property
    def lon_u(self) -> np.ndarray:
        return self.uv['lon_u'][:, 0].data

    @property
    def lat_u(self) -> np.ndarray:
        return self.uv['lat_u'][0, :].data

    @property
    def lon_v(self) -> np.ndarray:
        return self.uv['lon_v'][:, 0].data

    @property
    def lat_v(self) -> np.ndarray:
        return self.uv['lat_v'][0, :].data

    @property
    def h(self):
        if not hasattr(self, '_h'):
            if self._h_file is None:
                self._h_file = os.getenv('TPXO_ELEVATION')
                if self._h_file is None:
                    self._h_file = pathlib.Path(
                        appdirs.user_data_dir('tpxo')) / TPXO_ELEVATION
            if not self._h_file.exists():
                raise_missing_file(self._h_file, TPXO_ELEVATION)
            self._h = Dataset(self._h_file)
        return self._h

    @property
    def uv(self):
        if not hasattr(self, '_uv'):
            if self._u_file is None:
                self._u_file = os.getenv('TPXO_VELOCITY')
                if self._u_file is None:
                    self._u_file = pathlib.Path(
                        appdirs.user_data_dir('tpxo')) / TPXO_VELOCITY
            if not self._u_file.exists():
                raise_missing_file(self._u_file, TPXO_VELOCITY)
            self._uv = Dataset(self._u_file)
        return self._uv

    def _get_interpolation(self, phys_var, ncvar, constituent, vertices):
        lower_c = [c.lower() for c in self.constituents]
        if phys_var == 'elevation':
            ncarray = self.h
            self.x = self.lon_z
            self.y = self.lat_z
        elif ncvar == 'ua' or ncvar == 'up':
            ncarray = self.uv
            self.x = self.lon_u
            self.y = self.lat_u
        elif ncvar == 'va' or ncvar == 'vp':
            ncarray = self.uv
            self.x = self.lon_v
            self.y = self.lat_v

        zi = ncarray[ncvar][
            lower_c.index(constituent.lower()), :, :]
        xo = np.asarray(
            [x + 360. if x < 0. else x for x in vertices[:, 0]]).flatten()
        yo = vertices[:, 1].flatten()
        xi, yi = np.meshgrid(self.x, self.y, indexing='ij')
        xi = xi.flatten()
        yi = yi.flatten()
        zi = zi.flatten()
        dx = np.mean(np.diff(self.x))
        dy = np.mean(np.diff(self.y))
        # buffer the bbox by 2 difference units
        mask1 = np.logical_and(
            np.logical_and(
                xi >= np.min(xo) - 2 * dx,
                xi <= np.max(xo) + 2 * dx
            ),
            np.logical_and(
                yi >= np.min(yo) - 2 * dy,
                yi <= np.max(yo) + 2 * dy
            )
        )
        # remove junk values from input array
        mask2 = np.ma.masked_where(zi != 0., zi)
        iidx = np.where(np.logical_and(mask1, mask2))
        values = griddata(
            (xi[iidx], yi[iidx]),
            zi[iidx],
            (xo, yo),
            method='linear',
            fill_value=np.nan,
        )
        nan_idxs = np.where(np.isnan(values))
        values[nan_idxs] = griddata(
            (xi[iidx], yi[iidx]),
            zi[iidx],
            (xo[nan_idxs], yo[nan_idxs]),
            method='nearest',
        )
        return values
