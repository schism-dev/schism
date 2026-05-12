from collections import OrderedDict
from datetime import datetime, timedelta, timezone
from enum import Enum
import logging
from typing import Callable, Union

import numpy as np

from pyschism.forcing.bctides.tpxo import TPXO
from pyschism.forcing.bctides.fes2014 import FES2014
from pyschism.forcing.bctides.hamtide import HAMTIDE

_logger = logging.getLogger(__name__)

MAJOR_CONSTITUENTS = ('Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2')
MINOR_CONSTITUENTS = ('Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1')
ALL_CONSTITUENTS = MAJOR_CONSTITUENTS + MINOR_CONSTITUENTS


class TidalDatabase(Enum):
    TPXO = TPXO
    FES2014 = FES2014 
    HAMTIDE = HAMTIDE

    def _missing_(self, name):
        raise ValueError(f"{name} is not a valid tidal database.")


class ActiveConstituents:

    def __get__(self, obj, val):
        active_constituents = obj.__dict__.get('active_constituents')
        if active_constituents is None:
            active_constituents = OrderedDict()
            obj.__dict__['active_constituents'] = active_constituents
        return active_constituents

    def __set__(self, obj, val):
        return obj.__dict__['active_constituents']


class Tides:

    _active_constituents = ActiveConstituents()
    major_constituents = MAJOR_CONSTITUENTS
    minor_constituents = MINOR_CONSTITUENTS
    constituents = ALL_CONSTITUENTS
    _nodal_factors = {}
    _earth_equilibrium_arguments = {}
    _tidal_species_types = {}
    _amplitude_constants = {}
    _amplitudes = {
            'elevation': {},
            'velocity': {},
        }

    _phases = {
        'elevation': {},
        'velocity': {},
    }

    def __init__(
            self,
            tidal_database: Union[str, TidalDatabase] = TidalDatabase.TPXO,
            constituents='all',
    ):
        """Main class for requesting tidal boundary forcing for a SCHISM run.

        The user should add tidal requests to this class via it's methods.

        Args:
            elevation (optional): Apply elevation forcing to boundary, defaults
                to True.
            velocity (optional): Apply velocity forcing to boundary (currently
                disabled), defaults to False.
            database (optional): Tidal database to use in order to obtain
                boundary initial conditions, defaults to TidalDatabase.TPXO
        """
        self.tidal_database = tidal_database

        if isinstance(constituents, str):
            constituents = [constituents]

        constituents = list(constituents) if constituents is not None else []

        if 'all' in constituents and len(constituents) > 1:
            raise ValueError('When using all, must only pass one')

        elif 'major' in constituents and len(constituents) > 1:
            raise ValueError('When using major, must only pass one')

        if 'all' in constituents:
            self.use_all()

        elif 'major' in constituents:
            self.use_major()

        else:
            for constituent in constituents:
                self.use_constituent(constituent)

    def __iter__(self):
        for constituent in self.active_constituents:
            yield (constituent, self.get_tidal_constituent(constituent))

    def __len__(self):
        return len(self.active_constituents)

    def __call__(self, start_date: datetime, rnday: Union[int, datetime],
                 constituent: str):
        """Returns the tidal characteristic values for a given time period.
        """
        _logger.debug(
            f'Computing tidal factors for start_date={start_date} and '
            f'constituent {constituent}')
        return (self.get_tidal_species_type(constituent),
                self.get_tidal_potential_amplitude(constituent),  # TPK
                self.get_orbital_frequency(constituent),  # FF*
                self.get_nodal_factor(start_date, rnday, constituent),  # Amig*
                self.get_greenwich_factor(start_date, rnday, constituent))  # FACE* # noqa:E501

    def get_elevation(self, constituent, vertices):

        if constituent.lower() == 'z0':
            return np.full((vertices.shape[0],), float(self._Z0)), \
                    np.full((vertices.shape[0],), 0.)

        if constituent in self._amplitudes['elevation']:
            amps = self._amplitudes['elevation'][constituent]
            phases = self._phases['elevation'][constituent]
            if isinstance(amps, (float, int)):
                amps = np.full((vertices.shape[0],), float(amps))
            if isinstance(phases, (float, int)):
                phases = np.full((vertices.shape[0],), float(phases))
            return amps, phases
        
        # if vertices is None:
            # raise ValueError(
            #     'Argument vertices must not be None for tidal constituent '
            #     f'{constituent}.')
        return self.tidal_database.get_elevation(constituent, vertices)

    def get_velocity(self, constituent, vertices):
        if constituent.lower() == 'z0':
            return tuple(4*[np.full((vertices.shape[0],), 0.)])
        
        if constituent in self._amplitudes['velocity']:
            amps = self._amplitudes['velocity'][constituent]
            phases = self._phases['velocity'][constituent]
            if isinstance(amps, (float, int)):
                amps = np.full((vertices.shape[0],), float(amps))
            if isinstance(phases, (float, int)):
                phases = np.full((vertices.shape[0],), float(phases))
            return amps, phases

        # if vertices is None:
            # raise ValueError(
            #     'Argument vertices must not be None for tidal constituent '
            #     f'{constituent}.')
        return self.tidal_database.get_velocity(constituent, vertices)

    def use_all(self, potential=True, forcing=True):
        for constituent in self.tidal_database.constituents:
            if constituent not in self.tidal_potential_amplitudes:
                potential = False
            self._active_constituents[constituent] = {
                'potential': potential,
                'forcing': forcing
            }

    def use_major(self, potential=True, forcing=True):
        for constituent in self.tidal_database.constituents:
            if constituent in self.major_constituents:
                self.use_constituent(constituent, potential, forcing)

    def use_constituent(self, constituent, potential=True, forcing=True):
        if constituent not in self.constituents:
            raise ValueError(f"Constituent must be one of {self.constituents}")
        if constituent not in self.tidal_potential_amplitudes:
            potential = False
        self._active_constituents[constituent] = {
            "potential": potential,
            "forcing": forcing,
        }

    def drop_constituent(self, constituent):
        if constituent not in self.active_constituents:
            raise ValueError("Argument constituent must be one of "
                             f"{self.active_constituents}")
        self._active_constituents.pop(constituent)
        
    def add_constituent(
            self,
            name: str,
            angular_frequency,
            nodal_factor = 1.,
            earth_equilibrium_argument: float = 0.,
            amplitude_constant: float = None,
            tidal_species_type: int = None,
            elevation_amplitude: dict = None,
            elevation_phase: dict = 0.,
            velocity_amplitude: dict = None,
            velocity_phase = 0.,
    ):
        # self._amplitudes['elevation'].update({name: amplitude})
        self.orbital_frequencies.update({name: angular_frequency})
        # self.
        self._nodal_factors.update({name: nodal_factor})
        self._earth_equilibrium_arguments.update({name: earth_equilibrium_argument})
        self._tidal_species_types.update({name: tidal_species_type})
        self._amplitude_constants.update({name: amplitude_constant})  
        self._amplitudes['elevation'].update({name: elevation_amplitude})
        self._amplitudes['velocity'].update({name: velocity_amplitude})
        self._phases['elevation'].update({name: elevation_phase})
        self._phases['velocity'].update({name: velocity_phase})
        self._active_constituents[name] = {
                "potential": False,
                "forcing": True,
            }

    def add_Z0(self, value):
        self._active_constituents['Z0'] = {
            "potential": False,
            "forcing": True,
        }
        self._Z0 = float(value)

    def get_active_constituents(self):
        return list(self.active_constituents.keys())

    def get_active_forcing_constituents(self):
        return [const for const, _ in self.active_constituents.items()
                if _['forcing'] is True]

    def get_active_potential_constituents(self):
        return [const for const, _ in self.active_constituents.items()
                if _['potential'] is True]

    def get_tidal_potential_amplitude(self, constituent):
        if constituent in self.tidal_potential_amplitudes:
            return self.tidal_potential_amplitudes[constituent]

    def get_tidal_species_type(self, constituent):
        if constituent in self.tidal_species_type:
            return self.tidal_species_type[constituent]

    def get_orbital_frequency(self, constituent):
        # if constituent in self.
        return self.orbital_frequencies[constituent]

    def get_initial_conditions(self, constituent, vertices):
        return self.initial_conditions(constituent, vertices)

    def set_Z0(self, Z0):
        self._active_constituents['Z0'] = {
            'potential': False,
            'forcing': True
        }
        self._Z0 = Z0

    def _manage_dates(f: Callable):
        def decorator(self, start_date, rnday, constituent):
            val = f(self, start_date, rnday, constituent)
            if hasattr(self, 'start_date_utc'):
                del(self.start_date_utc)
            if hasattr(self, 'end_date_utc'):
                del(self.end_date_utc)
            return val
        return decorator

    @_manage_dates
    def get_nodal_factor(self, start_date: datetime,
                         rnday: Union[float, timedelta],
                         constituent: str):
        if constituent not in ALL_CONSTITUENTS:
            return self._nodal_factors[constituent]
        if start_date.tzinfo is not None and \
                start_date.tzinfo.utcoffset(start_date) is not None:
            self.start_date_utc = start_date.astimezone(timezone(timedelta(0)))
        else:
            self.start_date_utc = start_date

        if not isinstance(rnday, timedelta):
            rnday = timedelta(days=rnday)

        self.end_date_utc = self.start_date_utc + rnday
        if constituent == "M2":
            return self.EQ78
        elif constituent == "S2":
            return 1.0
        elif constituent == "N2":
            return self.EQ78
        elif constituent == "K1":
            return self.EQ227
        elif constituent == "M4":
            return (self.EQ78)**2.
        elif constituent == "O1":
            return self.EQ75
        elif constituent == "M6":
            return (self.EQ78)**3.
        elif constituent == "MK3":
            return self.EQ78*self.EQ227
        elif constituent == "S4":
            return 1.0
        elif constituent == "MN4":
            return (self.EQ78)**2.
        elif constituent == "Nu2":
            return self.EQ78
        elif constituent == "S6":
            return 1.0
        elif constituent == "MU2":
            return self.EQ78
        elif constituent == "2N2":
            return self.EQ78
        elif constituent == "OO1":
            return self.EQ77
        elif constituent == "lambda2":
            return self.EQ78
        elif constituent == "S1":
            return 1.0
        elif constituent == "M1":
            return self.EQ207
        elif constituent == "J1":
            return self.EQ76
        elif constituent == "Mm":
            return self.EQ73
        elif constituent == "Ssa":
            return 1.0
        elif constituent == "Sa":
            return 1.0
        elif constituent == "Msf":
            return self.EQ78
        elif constituent == "Mf":
            return self.EQ74
        elif constituent == "RHO":
            return self.EQ75
        elif constituent == "Q1":
            return self.EQ75
        elif constituent == "T2":
            return 1.0
        elif constituent == "R2":
            return 1.0
        elif constituent == "2Q1":
            return self.EQ75
        elif constituent == "P1":
            return 1.0
        elif constituent == "2SM2":
            return self.EQ78
        elif constituent == "M3":
            return self.EQ149
        elif constituent == "L2":
            return self.EQ215
        elif constituent == "2MK3":
            return self.EQ227*self.EQ78**2
        elif constituent == "K2":
            return self.EQ235
        elif constituent == "M8":
            return self.EQ78**4
        elif constituent == "MS4":
            return self.EQ78
        if constituent == "Z0":
            return 1.
        else:
            msg = f'Unrecognized constituent {constituent}'
            raise TypeError(msg)

    def _normalize_to_360(f: Callable):
        def decorator(self, start_date, rnday, constituent):
            return f(self, start_date, rnday, constituent) % 360.
        return decorator

    @_manage_dates
    @_normalize_to_360
    def get_greenwich_factor(self, start_date: datetime,
                             rnday: Union[float, timedelta],
                             constituent: str):
        if constituent in self._earth_equilibrium_arguments:
            return self._earth_equilibrium_arguments[constituent]
        if start_date.tzinfo is not None and \
                start_date.tzinfo.utcoffset(start_date) is not None:
            self.start_date_utc = start_date.astimezone(timezone(timedelta(0)))
        else:
            self.start_date_utc = start_date

        if not isinstance(rnday, timedelta):
            rnday = timedelta(days=rnday)

        self.end_date_utc = self.start_date_utc + rnday
        if constituent == "M2":
            return 2.*(self.DT-self.DS+self.DH)+2.*(self.DXI-self.DNU)
        elif constituent == "S2":
            return 2.*self.DT
        elif constituent == "N2":
            return 2.*(self.DT+self.DH)-3.*self.DS+self.DP+2. \
                * (self.DXI-self.DNU)
        elif constituent == "K1":
            return self.DT+self.DH-90.-self.DNUP
        elif constituent == "M4":
            return 4.*(self.DT-self.DS+self.DH)+4.*(self.DXI-self.DNU)
        elif constituent == "O1":
            return self.DT-2.*self.DS+self.DH+90.+2.*self.DXI-self.DNU
        elif constituent == "M6":
            return 6.*(self.DT-self.DS+self.DH)+6.*(self.DXI-self.DNU)
        elif constituent == "MK3":
            return 3.*(self.DT+self.DH)-2.*self.DS-90.+2.*(self.DXI-self.DNU) \
                - self.DNUP
        elif constituent == "S4":
            return 4.*self.DT
        elif constituent == "MN4":
            return 4.*(self.DT+self.DH)-5.*self.DS+self.DP+4.\
                * (self.DXI-self.DNU)
        elif constituent == "Nu2":
            return 2.*self.DT-3.*self.DS+4.*self.DH-self.DP+2. \
                * (self.DXI-self.DNU)
        elif constituent == "S6":
            return 6.*self.DT
        elif constituent == "MU2":
            return 2.*(self.DT+2.*(self.DH-self.DS))+2.*(self.DXI-self.DNU)
        elif constituent == "2N2":
            return 2.*(self.DT-2.*self.DS+self.DH+self.DP)+2. \
                * (self.DXI-self.DNU)
        elif constituent == "OO1":
            return self.DT+2.*self.DS+self.DH-90.-2.*self.DXI-self.DNU
        elif constituent == "lambda2":
            return 2.*self.DT-self.DS+self.DP+180.+2.*(self.DXI-self.DNU)
        elif constituent == "S1":
            return self.DT
        elif constituent == "M1":
            return self.DT-self.DS+self.DH-90.+self.DXI-self.DNU+self.DQ
        elif constituent == "J1":
            return self.DT+self.DS+self.DH-self.DP-90.-self.DNU
        elif constituent == "Mm":
            return self.DS-self.DP
        elif constituent == "Ssa":
            return 2.*self.DH
        elif constituent == "Sa":
            return self.DH
        elif constituent == "Msf":
            return 2.*(self.DS-self.DH)
        elif constituent == "Mf":
            return 2.*self.DS-2.*self.DXI
        elif constituent == "RHO":
            return self.DT+3.*(self.DH-self.DS)-self.DP+90.+2.\
                * self.DXI-self.DNU
        elif constituent == "Q1":
            return self.DT-3.*self.DS+self.DH+self.DP+90.+2.*self.DXI-self.DNU
        elif constituent == "T2":
            return 2.*self.DT-self.DH+self.DP1
        elif constituent == "R2":
            return 2.*self.DT+self.DH-self.DP1+180.
        elif constituent == "2Q1":
            return self.DT-4.*self.DS+self.DH+2.*self.DP+90.+2.*self.DXI \
                - self.DNU
        elif constituent == "P1":
            return self.DT-self.DH+90.
        elif constituent == "2SM2":
            return 2.*(self.DT+self.DS-self.DH)+2.*(self.DNU-self.DXI)
        elif constituent == "M3":
            return 3.*(self.DT-self.DS+self.DH)+3.*(self.DXI-self.DNU)
        elif constituent == "L2":
            return 2.*(self.DT+self.DH)-self.DS-self.DP+180.+2.\
                * (self.DXI-self.DNU)-self.DR
        elif constituent == "2MK3":
            return 3.*(self.DT+self.DH)-4.*self.DS+90.+4.*(self.DXI-self.DNU) \
                + self.DNUP
        elif constituent == "K2":
            return 2.*(self.DT+self.DH)-2.*self.DNUP2
        elif constituent == "M8":
            return 8.*(self.DT-self.DS+self.DH)+8.*(self.DXI-self.DNU)
        elif constituent == "MS4":
            return 2.*(2.*self.DT-self.DS+self.DH)+2.*(self.DXI-self.DNU)
        elif constituent == "Z0":
            return 0
        else:
            msg = f'Unrecognized constituent {constituent}'
            raise TypeError(msg)

    def get_lunar_node(self):
        return (259.1560564 - 19.328185764 * self.DYR - .0529539336 * self.DDAY
                - .0022064139 * self.hour_middle)

    def get_lunar_perigee(self):
        return (334.3837214 + 40.66246584 * self.DYR + .111404016 * self.DDAY
                + .004641834 * self.hour_middle)

    def get_lunar_mean_longitude(self):
        return (277.0256206 + 129.38482032 * self.DYR + 13.176396768
                * self.DDAY + .549016532 * self.start_date_utc.hour)

    def get_solar_perigee(self):
        return (281.2208569 + .01717836 * self.DYR + .000047064 * self.DDAY
                + .000001961 * self.start_date_utc.hour)

    def get_solar_mean_longitude(self):
        return (280.1895014 - .238724988 * self.DYR + .9856473288 * self.DDAY
                + .0410686387 * self.start_date_utc.hour)

    @property
    def EQ73(self):
        """ """
        return (2./3.-np.sin(self.I)**2)/.5021

    @property
    def EQ74(self):
        """ """
        return np.sin(self.I)**2/.1578

    @property
    def EQ75(self):
        """ """
        return np.sin(self.I)*np.cos(self.I/2.)**2/.37988

    @property
    def EQ76(self):
        """ """
        return np.sin(2.*self.I)/.7214

    @property
    def EQ77(self):
        """ """
        return np.sin(self.I)*np.sin(self.I/2.)**2/.0164

    @property
    def EQ78(self):
        """ """
        return (np.cos(self.I/2)**4)/.91544

    @property
    def EQ149(self):
        """ """
        return np.cos(self.I/2.)**6/.8758

    @property
    def EQ197(self):
        """ """
        return np.sqrt(2.310+1.435*np.cos(2.*(self.P - self.XI)))

    @property
    def EQ207(self):
        """ """
        return self.EQ75*self.EQ197

    @property
    def EQ213(self):
        """ """
        return np.sqrt(1.-12.*np.tan(self.I/2.)**2*np.cos(2.*self.P)
                       + 36.*np.tan(self.I/2.)**4)

    @property
    def EQ215(self):
        """ """
        return self.EQ78*self.EQ213

    @property
    def EQ227(self):
        """ """
        return np.sqrt(.8965*np.sin(2.*self.I)**2+.6001*np.sin(2.*self.I)
                       * np.cos(self.NU)+.1006)

    @property
    def EQ235(self):
        """ """
        return .001+np.sqrt(19.0444*np.sin(self.I)**4+2.7702*np.sin(self.I)**2
                            * np.cos(2.*self.NU)+.0981)

    @property
    def active_constituents(self):
        return self._active_constituents.copy()

    @property
    def all_constituents(self):
        if isinstance(self.tidal_database, HAMTIDE):
            return self.major_constituents
        elif isinstance(self.tidal_database, TPXO) or isinstance(self.tidal_database, FES2014):
            return (*self.major_constituents, *self.minor_constituents)
        else:
            raise NotImplementedError(
                'Unhandled forcing database instance of type '
                f'{type(self.tidal_database)}')

    orbital_frequencies = {
        'M4':      0.0002810378050173,
        'M6':      0.0004215567080107,
        'MK3':     0.0002134400613513,
        'S4':      0.0002908882086657,
        'MN4':     0.0002783986019952,
        'S6':      0.0004363323129986,
        'M3':      0.0002107783537630,
        '2MK3':    0.0002081166466594,
        'M8':      0.0005620756090649,
        'MS4':     0.0002859630068415,
        'M2':      0.0001405189025086,
        'S2':      0.0001454441043329,
        'N2':      0.0001378796994865,
        'Nu2':     0.0001382329037065,
        'MU2':     0.0001355937006844,
        '2N2':     0.0001352404964644,
        'lambda2': 0.0001428049013108,
        'T2':      0.0001452450073529,
        'R2':      0.0001456432013128,
        '2SM2':    0.0001503693061571,
        'L2':      0.0001431581055307,
        'K2':      0.0001458423172006,
        'K1':      0.0000729211583579,
        'O1':      0.0000675977441508,
        'OO1':     0.0000782445730498,
        'S1':      0.0000727220521664,
        'M1':      0.0000702594512543,
        'J1':      0.0000755603613800,
        'RHO':     0.0000653117453487,
        'Q1':      0.0000649585411287,
        '2Q1':     0.0000623193381066,
        'P1':      0.0000725229459750,
        'Mm':      0.0000026392030221,
        'Ssa':     0.0000003982128677,
        'Sa':      0.0000001991061914,
        'Msf':     0.0000049252018242,
        'Mf':      0.0000053234146919,
        'Z0':      0}

    tidal_potential_amplitudes = {
        'M2': 0.242334,
        'S2': 0.112841,
        'N2': 0.046398,
        'K2': 0.030704,
        'K1': 0.141565,
        'O1': 0.100514,
        'P1': 0.046843,
        'Q1': 0.019256,
        'Z0': 0}

    tidal_species_type = {
        'M2': 2,
        'S2': 2,
        'N2': 2,
        'K2': 2,
        'K1': 1,
        'O1': 1,
        'P1': 1,
        'Q1': 1,
        'Z0': 0}

    @property
    def hour_middle(self):
        return self.start_date_utc.hour + (
            (self.end_date_utc - self.start_date_utc).total_seconds()
            / 3600 / 2)

    @property
    def I(self):  # noqa:E741,E743
        return np.arccos(.9136949-.0356926*np.cos(self.N))

    @property
    def N(self):
        return np.deg2rad(self.DN)

    @property
    def DN(self):
        return self.get_lunar_node()

    @property
    def DYR(self):
        return self.start_date_utc.year - 1900.

    @property
    def DDAY(self):
        return (self.start_date_utc.timetuple().tm_yday
                + int((self.start_date_utc.year-1901.)/4.) - 1)

    @property
    def NU(self):
        return np.arcsin(.0897056*np.sin(self.N)/np.sin(self.I))

    @property
    def DT(self):
        return (180.+self.start_date_utc.hour*(360./24))

    @property
    def DS(self):
        return self.get_lunar_mean_longitude()

    @property
    def DP(self):
        return self.get_lunar_perigee()

    @property
    def P(self):
        return np.deg2rad(self.DP)

    @property
    def DH(self):
        return self.get_solar_mean_longitude()

    @property
    def DP1(self):
        return self.get_solar_perigee()  # HR

    @property
    def DNU(self):
        return np.rad2deg(self.NU)

    @property
    def XI(self):
        return self.N-2.*np.arctan(.64412*np.tan(self.N/2)) - self.NU

    @property
    def DXI(self):
        return np.rad2deg(self.XI)

    @property
    def NUP(self):
        return np.arctan(np.sin(self.NU) / (
            np.cos(self.NU)+.334766/np.sin(2.*self.I)))

    @property
    def DNUP(self):
        return np.rad2deg(self.NUP)

    @property
    def DPC(self):
        return self.DP - self.DXI

    @property
    def PC(self):
        return np.deg2rad(self.DPC)

    @property
    def R(self):
        return (np.arctan(np.sin(2.*self.PC)
                          / ((1./6.)*(1./np.tan(.5*self.I))**2 - np.cos(2.*self.PC))))

    @property
    def DR(self):
        return np.rad2deg(self.R)

    @property
    def NUP2(self):
        return (np.arctan(np.sin(2.*self.NU)
                          / (np.cos(2.*self.NU)+.0726184 / np.sin(self.I)**2))/2.)

    @property
    def DNUP2(self):
        return np.rad2deg(self.NUP2)

    @property
    def Q(self):
        return np.arctan2((5.*np.cos(self.I)-1.)*np.sin(self.PC),
                          (7.*np.cos(self.I)+1.)*np.cos(self.PC))

    @property
    def DQ(self):
        return np.rad2deg(self.Q)

    @property
    def ntip(self):
        return len(self.get_active_potential_constituents())

    @property
    def cutoff_depth(self):
        if self.ntip == 0:
            return 0
        try:
            return self.__cutoff_depth
        except AttributeError:
            return 50.

    @cutoff_depth.setter
    def cutoff_depth(self, cutoff_depth):
        assert isinstance(cutoff_depth, (int, float))
        self.__cutoff_depth = cutoff_depth

    @property
    def nbfr(self):
        return len(self.get_active_forcing_constituents())

    @property
    def tidal_database(self):
        return self._tidal_database

    @tidal_database.setter
    def tidal_database(self, database: str):
        if isinstance(database, str):
            database = TidalDatabase[database.upper()].value()
        if isinstance(database, TidalDatabase):
            database = database.value()
        self._tidal_database = database
