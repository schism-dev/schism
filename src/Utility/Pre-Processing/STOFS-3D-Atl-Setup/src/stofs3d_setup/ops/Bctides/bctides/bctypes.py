from abc import abstractmethod, abstractproperty

from pyschism.forcing.base import ModelForcing


class BoundaryForcing(ModelForcing):

    @abstractmethod
    def get_boundary_string(self, boundary) -> str:
        pass

    @abstractproperty
    def forcing_digit(self):
        """ bctype int value """
        raise NotImplementedError


Bctype = BoundaryForcing
