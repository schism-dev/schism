from abc import ABC, abstractmethod


class TidalDataProvider(ABC):

    @abstractmethod
    def get_elevation(self, constituent, vertices):
        """Returns amplitude and phase of constituent at given vertices."""

    @abstractmethod
    def get_velocity(self, constituent, vertices):
        """Returns amplitude and phase of constituent at given vertices."""

    @property
    @abstractmethod
    def constituents(self):
        """Returns list of constituents available in database."""
