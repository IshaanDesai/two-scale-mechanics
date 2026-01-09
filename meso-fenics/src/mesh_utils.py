import numpy as np

from .util import Registry

class Dimensions:
    def __init__(self, low, lengths, elements):
        self._low = low
        self._lengths = lengths
        self._elements = elements

    @property
    def low(self): return self._low
    @property
    def lengths(self): return self._lengths
    @property
    def elements(self): return self._elements
    @property
    def x0(self): return self._low[0]
    @property
    def y0(self): return self._low[1]
    @property
    def z0(self): return self._low[2]
    @property
    def lx(self): return self._lengths[0]
    @property
    def ly(self): return self._lengths[1]
    @property
    def lz(self): return self._lengths[2]
    @property
    def nx(self): return self._elements[0]
    @property
    def ny(self): return self._elements[1]
    @property
    def nz(self): return self._elements[2]

    def get(self):
        return self.x0, self.y0, self.z0, self.lx, self.ly, self.lz, self.nx, self.ny, self.nz


class Locators:
    DIRICLET_TAG :int = 2
    NEUMANN_TAG  :int = 3

    FUNCTIONS :Registry = Registry()

    @staticmethod
    @FUNCTIONS.register
    def plane_xy_low(u): return np.isclose(u[2], 0)
    @staticmethod
    @FUNCTIONS.register
    def plane_xz_low(u): return np.isclose(u[1], 0)
    @staticmethod
    @FUNCTIONS.register
    def plane_yz_low(u): return np.isclose(u[0], 0)
    @staticmethod
    @FUNCTIONS.register
    def plane_xy_high(u): return np.isclose(u[2], 1)
    @staticmethod
    @FUNCTIONS.register
    def plane_xz_high(u): return np.isclose(u[1], 1)
    @staticmethod
    @FUNCTIONS.register
    def plane_yz_high(u): return np.isclose(u[0], 1)

    @staticmethod
    def from_name(name): return Locators.FUNCTIONS.get_by_name(name)

    @staticmethod
    def is_name_valid(name): return Locators.FUNCTIONS.is_name_valid(name)