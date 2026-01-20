import numpy as np

from .util import Registry

class Dimensions:
    """
    Container for mesh dimension information.

    Stores the lower corner coordinates, lengths, and number of elements
    in each spatial direction.

    Parameters
    ----------
    low : array_like
        Lower corner coordinates [x0, y0, z0].
    lengths : array_like
        Lengths in each direction [lx, ly, lz].
    elements : array_like
        Number of elements in each direction [nx, ny, nz].

    Attributes
    ----------
    low : array_like
        Lower corner coordinates.
    lengths : array_like
        Lengths in each direction.
    elements : array_like
        Number of elements in each direction.
    """
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
        """
        Get all dimension parameters as a tuple.

        Returns
        -------
        tuple
            (x0, y0, z0, lx, ly, lz, nx, ny, nz)
        """
        return self.x0, self.y0, self.z0, self.lx, self.ly, self.lz, self.nx, self.ny, self.nz


class Locators:
    """
    Collection of locator functions for identifying boundary facets.

    Provides predefined locator functions for identifying facets on
    axis-aligned planes of a normalized unit cube domain.

    Attributes
    ----------
    DIRICLET_TAG : int
        Tag value for Dirichlet boundary facets (2).
    NEUMANN_TAG : int
        Tag value for Neumann boundary facets (3).
    FUNCTIONS : Registry
        Registry of available locator functions.
    """
    DIRICLET_TAG :int = 2
    NEUMANN_TAG  :int = 3

    FUNCTIONS :Registry = Registry()

    @staticmethod
    @FUNCTIONS.register
    def plane_xy_low(u):
        """
        Locate facets on the low XY plane (z=0).

        Parameters
        ----------
        u : numpy.ndarray
            Array of coordinates with shape (3, n_points).

        Returns
        -------
        numpy.ndarray
            Boolean array indicating which points are on the plane.
        """
        return np.isclose(u[2], 0)
    @staticmethod
    @FUNCTIONS.register
    def plane_xz_low(u):
        """
        Locate facets on the low XZ plane (y=0).

        Parameters
        ----------
        u : numpy.ndarray
            Array of coordinates with shape (3, n_points).

        Returns
        -------
        numpy.ndarray
            Boolean array indicating which points are on the plane.
        """
        return np.isclose(u[1], 0)
    @staticmethod
    @FUNCTIONS.register
    def plane_yz_low(u):
        """
        Locate facets on the low YZ plane (x=0).

        Parameters
        ----------
        u : numpy.ndarray
            Array of coordinates with shape (3, n_points).

        Returns
        -------
        numpy.ndarray
            Boolean array indicating which points are on the plane.
        """
        return np.isclose(u[0], 0)
    @staticmethod
    @FUNCTIONS.register
    def plane_xy_high(u):
        """
        Locate facets on the high XY plane (z=1).

        Parameters
        ----------
        u : numpy.ndarray
            Array of coordinates with shape (3, n_points).

        Returns
        -------
        numpy.ndarray
            Boolean array indicating which points are on the plane.
        """
        return np.isclose(u[2], 1)
    @staticmethod
    @FUNCTIONS.register
    def plane_xz_high(u):
        """
        Locate facets on the high XZ plane (y=1).

        Parameters
        ----------
        u : numpy.ndarray
            Array of coordinates with shape (3, n_points).

        Returns
        -------
        numpy.ndarray
            Boolean array indicating which points are on the plane.
        """
        return np.isclose(u[1], 1)
    @staticmethod
    @FUNCTIONS.register
    def plane_yz_high(u):
        """
        Locate facets on the high YZ plane (x=1).

        Parameters
        ----------
        u : numpy.ndarray
            Array of coordinates with shape (3, n_points).

        Returns
        -------
        numpy.ndarray
            Boolean array indicating which points are on the plane.
        """
        return np.isclose(u[0], 1)

    @staticmethod
    def from_name(name):
        """
        Retrieve a locator function by name.

        Parameters
        ----------
        name : str
            Name of the locator function.

        Returns
        -------
        callable
            The locator function.
        """
        return Locators.FUNCTIONS.get_by_name(name)

    @staticmethod
    def is_name_valid(name):
        """
        Check if a locator name is valid.

        Parameters
        ----------
        name : str
            Name to check.

        Returns
        -------
        bool
            True if the name is a valid locator, False otherwise.
        """
        return Locators.FUNCTIONS.is_name_valid(name)