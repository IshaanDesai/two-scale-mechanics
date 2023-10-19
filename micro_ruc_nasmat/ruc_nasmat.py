"""
Solve RVE using NASMAT: https://software.nasa.gov/software/LEW-20244-1
"""
import numpy as np
from ctypes import cdll, byref, c_double, POINTER


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        """
        self._dims = 3
        self._sim_id = sim_id

        self._nasmat = cdll.LoadLibrary('./np_emulator.so')

    def solve(self, strains, dt):
        assert dt != 0

        strain_values = np.zeros((6))
        for i in range(3):
            strain_values[i] = strains["strains1to3"][i]
            strain_values[i + 3] = strains["strains4to6"][i]

        strain_values_ = strains.ctypes.data_as(POINTER(c_double))

        stresses = np.arange(0, 6, dtype=np.float64)
        stresses_ = stresses.ctypes.data_as(POINTER(c_double))

        self._nasmat.__mod_nasmat_MOD_nasmat3d_solve.argtypes = [POINTER(c_double)]
        self._nasmat.__mod_nasmat_MOD_nasmat3d_solve.restype = c_double

        self._nasmat.__mod_nasmat_MOD_nasmat3d_solve(strain_values_, stresses_)

        return {"stresses1to3": stresses[0:3], "stresses4to6": stresses[3:6]}
