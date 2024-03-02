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

        self.nasmat_lib = cdll.LoadLibrary('./lib_nasmat.so')
        self.homogenize = getattr(self.nasmat_lib, "__mod_precice_adapter_nasmat_MOD_precice_nasmat_homogenization")

    def solve(self, macro_data, dt):
        assert dt != 0

        strain_values = np.zeros((6))
        for i in range(3):
            strain_values[i] = macro_data["strains1to3"][i]
            strain_values[i + 3] = macro_data["strains4to6"][i]


        cmat = np.arange(0, 36, dtype=np.float64)
        cmat_ = cmat.ctypes.data_as(POINTER(c_double))
        self.homogenize.restype = c_double
        self.homogenize( cmat_)

        stresses = np.zeros((6))
        cmat = np.ctypeslib.as_array(cmat_, shape=(6, 6))
        stresses = np.dot(cmat, strain_values)
        return {"stresses1to3": stresses[0:3], "stresses4to6": stresses[3:6]}



def main():
    sim = MicroSimulation(1)
    strains = {"strains1to3": [0.1, 0.2, 0.3], "strains4to6": [0.4, 0.5, 0.6]}
    dt = 0.1
    print(sim.solve(strains, dt))

if __name__ == "__main__":
    main()
