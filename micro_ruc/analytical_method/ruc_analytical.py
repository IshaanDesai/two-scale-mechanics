"""
Python script to replicate the behavior of a RUC. Actually a micro simulation is not solved,
but instead a simple material model is resolved.
The material model is written by Minh Hoang Nguyen.
"""
import numpy as np


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        """
        self._dims = 3
        self._sim_id = sim_id

        # Parameters for material model written by Minh Hoang Nguyen
        self._E11 = 1.6e5
        self._E22 = 9.0e3
        self._E33 = self._E22
        self._nu12 = 0.32
        self._nu23 = 0.44
        self._nu13 = self._nu12
        self._G12 = 5.7e3
        self._G13 = self._G12
        self._G23  = self._E22 / 2.0 / (1.0 + self._nu23)
        self._nu31 = self._nu13 * self._E33 / self._E11
        self._nu21 = self._nu12 * self._E22 / self._E11
        self._nu31 = self._nu13 * self._E33 / self._E11
        self._nu32 = self._nu23 * self._E33 / self._E22

        self._Q = np.zeros((6,6))

    def solve(self, strains, dt):
        assert dt != 0

        strain_values = np.zeros((6))
        for i in range(3):
            strain_values[i] = strains["strains1to3"][i]
            strain_values[i + 3] = strains["strains4to6"][i]

        # Material model from VUMAT.f written by Minh Hoang Nguyen
        temp = self._nu12 * self._nu21 + self._nu23 * self._nu32 + self._nu13 * self._nu31 + 2.0 * self._nu21 * self._nu32 * self._nu13

        self._Q[:,:] = 0

        self._Q[0,0] = (1.0 - self._nu23 * self._nu32) * self._E11 / (1.0 - temp)
        self._Q[1,1] = (1.0 - self._nu13 * self._nu31) * self._E22 / (1.0 - temp)
        self._Q[2,2] = (1.0 - self._nu12 * self._nu21) * self._E33 / (1.0 - temp)
        self._Q[0,1] = (self._nu21 + self._nu31 * self._nu23) * self._E11 / (1.0 - temp)
        self._Q[1,0] = self._Q[0,1]
        self._Q[0,2] = (self._nu31 + self._nu21 * self._nu32) * self._E11 / (1.0 - temp)
        self._Q[2,0] = self._Q[0,2]
        self._Q[1,2] = (self._nu32 + self._nu12 * self._nu31) * self._E22 / (1.0 - temp)
        self._Q[2,1] = self._Q[1,2]
        self._Q[3,3] = self._G12
        self._Q[4,4] = self._G23
        self._Q[5,5] = self._G13

        stresses = np.matmul(self._Q, strain_values)

        return {"stresses1to3": stresses[0:3], "stresses4to6": stresses[3:6]}
