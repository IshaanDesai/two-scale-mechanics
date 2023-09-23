"""
Dummy Python script to replicate the behavior of a RUC which does nothing.
"""


class MicroSimulation:

    def __init__(self):
        """
        Constructor of MicroSimulation class.
        """
        self._dims = 3
        self._stresses = None
        self._strains = None
        self._state = None

    def solve(self, strains, dt):
        assert dt != 0
        self._strains = strains["strains"]
        for d in range(self._dims):
            self._stresses.append(self._strains[d] * 0.1)

        return {"stresses": self._stresses.copy()}

    def set_state(self, strains):
        self._strains = strains

    def get_state(self):
        return self._strains
