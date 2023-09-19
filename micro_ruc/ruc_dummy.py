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

    def solve(self, macro_data, dt):
        assert dt != 0
        self._micro_vector_data = []
        self._micro_scalar_data = macro_data["macro-scalar-data"] + 1
        for d in range(self._dims):
            self._micro_vector_data.append(macro_data["macro-vector-data"][d] + 1)

        return {"micro-scalar-data": self._micro_scalar_data.copy(),
                "micro-vector-data": self._micro_vector_data.copy()}

    def set_state(self, state):
        self._state = state

    def get_state(self):
        return self._state
