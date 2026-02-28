import json
import numpy as np
import os

class MicroSimulation:
    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        """
        self._sim_id = None
        self._tan = None

        if sim_id >= 0:
            self._sim_id = sim_id
            self._load_input_file()

    def _load_input_file(self):
        found_input1 = os.path.exists("./input0.json")
        found_lin_np = os.path.exists("./lin_np_input.json")
        inp_path = None
        if found_input1: inp_path = "./input0.json"
        elif found_lin_np: inp_path = "./lin_np_input.json"
        else: raise RuntimeError("No input file found")

        with open(inp_path, "r") as infile:
            data = json.load(infile)
        self._tan = np.asarray(data['zero_tangent'])

    def solve(self, macro_data, dt):
        eps = np.zeros(6)
        eps[0:3] = macro_data["strains1to3"]
        eps[3:6] = macro_data["strains4to6"]

        sig = np.matmul(self._tan, eps)
        if np.allclose(sig, 0): sig += 1e-16 * np.random.uniform(0, 1, 6)

        t = self._tan
        result = {
            "stresses1to3" : sig[0:3].tolist(),
            "stresses4to6" : sig[3:6].tolist(),
            "cmat1"        : [t[0, 0], t[0, 1], t[0, 2]],
            "cmat2"        : [t[0, 3], t[0, 4], t[0, 5]],
            "cmat3"        : [t[1, 1], t[1, 2], t[1, 3]],
            "cmat4"        : [t[1, 4], t[1, 5], t[2, 2]],
            "cmat5"        : [t[2, 3], t[2, 4], t[2, 5]],
            "cmat6"        : [t[3, 3], t[3, 4], t[3, 5]],
            "cmat7"        : [t[4, 4], t[4, 5], t[5, 5]],
        }
        return result

    def set_state(self, state):
        self._tan = state["tan"]
        self._sim_id = state["gid"]

    def get_state(self):
        return {"gid": self._sim_id, "tan": self._tan}

    def get_global_id(self):
        return self._sim_id