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
        found_lin_np = os.path.exists("./large_lin_np_input.json")
        inp_path = None
        if found_input1: inp_path = "./input0.json"
        elif found_lin_np: inp_path = "./large_lin_np_input.json"
        else: raise RuntimeError("No input file found")

        with open(inp_path, "r") as infile:
            data = json.load(infile)
        self._tan = np.asarray(data['zero_tangent'])

    def solve(self, macro_data, dt):
        eps = np.zeros(9)
        eps[0:3] = macro_data["strains1to3"]
        eps[3:6] = macro_data["strains4to6"]
        eps[6:9] = macro_data["strains7to9"]

        sig = np.matmul(self._tan, eps)
        if np.allclose(sig, 0): sig += 1e-16 * np.random.uniform(0, 1, 9)

        t = self._tan
        result = {
            "stresses1to3" : sig[0:3].tolist(),
            "stresses4to6" : sig[3:6].tolist(),
            "stresses7to9" : sig[6:9].tolist(),
            "cmat1" : self._tan[0, 0:3].tolist(),
            "cmat2" : self._tan[0, 3:6].tolist(),
            "cmat3" : self._tan[0, 6:9].tolist(),
            "cmat4" : self._tan[1, 0:3].tolist(),
            "cmat5" : self._tan[1, 3:6].tolist(),
            "cmat6" : self._tan[1, 6:9].tolist(),
            "cmat7" : self._tan[2, 0:3].tolist(),
            "cmat8" : self._tan[2, 3:6].tolist(),
            "cmat9" : self._tan[2, 6:9].tolist(),
            "cmat10" : self._tan[3, 0:3].tolist(),
            "cmat11" : self._tan[3, 3:6].tolist(),
            "cmat12" : self._tan[3, 6:9].tolist(),
            "cmat13" : self._tan[4, 0:3].tolist(),
            "cmat14" : self._tan[4, 3:6].tolist(),
            "cmat15" : self._tan[4, 6:9].tolist(),
            "cmat16" : self._tan[5, 0:3].tolist(),
            "cmat17" : self._tan[5, 3:6].tolist(),
            "cmat18" : self._tan[5, 6:9].tolist(),
            "cmat19" : self._tan[6, 0:3].tolist(),
            "cmat20" : self._tan[6, 3:6].tolist(),
            "cmat21" : self._tan[6, 6:9].tolist(),
            "cmat22" : self._tan[7, 0:3].tolist(),
            "cmat23" : self._tan[7, 3:6].tolist(),
            "cmat24" : self._tan[7, 6:9].tolist(),
            "cmat25" : self._tan[8, 0:3].tolist(),
            "cmat26" : self._tan[8, 3:6].tolist(),
            "cmat27" : self._tan[8, 6:9].tolist(),
        }
        return result

    def set_state(self, state):
        self._tan = state["tan"]
        self._sim_id = state["gid"]

    def get_state(self):
        return {"gid": self._sim_id, "tan": self._tan}

    def get_global_id(self):
        return self._sim_id