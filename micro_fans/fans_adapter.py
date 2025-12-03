"""
Micro simulation
In this script we solve a dummy micro problem to just show the working of the macro-micro coupling
"""
import json
import copy
import subprocess
import os
import h5py

class MicroSimulation:
    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        """
        self._sim_id = sim_id
        self._json_config_template = MicroSimulation.load_config_template("./fans_adapter_input.json")
        try:
            if not os.path.exists('./tmp'): os.makedirs('./tmp')
        except OSError:
            pass

        self._state = None

    @staticmethod
    def load_config_template(path):
        try:
            with open(path, "r") as config_file:
                return json.load(config_file)
        except FileNotFoundError:
            raise FileNotFoundError("Config file not found")

    def get_local_input_file_name(self):
        return f"tmp/fans_adapter_input_{self._sim_id}.json"

    def get_local_output_file_name(self):
        return f"tmp/fans_adapter_output_{self._sim_id}.h5"

    def write_local_input(self, strains):
        json_copy = copy.deepcopy(self._json_config_template)
        json_copy["macroscale_loading"][0][0] = strains

        try:
            with open(f"./{self.get_local_input_file_name()}", "w") as config_file:
                json.dump(json_copy, config_file)
        except FileNotFoundError:
            raise FileNotFoundError("Input filepath invalid")

    def load_local_output(self):
        stresses = [0] * 6
        tangents = [0] * 21

        c_tensor = None
        with h5py.File(f"./{self.get_local_output_file_name()}", 'r') as h5:
            c_tensor = h5["sphere"]['32x32x32']['ms_results']['load0']['time_step0']['homogenized_tangent'][...]
            stresses[:] = h5["sphere"]['32x32x32']['ms_results']['load0']['time_step0']['stress_average'][...]

        tangents[0  :6] = c_tensor[0,  :]
        tangents[6 :11] = c_tensor[1, 1:]
        tangents[11:15] = c_tensor[2, 2:]
        tangents[15:18] = c_tensor[3, 3:]
        tangents[18:20] = c_tensor[4, 4:]
        tangents[20   ] = c_tensor[5,  5]

        return stresses, tangents

    def solve(self, macro_data, dt):
        assert dt != 0

        strains = [0] * 6
        strains[0:3] = macro_data['strains1to3']
        #strains[3:6] = macro_data['strains4to6']
        # FANS uses inverted order for off-diag elements
        strains[3:6] = macro_data['strains4to6'][::-1]
        self.write_local_input(strains)
        subprocess.call(f"FANS ./{self.get_local_input_file_name()} ./{self.get_local_output_file_name()}", shell=True)
        stresses, tangents = self.load_local_output()

        return {
            "stresses1to3": stresses[0:3],
            "stresses4to6": stresses[3:6][::-1],
            "cmat1": tangents[0:3],
            "cmat2": tangents[3:6][::-1],
            "cmat3": [tangents[i] for i in [6, 7, 10]],
            "cmat4": [tangents[i] for i in [9, 8, 11]],
            "cmat5": tangents[12:15][::-1],
            "cmat6": [tangents[i] for i in [20, 19, 17]],
            "cmat7": [tangents[i] for i in [18, 16, 15]],
        }

    def set_state(self, state):
        self._state = state

    def get_state(self):
        return self._state
