config = {
    "microstructure": {
        "filepath": "sphere32.h5",
        "datasetname": "/sphere/32x32x32/ms",
        "L": [1.0, 1.0, 1.0]
    },

    "problem_type": "mechanical",
    "strain_type": "large",
    "materials": [
        {
            "phases": [0,1],
            "matmodel": "CompressibleNeoHookean",
            "material_properties": {
                "bulk_modulus": [62.5000, 222.222],
                "shear_modulus": [28.8462, 166.6667]
            }
        }
    ],

    "FE_type": "HEX8",
    "method": "cg",
    "error_parameters":{
        "measure": "Linfinity",
        "type": "absolute",
        "tolerance": 1e-10
    },
    "n_it": 100,
    "macroscale_loading": [
        [[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]]
    ],
    "results": ["stress_average", "strain_average", "absolute_error", "homogenized_tangent",
                "microstructure", "displacement", "displacement_fluctuation", "stress", "strain"]
}

import os
if not os.path.exists("./sphere32.h5"):
    raise IOError("File 'sphere32.h5' does not exist")

import json
with open("input_large.json", "w") as json_file:
    json.dump(config, json_file)

import subprocess
proc = subprocess.run("mpiexec -n 8 FANS input_large.json output_large.h5", shell=True)
proc.check_returncode()

import h5py
with h5py.File("output_large.h5", "r") as h5file:
    tan = h5file["/sphere/32x32x32/ms_results/load0/time_step0/homogenized_tangent"][...]

model_input = {
    "zero_tangent": tan.tolist(),
}

with open("large_lin_np_input.json", "w") as json_file:
    json.dump(model_input, json_file)

os.remove("input_large.json")
os.remove("sphere32.h5")
os.remove("output_large.h5")