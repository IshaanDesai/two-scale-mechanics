import subprocess
import os
import json
import h5py
import numpy as np

from matplotlib import pyplot as plt

def gen_input_compressible():
    data = {
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
        "macroscale_loading": [],
        "results": ["stress_average", "strain_average"]
    }
    return data

def gen_input_linear():
    data = {
        "microstructure": {
            "filepath": "sphere32.h5",
            "datasetname": "/sphere/32x32x32/ms",
            "L": [1.0, 1.0, 1.0]
        },
        "problem_type": "mechanical",
        "matmodel": "LinearElasticIsotropic",
        "material_properties":{
            "bulk_modulus": [62.5000, 222.222],
            "shear_modulus": [28.8462, 166.6667]
        },
        "method": "cg",
        "error_parameters":{
            "measure": "Linfinity",
            "type": "absolute",
            "tolerance": 1e-10
        },
        "n_it": 100,
        "macroscale_loading":   [],
        "results": ["stress_average", "strain_average"]
    }
    return data

def attach_loadings(input_data, num_loadings):
    is_small_strain = True
    if "strain_type" in input_data:
        is_small_strain = input_data["strain_type"] == "small"

    loadings = [None] * num_loadings
    main_loads = np.linspace(1e-0, 1e+2, num_loadings)
    for i in range(num_loadings):
        if is_small_strain: loadings[i] = [[main_loads[i], 0.0, 0.0, 0.0, 0.0, 0.0]]
        else:
            step_data = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5]
            steps = [None] * len(step_data)
            for j in range(len(step_data)):
                step_load = np.zeros(9)
                step_load[0] = main_loads[i]
                step_load[4] = 1.0
                step_load[8] = 1.0
                step_load[1] = step_data[j]
                steps[j] = step_load.tolist()
            loadings[i] = steps
    input_data["macroscale_loading"] = loadings

def write_input(input_data):
    with open("input.json", "w") as f:
        json.dump(input_data, f, indent=4)





#if not os.path.exists('./input.json'):
#    print("Please place desired input file in this directory. Then run again.")
input_data = gen_input_compressible()
attach_loadings(input_data, 100)
write_input(input_data)

if not os.path.exists('./output.h5'):
    subprocess.run(
        "mpiexec -n 8 FANS ./input.json ./output.h5", shell=True
    )

with open("./input.json") as f:
    micro_config = json.load(f)
dset = micro_config["microstructure"]["datasetname"]
result_path = f"{dset}_results"

with h5py.File("./output.h5", "r") as f:
    num_loads = len(f[result_path].keys())

    loads = [None] * num_loads
    for load_id in range(num_loads):
        step_ds = f[f"{result_path}/load{load_id}"]
        num_steps = len(step_ds.keys())
        steps = [None] * num_steps
        for step_id in range(num_steps):
            result_ds = step_ds[f"time_step{step_id}"]
            results = {k:result_ds[k][...] for k in result_ds.keys()}
            steps[step_id] = results
        loads[load_id] = steps

plt.figure(figsize=(12, 8))
for load in loads:
    end_step = load[-1]
    stress = np.array(end_step["stress_average"])
    strain = np.array(end_step["strain_average"])
    plt.scatter(strain[0], stress[0])
plt.yscale("log")
plt.xscale("log")
plt.show()
