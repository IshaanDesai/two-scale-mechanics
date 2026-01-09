import h5py
import numpy as np
import json
import subprocess

def swap(arr, i, j):
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp

def gen_fans_input(strain_data):
    fans_input = dict()
    fans_input["microstructure"] = {
        "filepath": "pseudo-coupling/sphere32.h5",
        "datasetname": "/sphere/32x32x32/ms",
        "L": [1.0, 1.0, 1.0]
    }
    fans_input["problem_type"] = "mechanical"
    fans_input["strain_type"] = "small"
    fans_input["matmodel"] = "LinearElasticIsotropic"
    fans_input["material_properties"] = {
        "bulk_modulus": [62.5000, 222.222],
        "shear_modulus": [28.8462, 166.6667]
    }
    fans_input["method"] = "cg"
    fans_input["error_parameters"] = {
        "measure": "Linfinity",
        "type": "absolute",
        "tolerance": 1e-10
    }
    fans_input["n_it"] = 100
    fans_input["results"] = ["homogenized_tangent", "stress_average"]

    strain_data_view = strain_data.reshape(-1, 6)
    load_cases = []
    for idx in range(strain_data_view.shape[0]):
        arr = [strain_data_view[idx, i] for i in range(6)]
        swap(arr, 3, 5)
        load_cases.append([arr])
    fans_input["macroscale_loading"] = load_cases
    num_loads = len(load_cases)

    with open('./pseudo-coupling/fans-input.json', 'w') as f:
        json.dump(fans_input, f)

    return num_loads

def process_outputs(num_loads):
    stresses = np.zeros((num_loads, 6))
    tangents = np.zeros((num_loads, 6, 6))
    selection = np.array([
        [0, 1, 2, 3, 4, 5],
        [1, 6, 7, 8, 9, 10],
        [2, 7, 11, 12, 13, 14],
        [3, 8, 12, 15, 16, 17],
        [4, 9, 13, 16, 18, 19],
        [5, 10, 14, 17, 19, 20],
    ])

    with h5py.File(f"./pseudo-coupling/fans-output.h5", 'r') as h5:
        for i in range(num_loads):
            c_tensor = h5["sphere"]['32x32x32']['ms_results'][f"load{i}"]['time_step0']['homogenized_tangent'][...]
            stress = h5["sphere"]['32x32x32']['ms_results'][f"load{i}"]['time_step0']['stress_average'][...]
            swap(stress, 3, 5)

            c_buffer = np.zeros((21,))
            c_buffer[0:6] = c_tensor[0, :]
            c_buffer[6:11] = c_tensor[1, 1:]
            c_buffer[11:15] = c_tensor[2, 2:]
            c_buffer[15:18] = c_tensor[3, 3:]
            c_buffer[18:20] = c_tensor[4, 4:]
            c_buffer[20] = c_tensor[5, 5]
            swap(c_buffer, 3, 5)
            swap(c_buffer, 8, 10)
            swap(c_buffer, 12, 14)
            swap(c_buffer, 15, 20)
            swap(c_buffer, 16, 19)
            swap(c_buffer, 17, 18)

            tangents[i, :, :] = c_buffer[selection]
            stresses[i, :] = stress[:]

    return stresses, tangents

if __name__ == '__main__':
    with h5py.File("../output/bar_strain.h5", "r") as f:
        strain_data = f["strain_data"][:]

    num_loads = gen_fans_input(strain_data)

    subprocess.call(f"FANS ./pseudo-coupling/fans-input.json ./pseudo-coupling/fans-output.h5", shell=True)

    stresses, tangents = process_outputs(num_loads)

    with h5py.File("pseudo-input.h5", 'w') as f:
        f.create_dataset("stress_data", data=stresses.flatten())
        f.create_dataset("tan_data",    data=tangents.flatten())
