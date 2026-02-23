import importlib as ipl
import json
import os
import sys
import numpy as np
import shutil

sys.path.append(os.getcwd())

def write_json(data, filename):
    with open(filename, "w") as outfile:
        json.dump(data, outfile, indent=4)

def load_json(filename):
    with open(filename, "r") as infile:
        data = json.load(infile)
    return data

def load_backend_class(path_to_micro_file):
    CLS_NAME = "MicroSimulation"
    return getattr(ipl.import_module(path_to_micro_file, CLS_NAME), CLS_NAME)

def compute_zero_load_tangent():
    from mpi4py import MPI

    cls = load_backend_class("PyFANS")
    zero_load_input = {
        "strains1to3": [0.0, 0.0, 0.0],
        "strains4to6": [0.0, 0.0, 0.0]
    }

    sim = cls(0)
    zero_load_output = sim.solve(zero_load_input, 0.0)
    MPI.COMM_WORLD.Barrier()

    if MPI.COMM_WORLD.Get_rank() == 0:
        tan_buffer = np.zeros((7, 3))
        for i in range(7):
            tan_buffer[i, :] = zero_load_output[f"cmat{i+1}"]
        broad_cast_array = np.array([
            [ 0,  1,  2,  3,  4,  5],
            [ 1,  6,  7,  8,  9, 10],
            [ 2,  7, 11, 12, 13, 14],
            [ 3,  8, 12, 15, 16, 17],
            [ 4,  9, 13, 16, 18, 19],
            [ 5, 10, 14, 17, 19, 20]
        ])

        tan = tan_buffer.flatten()[broad_cast_array].tolist()

        output = {"zero_tangent": tan}
        write_json(output, "./lin_np_input.json")
        shutil.copy("../micro-fans-util/model_lin_np.py", "./model_lin.py")

    MPI.COMM_WORLD.Barrier()



