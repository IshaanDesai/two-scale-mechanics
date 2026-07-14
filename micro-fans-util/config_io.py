import json
import copy
from typing import Optional


def write_json(data, filename):
    with open(filename, "w") as outfile:
        json.dump(data, outfile, indent=4)


def load_json(filename):
    with open(filename, "r") as infile:
        data = json.load(infile)
    return data


def add_ADA(config):
    config_ada = config
    config_ada["simulation_params"]["adaptivity"] = True
    config_ada["simulation_params"]["adaptivity_settings"] = {
        "type": "global",
        "data": config_ada["coupling_params"]["read_data_names"],
        "history_param": 0.5,
        "coarsening_constant": 0.3,
        "refining_constant": 0.4,
        "every_implicit_iteration": True,
        "output_type": "all",
        "output_n": 1,
        "mappings": [
            {
                "src_fields": config_ada["coupling_params"]["read_data_names"],
                "dst_fields": config_ada["coupling_params"]["write_data_names"],
                "interp_id": "adaptivity",
            }
        ],
    }
    config_ada["simulation_params"]["interpolation_configs"] = [
        {
            "type": "RBF",
            "id": "adaptivity",
            "rbf_config": {"basis": {"type": "c6"}, "n_neighbors": 50},
            "domain_config": {
                "max_filling": 8,
                "coarsening_factor": 2,
                "projection": {"type": "std", "target_dims": 3},
            },
        }
    ]
    return config_ada


def add_MADA(config, dim_mada: Optional[int] = None):
    micro_names = ["PyFANS0", "PyFANS1", "PyFANS2"]
    if dim_mada is None:
        dim_mada = 3
    if dim_mada < 2 or dim_mada > 3:
        raise ValueError("dim_mada must be 2 or 3")

    config_mada = config
    config_mada["simulation_params"]["model_adaptivity"] = True
    config_mada["simulation_params"]["model_adaptivity_settings"] = {
        "switching_function": "mada_switcher"
    }
    config_mada["micro_file_names"] = [micro_names[idx] for idx in range(dim_mada)]
    return config_mada


def add_stateless(config, has_mada):
    dim = len(config["micro_file_names"])
    config["micro_stateless_flags"] = [True] * dim
    return config


def gen_config(
    num_mm_ranks,
    num_workers,
    use_slurm,
    mpi_impl,
    decomp_dim,
    target_configs,
    dim_mada: Optional[int] = None,
):
    base_config = load_json("micro-manager-pyfans-config.json")
    base_config["simulation_params"]["decomposition"][decomp_dim] = num_mm_ranks
    base_config["tasking"]["num_workers"] = num_workers
    base_config["tasking"]["is_slurm"] = use_slurm
    base_config["tasking"]["mpi_impl"] = mpi_impl

    for use_ada, use_mada, use_stateless in target_configs:
        config = copy.deepcopy(base_config)
        if use_ada:
            add_ADA(config)
        if use_mada:
            add_MADA(config, dim_mada)
        if use_stateless:
            add_stateless(config, use_mada)

        file_name = "micro-manager-pyfans-config"
        if use_ada:
            file_name += "-ada"
        if use_mada:
            if dim_mada is not None:
                file_name += f"-mada{dim_mada}"
            else:
                file_name += "-mada"
        if use_stateless:
            file_name += "-stateless"
        file_name += ".json"
        write_json(config, file_name)


USE_ADA = True
NO_ADA = False
USE_MADA = True
NO_MADA = False
USE_STATELESS = True
NO_STATELESS = False

MPI_INTEL = "intel"
MPI_OPEN = "open"

USE_SLURM = True
NO_SLURM = False
