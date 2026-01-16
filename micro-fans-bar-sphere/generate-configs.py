import json
import copy


def write_json(data, filename):
    with open(filename, 'w') as outfile:
        json.dump(data, outfile, indent=4)


def load_json(filename):
    with open(filename, 'r') as infile:
        data = json.load(infile)
    return data


def add_ADA(config):
    config_ada = config
    config_ada['simulation_params']['adaptivity'] = True
    config_ada['simulation_params']['adaptivity_settings'] = {
        "type": "local",
        "data": config_ada['coupling_params']['read_data_names'],
        "history_param": 0.5,
        "coarsening_constant": 0.3,
        "refining_constant": 0.4,
        "every_implicit_iteration": True,
        "output_cpu_time": True
    }
    return config_ada

def add_MADA(config):
    config_mada = config
    config_mada['simulation_params']['model_adaptivity'] = True
    config_mada['simulation_params']['model_adaptivity_settings'] = {
        "micro_file_names": ["PyFANS", "PyFANS", "PyFANS"],
        "switching_function": "mada_switcher"
    }
    return config_mada

def add_stateless(config, has_mada):
    if has_mada:
        config['simulation_params']['model_adaptivity_settings']["micro_stateless"] = [True, True, True]
    else:
        config['micro_stateless'] = True
    return config

if __name__ == '__main__':
    num_mm_ranks = 1
    num_workers = 2

    base_config = load_json('micro-manager-pyfans-config.json')

    target_configs = [
        # ADA MADA Stateless
        (True, False, False),
        (True, True, False),
        (True, True, True),
        (True, False, True),
        (False, True, False),
        (False, True, True),
        (False, False, True),
    ]

    for use_ada, use_mada, use_stateless in target_configs:
        config = copy.deepcopy(base_config)
        if use_ada: add_ADA(config)
        if use_mada: add_MADA(config)
        if use_stateless: add_stateless(config, use_mada)

        file_name = 'micro-manager-pyfans-config'
        if use_ada: file_name += '-ada'
        if use_mada: file_name += '-mada'
        if use_stateless: file_name += '-stateless'
        file_name += '.json'
        write_json(config, file_name)
