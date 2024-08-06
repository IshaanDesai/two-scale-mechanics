# Fourier Accelerated Nodal Solvers (FANS)

## Install FANS

Get the FANS source code from the University of Stuttgart [FANS repository](https://github.tik.uni-stuttgart.de/DAE/FANS) and move into the `FANS` directory. Go to the FANS directory and switch to the `micro_manager_fans` branch. Build FANS by following the instructions in its [README](FANS/README.md).

## Build the FANS python bindings

Go to the `python/` folder in this repository

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

## Configure the FANS simulation for the Micro Manager

FANS is configured in a JSON configuration file. The path to said file is read from a `input.json` file.

The name `input.json` is fixed and it contains the path to the input file and output file path.

```json
{
    "output_path": "result.h5",
    "input_file": "input_files/test_MechLinear.json"
}
```

The output path is required by FANS but not used in the preCICE coupling. The input path is the path to the FANS input file where FANS internal properties like material parameters of the micro structure and the path to the micro structure are set.

## Run the FANS simulation

There are two shared libraries available as FANS python bindings that can be used to run FANS

- `MicroFANS` for mechanical problems. Corresponding Micro Manager configuration file is `micro-manager-config-mesh.json`.
- `MicroFANSTHERMAL` for thermal problems. Corresponding Micro Manager configuration file is `micro-manager-config-thermal.json`.

```bash
cd ../../test
micro-manager-precice micro-manager-config-mech.json
```

to run a mechanical simulation
