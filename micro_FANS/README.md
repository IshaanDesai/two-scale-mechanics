# Fourier Accelerated Nodal Solvers (FANS)

## Install FANS

Clone [FANS](https://github.com/DataAnalyticsEngineering/FANS) and build it as a library.

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

When FANS is built as a library, a file `libFANS.so` is generated. Point to this file in the `micro_file_name` entry in the `micro-manager-config.json` configuration.

Run the Micro Manager

```bash
micro-manager-precice micro-manager-config.json
```
