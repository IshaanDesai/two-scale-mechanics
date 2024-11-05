# Fourier Accelerated Nodal Solvers (FANS)

## Install FANS

Clone [FANS](https://github.com/DataAnalyticsEngineering/FANS) and build [pyFANS](https://github.com/DataAnalyticsEngineering/FANS/tree/develop/pyfans).

## Configure the FANS simulation for the Micro Manager

FANS is configured in a JSON configuration file. The name `input.json` is fixed. More information about this file is in the [FANS configuration documentation](https://github.com/DataAnalyticsEngineering/FANS?tab=readme-ov-file#input-file-format).

## Run the FANS simulation

Point to the `pyFANS.so` file in the `micro_file_name` entry in the `micro-manager-config.json` configuration.

Run the Micro Manager

```bash
micro-manager-precice micro-manager-config.json
```
