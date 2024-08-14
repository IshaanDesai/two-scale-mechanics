# Micro-scale simulation using NASMAT

This code solves a single fiber micro simulation using [NASMAT](https://software.nasa.gov/software/LEW-20244-1).

## Running

Run the Micro Manager in the following way

```bash
micro-manager-precice micro-manager-config.json
```

Run the Micro Manager in parallel in the following way

```bash
mpiexec -n <number-of-procs> micro-manager-precice micro-manager-config.json
```

by adding the number of processors you wish to solve with.
