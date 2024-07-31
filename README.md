# composite-multiscale

Two-scale coupled simulation of a composite structure using the preCICE coupling library. One meso-scale simulation is coupled to many micro-scale simulations. Both the scales are solved using a range of solvers.

## Setup

The meso-scale model is a 3D beam structure which is being axially loaded. The micro-scale model is a 3D single fibre structure.

## Solvers

### Meso scale

- Dummy solver written in Python (only for testing)
- [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/)

### Micro scale

- [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/)
- [FANS](https://github.tik.uni-stuttgart.de/DAE/FANS)

#### Fourier Accelerated Nodal Solvers (FANS)

##### Install FANS

Get the FANS source code from the University of Stuttgart [FANS repository](https://github.tik.uni-stuttgart.de/DAE/FANS) and move into the `FANS` directory.

```bash
git clone https://github.tik.uni-stuttgart.de/DAE/FANS.git
```

and switch to the `micro_manager_fans` branch

See the FANS [README](FANS/README.md) for more information on the requirements of FANS.

To build and install FANS, run

```bash
mkdir build
cd build
cmake ..
sudo make install
```

##### Build the FANS python bindings

Move into the `composite-multiscale` directory and run

```bash
cd micro_FANS/python
mkdir build
cd build
cmake ..
cmake --build . -j
```

##### Configure the FANS simulation for the Micro Manager

FANS is configured in a JSON configuration file. The path to said file is read from a `input.json` file.

The name `input.json` is fixed and it contains the path to the input file and output file path.

```json
{
    "output_path": "result.h5",
    "input_file": "input_files/test_MechLinear.json"
}
```

The output path is required by FANS but not used in the preCICE coupling.
The input path is the path to the FANS input file where FANS internal properties like material parameters of the micro structure and the path to the micro structure are set.

##### Run the FANS simulation

There are two shared libraries available as FANS python bindings that can be used to run FANS with preCICE:

- `MicroFANS` for mechanical problems
- `MicroFANSTHERMAL` for thermal problems

```bash
cd ../../test
micro-manager-precice micro-manager-config-mech.json
```

## Dependencies

Apart from the dependencies of the solvers, the case setup itself relies on the following software:

- preCICE: See the scripts to [build_precice with MPI](build_scripts/build_precice_v3_with_Intel_MPI.sh) and [wihtout MPI](build_scripts/build_precice_v3_without_MPI.sh).
- Micro Manager: See the script to [build the Micro Manager](build_scripts/build_micro_manager.sh).

## Running the simulation

Files relevant to solvers on the meso scale are located in the folders with names `meso_*`, and similar solvers for micro scale are located in folders with names `micro_*`. Each folder has a run script, to be run directly to start the solver.

## Running simulations on a cluster

The ABAQUS-ABAQUS case was originally designed to be run on the [Great Lakes HPC cluster](https://arc.umich.edu/greatlakes/) at the University of Michigan. All job scripts are intended to work on any cluster with a SLURM system. Be aware that minor modifications may be necessary to the job scripts to be made usable on a cluster.
