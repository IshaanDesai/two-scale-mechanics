# composite-multiscale

Two-scale coupled simulation of a composite structure using the preCICE coupling library. One meso-scale simulation is coupled to many micro-scale simulations. Both the scales are solved using [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/).

## Setup

The meso-scale model is a 3D beam structure which is being axially loaded. The micro-scale model is a 3D single fibre structure.

## Solvers

### Meso scale

ABAQUS and Dummy *(more details coming soon)*

### Micro scale

ABAQUS and NASMAT *(more details coming soon)*

## Dependencies

* preCICE: See the scripts to [build_precice with MPI](build_scripts/build_precice_v3_with_Intel_MPI.sh) and [wihtout MPI](build_scripts/build_precice_v3_without_MPI.sh).

* Micro Manager: See the script to [build the Micro Manager](build_scripts/build_micro_manager.sh).

## Running the dummy simulation

Open two terminals. In the first terminal, go to `meso_dummy/` folder. Run the dummy solver:

```bash
python dummy.py
```

In the second terminal, run the Micro Manager:

```bash
./run_micromanager_nasmat.sh
```

## Running the 3D ply laminate simulation

This case was originally designed to be run on the [Great Lakes HPC cluster](https://arc.umich.edu/greatlakes/) at the University of Michigan. In principle the setup should work on any cluster which has access to ABAQUS licenses. To run the case, submit a job via the job script:

```bash
sbatch submit_job.sbat 
```

## Post-processing
