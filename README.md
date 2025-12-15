# composite-multiscale

Two-scale coupled simulation of a composite structure using the preCICE coupling library. One meso-scale simulation is coupled to many micro-scale simulations. Both the scales are solved using a range of solvers.

## Setups

### Single element case

TODO

### Cantilever beam

TODO

### Notch

TODO

## Solvers

### Meso scale

- Dummy solver written in Python (only for testing)
- [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/)
- [CalculiX](https://www.calculix.de/)

### Micro scale

- [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/)
- [FANS](https://github.com/DataAnalyticsEngineering/FANS)
- [NASMAT](https://software.nasa.gov/software/LEW-20244-1)

## Dependencies

Apart from the dependencies of the solvers, the case setup itself relies on the following software:

- preCICE: See the scripts to [build_precice with MPI](build_scripts/build_precice_v3_with_Intel_MPI.sh) and [without MPI](build_scripts/build_precice_v3_without_MPI.sh).
- Micro Manager: See the script to [build the Micro Manager](build_scripts/build_micro_manager.sh).

## Running the simulation

Files relevant to solvers on the meso scale are located in the folders with names `meso_*`, and similar solvers for micro scale are located in folders with names `micro_*`. Each folder has a README explaining the simulation setup and how to start the solver.

## Running simulations on a cluster

The ABAQUS-ABAQUS case was originally designed to be run on the [Great Lakes HPC cluster](https://arc.umich.edu/greatlakes/) at the University of Michigan. All job scripts are intended to work on any cluster with a SLURM system. Be aware that minor modifications may be necessary to the job scripts to be made usable on a cluster.

## Notations

The stress, strain, and stiffness tensor are represented in either the Voigt or the Mandel notation. Every solver uses either of the notation and its own ordering of elements. The README in each solver folder describes the respective notation and ordering.

### Coupling data structures

The stress and strain tensors are represented by $1\times6$ vectors, and the stiffness tensor is represented as a $6\times6$ matrix. The data structures used in the coupling are structured as follows

- `stresses1to3`: $\sigma_{1}, \sigma_{2}, \sigma_{3}$
- `stresses4to6`: $\sigma_{4}, \sigma_{5}, \sigma_{6}$
- `strains1to3`: $\varepsilon_{1},\varepsilon_{2},\varepsilon_{3}$
- `strains4to6`: $\varepsilon_{4},\varepsilon_{5},\varepsilon_{6}$
- `cmat1`: $C_{11}, C_{12}, C_{13}$
- `cmat2`: $C_{14}, C_{15}, C_{16}$
- `cmat3`: $C_{22}, C_{23}, C_{24}$
- `cmat4`: $C_{25}, C_{26}, C_{33}$
- `cmat5`: $C_{34}, C_{35}, C_{36}$
- `cmat6`: $C_{44}, C_{45}, C_{46}$
- `cmat7`: $C_{55}, C_{56}, C_{66}$

Different solvers can use different notations, for example, CalculiX uses the Voigt notation, but FANS uses the Mandel notation. The multiplying factors need to be adjusting either at the time of writing data or at the time of reading data.
