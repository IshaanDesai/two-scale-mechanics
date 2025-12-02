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

The stress and strain tensors, and the stiffness matrix are symmetric for most mechanics problems. Therefore, all the components need not be communicated over the coupling. What needs to be handled is the notation used by both the solvers being coupled. For example, CalculiX uses the Voigt notation to represent the stress, strain and stiffness tensors. FANS uses Mandel notation.

### Voigt notation

Stress tensor is represented as a vector

$$
\begin{bmatrix}
\sigma_{xx} & \sigma_{yy} & \sigma_{zz} & \sigma_{yz} & \sigma_{zx} & \sigma_{xy}
\end{bmatrix}
$$

Strain tensor is represented as a vector

$$
\begin{bmatrix}
\varepsilon_{xx} & \varepsilon_{yy} & \varepsilon_{zz} & 2\varepsilon_{yz} & 2\varepsilon_{zx} & 2\varepsilon_{xy}
\end{bmatrix}
$$

Stiffness tensor is represented as a matrix

$$
\begin{bmatrix}
C_{(xx)(xx)} & C_{(xx)(yy)} & C_{(xx)(zz)} & C_{(xx)(yz)} & C_{(xx)(zx)} & C_{(xx)(xy)}\\
C_{(yy)(xx)} & C_{(yy)(yy)} & C_{(yy)(zz)} & C_{(yy)(yz)} & C_{(yy)(zx)} & C_{(yy)(xy)}\\
C_{(zz)(xx)} & C_{(zz)(yy)} & C_{(zz)(zz)} & C_{(zz)(yz)} & C_{(zz)(zx)} & C_{(zz)(xy)}\\
C_{(yz)(xx)} & C_{(zx)(yy)} & C_{(zx)(zz)} & C_{(yz)(yz)} & C_{(yz)(zx)} & C_{(yz)(xy)}\\
C_{(zx)(xx)} & C_{(zx)(yy)} & C_{(zx)(zz)} & C_{(zx)(yz)} & C_{(zx)(zx)} & C_{(zx)(xy)}\\
C_{(xy)(xx)} & C_{(xy)(yy)} & C_{(xy)(zz)} & C_{(xy)(yz)} & C_{(xy)(zx)} & C_{(xy)(xy)}\\
\end{bmatrix}
$$

### Mandel notation

Stress tensor represented as a vector

$$
\begin{bmatrix}
\sigma_{xx} & \sigma_{yy} & \sigma_{zz} & \sqrt2\sigma_{yz} & \sqrt2\sigma_{zx} & \sqrt2\sigma_{xy}
\end{bmatrix}
$$

Strain tensor represented as a vector

$$
\begin{bmatrix}
\varepsilon_{xx} & \varepsilon_{yy} & \varepsilon_{zz} & \sqrt2\varepsilon_{yz} & \sqrt2\varepsilon_{zx} & \sqrt2\varepsilon_{xy}
\end{bmatrix}
$$

Stiffness tensor represented as a matrix

$$
\begin{bmatrix}
C_{(xx)(xx)} & C_{(xx)(yy)} & C_{(xx)(zz)} & \sqrt2C_{(xx)(yz)} & \sqrt2C_{(xx)(zx)} & \sqrt2C_{(xx)(xy)}\\
C_{(yy)(xx)} & C_{(yy)(yy)} & C_{(yy)(zz)} & \sqrt2C_{(yy)(yz)} & \sqrt2C_{(yy)(zx)} & \sqrt2C_{(yy)(xy)}\\
C_{(zz)(xx)} & C_{(zz)(yy)} & C_{(zz)(zz)} & \sqrt2C_{(zz)(yz)} & \sqrt2C_{(zz)(zx)} & \sqrt2C_{(zz)(xy)}\\
\sqrt2C_{(yz)(xx)} & \sqrt2C_{(zx)(yy)} & \sqrt2C_{(zx)(zz)} & 2C_{(yz)(yz)} & 2C_{(yz)(zx)} & 2C_{(yz)(xy)}\\
\sqrt2C_{(zx)(xx)} & \sqrt2C_{(zx)(yy)} & \sqrt2C_{(zx)(zz)} & 2C_{(zx)(yz)} & 2C_{(zx)(zx)} & 2C_{(zx)(xy)}\\
\sqrt2C_{(xy)(xx)} & \sqrt2C_{(xy)(yy)} & \sqrt2C_{(xy)(zz)} & 2C_{(xy)(yz)} & 2C_{(xy)(zx)} & 2C_{(xy)(xy)}\\
\end{bmatrix}
$$

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
