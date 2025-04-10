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

The stress tensor is

$$
\begin{pmatrix}
\sigma_{xx} & \sigma_{xy} & \sigma_{xz}\\
\sigma_{yx} & \sigma_{yy} & \sigma_{yz}\\
\sigma_{zx} & \sigma_{zy} & \sigma_{zz}\\
\end{pmatrix}
$$

The strain tensor is

$$
\begin{pmatrix}
\varepsilon_{xx} & \varepsilon_{xy} & \varepsilon_{xz}\\
\varepsilon_{yx} & \varepsilon_{yy} & \varepsilon_{yz}\\
\varepsilon_{zx} & \varepsilon_{zy} & \varepsilon_{zz}\\
\end{pmatrix}
$$

The stiffness matrix is

$$
\begin{pmatrix}
C_{xxxx} & C_{xxyy} & C_{xxzz} & C_{xxyz} & C_{xxxz} & C_{xxxy}\\
C_{xxyy} & C_{yyyy} & C_{yyzz} & C_{yyyz} & C_{yyxz} & C_{yyxy}\\
C_{xxzz} & C_{yyzz} & C_{zzzz} & C_{zzyz} & C_{zzxz} & C_{zzxy}\\
C_{xxyz} & C_{yyyz} & C_{zzyz} & C_{yzyz} & C_{yzxz} & C_{yzxy}\\
C_{xxxz} & C_{yyxz} & C_{zzxz} & C_{yzxz} & C_{xzxz} & C_{xzxy}\\
C_{xxxy} & C_{yyxy} & C_{zzxy} & C_{yzxy} & C_{xzxy} & C_{xyxy}\\
\end{pmatrix}
$$

### Voigt notation

Stress tensor represented as

$$\sigma = (\sigma_{xx},\sigma_{yy},\sigma_{zz},\sigma_{yz},\sigma_{xz},\sigma_{xy})$$

The vector form of the tensor is oftentimes represented with numeric indices

$$\sigma = (\sigma_{1},\sigma_{2},\sigma_{3},\sigma_{4},\sigma_{5},\sigma_{6})$$

Strain tensor represented as

$$\varepsilon = (\varepsilon_{xx},\varepsilon_{yy},\varepsilon_{zz},2\varepsilon_{yz},2\varepsilon_{xz},2\varepsilon_{xy})$$

The vector form of the tensor is oftentimes represented with numeric indices

$$\varepsilon = (\varepsilon_{1},\varepsilon_{2},\varepsilon_{3},\varepsilon_{4},\varepsilon_{5},\varepsilon_{6})$$

Stiffness matrix represented as

$$
C = (C_{xxxx}, C_{xxyy}, C_{xxzz}, C_{xxyz}, C_{xxxz}, C_{xxxy},\\
C_{yyyy}, C_{yyzz}, C_{yyyz}, C_{yyxz}, C_{yyxy},\\
C_{zzzz}, C_{xxyz}, C_{yyxz}, C_{zzxy},\\
C_{yzyz}, C_{yzxz}, C_{yzxy},\\
C_{xzxz}, C_{xzxy},\\
C_{xyxy})
$$

The vector form of the tensor is oftentimes represented with numeric indices

$$
C = (C_{11}, C_{12}, C_{13}, C_{14}, C_{15}, C_{16},\\
C_{22}, C_{23}, C_{24}, C_{25}, C_{26},\\
C_{33}, C_{34}, C_{35}, C_{36},\\
C_{44}, C_{45}, C_{46},\\
C_{55}, C_{56},\\
C_{66})
$$

### Mandel notation



### Coupling data structures:

- `stresses1to3`: $\sigma_{1}, \sigma_{2}, \sigma_{3}$
- `stresses4to6`: $\sigma_{4}, \sigma_{5}, \sigma_{6}$
- `strains1to3`: $\varepsilon_{1},\varepsilon_{2},\varepsilon_{3}$
- `strains4to6`: $\varepsilon_{4},\varepsilon_{5},\varepsilon_{6}$
- `cmat1`: $C_{11}, C_{12}, C_{13}$
- `cmat1`: $C_{14}, C_{15}, C_{16}$
- `cmat1`: $C_{22}, C_{23}, C_{24}$
- `cmat1`: $C_{25}, C_{26}, C_{33}$
- `cmat1`: $C_{34}, C_{35}, C_{36}$
- `cmat1`: $C_{44}, C_{45}, C_{46}$
- `cmat1`: $C_{55}, C_{56}, C_{66}$