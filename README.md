# composite-multiscale

Two-scale coupled simulation of a composite structure using the preCICE coupling library. One meso-scale simulation is coupled to many micro-scale simulations. Both the scales are solved using a range of solvers.

## Setup

The meso-scale model is a 3D beam structure which is being axially loaded. The micro-scale model is a 3D single fibre structure.

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

The stress tensor is

$$
\begin{bmatrix}
\sigma_{xx} & \sigma_{xy} & \sigma_{xz}\\
\sigma_{yx} & \sigma_{yy} & \sigma_{yz}\\
\sigma_{zx} & \sigma_{zy} & \sigma_{zz}\\
\end{bmatrix}
$$

Due to symmetry, the values of interest are

$$
\begin{bmatrix}
\sigma_{xx} & \sigma_{xy} & \sigma_{xz}\\
. & \sigma_{yy} & \sigma_{yz}\\
. & . & \sigma_{zz}\\
\end{bmatrix}
$$

The strain tensor is

$$
\begin{bmatrix}
\varepsilon_{xx} & \varepsilon_{xy} & \varepsilon_{xz}\\
\varepsilon_{yx} & \varepsilon_{yy} & \varepsilon_{yz}\\
\varepsilon_{zx} & \varepsilon_{zy} & \varepsilon_{zz}\\
\end{bmatrix}
$$

Due to symmetry, the values of interest are

$$
\begin{bmatrix}
\varepsilon_{xx} & \varepsilon_{xy} & \varepsilon_{xz}\\
. & \varepsilon_{yy} & \varepsilon_{yz}\\
. & . & \varepsilon_{zz}\\
\end{bmatrix}
$$

The stiffness matrix is

$$
\begin{bmatrix}
C_{xxxx} & C_{xxyy} & C_{xxzz} & C_{xxyz} & C_{xxxz} & C_{xxxy}\\
C_{xxyy} & C_{yyyy} & C_{yyzz} & C_{yyyz} & C_{yyxz} & C_{yyxy}\\
C_{xxzz} & C_{yyzz} & C_{zzzz} & C_{zzyz} & C_{zzxz} & C_{zzxy}\\
C_{xxyz} & C_{yyyz} & C_{zzyz} & C_{yzyz} & C_{yzxz} & C_{yzxy}\\
C_{xxxz} & C_{yyxz} & C_{zzxz} & C_{yzxz} & C_{xzxz} & C_{xzxy}\\
C_{xxxy} & C_{yyxy} & C_{zzxy} & C_{yzxy} & C_{xzxy} & C_{xyxy}\\
\end{bmatrix}
$$

Due to symmetry, the values of interest are

$$
\begin{bmatrix}
C_{xxxx} & C_{xxyy} & C_{xxzz} & C_{xxyz} & C_{xxxz} & C_{xxxy}\\
. & C_{yyyy} & C_{yyzz} & C_{yyyz} & C_{yyxz} & C_{yyxy}\\
. & . & C_{zzzz} & C_{zzyz} & C_{zzxz} & C_{zzxy}\\
. & . & . & C_{yzyz} & C_{yzxz} & C_{yzxy}\\
. & . & . & . & C_{xzxz} & C_{xzxy}\\
. & . & . & . & . & C_{xyxy}\\
\end{bmatrix}
$$

### Voigt notation

Stress tensor is represented as

$$
\begin{bmatrix}
\sigma_{xx} & \sigma_{xy} & \sigma_{xz}\\
. & \sigma_{yy} & \sigma_{yz}\\
. & . & \sigma_{zz}\\
\end{bmatrix}
$$

Strain tensor is represented as

$$
\begin{bmatrix}
\varepsilon_{xx} & 2\varepsilon_{xy} & 2\varepsilon_{xz}\\
. & \varepsilon_{yy} & 2\varepsilon_{yz}\\
. & . & \varepsilon_{zz}\\
\end{bmatrix}
$$

Stiffness matrix is represented as

$$
\begin{bmatrix}
C_{xxxx} & C_{xxyy} & C_{xxzz} & C_{xxyz} & C_{xxxz} & C_{xxxy}\\
. & C_{yyyy} & C_{yyzz} & C_{yyyz} & C_{yyxz} & C_{yyxy}\\
. & . & C_{zzzz} & C_{zzyz} & C_{zzxz} & C_{zzxy}\\
. & . & . & C_{yzyz} & C_{yzxz} & C_{yzxy}\\
. & . & . & . & C_{xzxz} & C_{xzxy}\\
. & . & . & . & . & C_{xyxy}\\
\end{bmatrix}
$$

$$
\begin{bmatrix}
C_{11} & C_{12} & C_{13} & C_{14} & C_{15} & C_{16}\\
. & C_{22} & C_{23} & C_{24} & C_{25} & C_{26}\\
. & . & C_{33} & C_{34} & C_{35} & C_{36}\\
. & . & . & C_{44} & C_{45} & C_{46}\\
. & . & . & . & C_{55} & C_{56}\\
. & . & . & . & . & C_{66}\\
\end{bmatrix}
$$

The matrix values are represented in a single vector.

### Mandel notation

Stress tensor represented as

$$
\begin{bmatrix}
\sigma_{xx} & \sqrt2\sigma_{xy} & \sqrt2\sigma_{xz}\\
. & \sigma_{yy} & \sqrt2\sigma_{yz}\\
. & . & \sigma_{zz}\\
\end{bmatrix}
$$

Strain tensor represented as

$$
\begin{bmatrix}
\varepsilon_{xx} & \sqrt2\varepsilon_{xy} & \sqrt2\varepsilon_{xz}\\
. & \varepsilon_{yy} & \sqrt2\varepsilon_{yz}\\
. & . & \varepsilon_{zz}\\
\end{bmatrix}
$$

Stiffness matrix represented as

$$
\begin{bmatrix}
C_{xxxx} & C_{xxyy} & C_{xxzz} & \sqrt2C_{xxyz} & \sqrt2C_{xxxz} & \sqrt2C_{xxxy}\\
. & C_{yyyy} & C_{yyzz} & \sqrt2C_{yyyz} & \sqrt2C_{yyxz} & \sqrt2C_{yyxy}\\
. & . & C_{zzzz} & \sqrt2C_{zzyz} & \sqrt2C_{zzxz} & \sqrt2C_{zzxy}\\
. & . & . & 2C_{yzyz} & 2C_{yzxz} & 2C_{yzxy}\\
. & . & . & . & 2C_{xzxz} & 2C_{xzxy}\\
. & . & . & . & . & 2C_{xyxy}\\
\end{bmatrix}
$$

### Coupling data structures

The coupling variables correspond to the following numeric notation specific quantities:

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

The user needs to take care that either both the macro and micro scale simulations use the same notation, or the change of notation is accounted for on one side. For example, CalculiX uses the Voigt notation, but FANS uses the Mandel notation.
