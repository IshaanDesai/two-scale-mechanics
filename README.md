# Two-scale Mechanics Problems Experimental Setups

This repository contains experimental case setups for two-scale coupled simulations of mechanical systems using the [preCICE coupling library](https://precice.org/) and the [Micro Manager](https://precice.org/tooling-micro-manager-overview.html).

Each folder contains a case setup for either a meso-scale problem or a micro-scale problem solved with an existing simulation software.
The folder name follows the nomenclature `scale-solver-geometry`.
For example, `meso-ccx-notch` contains the case files for a meso-scale problem solved using CalculiX on a notch geometry.
A meso-scale problem and a micro-scale problem need to be executed from their respective folders simultaneously.
Each folder has a README with instructions on how to run the case.

## Available Case Setups

### Meso-scale participants

| Folder | Solver | Geometry |
|---|---|---|
| `meso-abaqus-crossply` | Abaqus | Cross-ply laminate |
| `meso-ccx-notch` | CalculiX | Notch |
| `meso-ccx-single-element` | CalculiX | Single element |
| `meso-dummy` | Python dummy | Generic (for testing) |
| `meso-fenics` | FEniCSx | Bar and notch |

### Micro-scale participants

| Folder | Solver | Geometry |
|---|---|---|
| `micro-abaqus-sphere` | Abaqus | Sphere RUC |
| `micro-fans-bar-sphere` | FANS | Bar with spherical microstructure |
| `micro-fans-bar-voro` | FANS | Bar with Voronoi microstructure |
| `micro-fans-notch-sphere` | FANS | Notch with spherical microstructure |
| `micro-fans-notch-voro` | FANS | Notch with Voronoi microstructure |
| `micro-fans-rivet` | FANS | Rivet |
| `micro-nasmat-sphere` | NASMAT | Sphere RUC |

### Utilities and scripts

- `fenics-fans-scripts/` — helper scripts and job configurations for running FEniCSx–FANS cases on HPC systems
- `micro-fans-surrogate/` — surrogate model scripts for FANS-based micro simulations
- `micro-fans-util/` — utility scripts for generating FANS model configurations and sphere meshes

### preCICE configuration files

The top-level `precice-config-*.xml` files are shared preCICE configurations used by multiple case setups:

- `precice-config-for-fans.xml` — configuration for FANS-based micro participants
- `precice-config-for-nasmat.xml` — configuration for NASMAT-based micro participants
- `precice-config-fans-small-strain.xml` — FANS configuration for small-strain problems
- `precice-config-fans-large-strain.xml` — FANS configuration for large-strain problems

## Dependencies

Running these two-scale cases requires

- preCICE: [installation instructions](https://precice.org/installation-overview.html)
- Micro Manager: [installation instructions](https://precice.org/tooling-micro-manager-installation.html)

Each setup individually requires the relevant solver software to be installed.
These details are specified in the README file in the respective folders.

## Cleaning up

The `clean-all.sh` script in the root of the repository iterates over all subfolders and calls the `clean.sh` script in each one to remove generated output files.

## Running on HPC systems

Large problems and scaling studies with these case setups have been conducted on [SuperMUC-NG](https://doku.lrz.de/supermuc-ng-10745965.html).
For instructions and guidelines for running problems on SuperMUC-NG, see [these notes](./SuperMUC-NG-Notes.md).

## Contributors

The two-scale problems collected in this repository are the result of several collaborations.
Contributors are

- [Ishaan Desai](https://github.com/IshaanDesai)
- [Minh Hoang Nguyen](https://github.com/mhoangnUM)
- [Ibrahim Kaleel](https://github.com/kalupaika)
- [Torben Schiz](https://github.com/tjwsch)
- [Alex Hocks](https://github.com/Snapex2409)
