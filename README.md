# Two-scale Mechanics Experimental Setups

Experimental case setups of two-scale coupled simulations of mechanical systems using the [preCICE coupling library](https://precice.org/) and the [Micro Manager](https://precice.org/tooling-micro-manager-overview.html).

Simulation setup of each case is in a separate folder, from which it is to be executed. The nomenclature of the folder name is *scale_solver_geometry*. For example *meso_ccx_notch* means a meso-scale simulation with CalculiX using a notch geometry.

## Dependencies

All setups require

- [preCICE](https://precice.org/installation-overview.html)
- [Micro Manager](https://precice.org/tooling-micro-manager-installation.html)

Each setup individually requires the relevant software to be installed. These details are specified in the README file in the respective folders.

For known issues and tips regarding running on SuperMUC-NG see [here](./SuperMUC-NG-Notes.md).