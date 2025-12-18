# composite-multiscale

Case setups of two-scale coupled simulations of mechanical systems using the [preCICE coupling library](https://precice.org/) and the [Micro Manager](https://precice.org/tooling-micro-manager-overview.html).

Simulation setup of each case is in a separate folder, from which it is to be executed. The nomenclature of the folder name is *scale_solver_geometry*. For example *meso_ccx_notch* means a meso-scale simulation with CalculiX using a notch geometry.

## Dependencies

All setups require

- [preCICE](https://precice.org/installation-overview.html)
- [Micro Manager](https://precice.org/tooling-micro-manager-installation.html)

Each setup individually requires the relevant software to be installed. These details are specified in the README file in the respective folders.

## Available software

### Meso scale

- Dummy solver written in Python (only for testing)
- [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/)
- [CalculiX](https://www.calculix.de/)

### Micro scale

- [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/)
- [FANS](https://github.com/DataAnalyticsEngineering/FANS)
- [NASMAT](https://software.nasa.gov/software/LEW-20244-1)
