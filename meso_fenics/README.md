# Meso-scale problem solved using FEniCSx

This model is able to solve small or large strain solid mechanics systems using FenicsX.
Simulations can run on arbitrary meshes in .msh format. Supported simulation types are pure mesoscopic,
pseudo coupled or coupled.

## Dependencies
- fenicsxprecice
- numpy  >1, <2
- mpi4py~=4.1.1
- petsc4py~=3.24.1
- fenics-basix~=0.10.0
- fenics-dolfinx~=0.10.0
- fenics-ufl~=2025.2.0
- gmsh~=4.15.0

## How to run

Run the following command with `meso_fenics` as the current working directory.
```
python macro.py path-to-config.json
```
For more information about implementational details or how to create a config file see the [docs](doc/Overview.md).

# Known Issues

- When providing a mesh in .msh format, physical groups need to be defined such that gmsh will create a facet_tags object upon loading.