# Meso-scale problem solved using FEniCSx

Solver for small or large strain solid mechanics systems using FEniCSx. Simulations can run on arbitrary meshes in .msh format. Supported simulation types are pure mesoscopic, pseudo coupled or coupled.

## Dependencies

- fenicsxprecice (needs to be manually installed from source [fenicsx-adapter](https://github.com/precice/fenicsx-adapter))
- numpy  >1, <2
- mpi4py~=4.1.1
- petsc4py~=3.24.1
- fenics-basix~=0.10.0
- fenics-dolfinx~=0.10.0
- fenics-ufl~=2025.2.0
- gmsh~=4.15.0

## How to run

The `run.sh` script is used to run meso-scale problems. For example, to run the bar geometry, run:

```bash
./run.sh -problem bar
```

For more information about implementation details or how to create a config file see the [docs](doc/Overview.md).

## Known Issues

When providing a mesh in .msh format, physical groups need to be defined such that gmsh will create a `facet_tags` object upon loading.
