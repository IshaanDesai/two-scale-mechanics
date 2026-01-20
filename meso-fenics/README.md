# Meso-scale problem solved using FEniCSx

Solver for small or large strain solid mechanics systems using FEniCSx. Simulations can run on arbitrary meshes in .msh format. Supported simulation types are pure mesoscopic, pseudo coupled or coupled.

## Dependencies

The dependencies are

- [dolfinx](https://github.com/FEniCS/dolfinx)
- fenicsxprecice (the *develop* branch of the repository [fenicsx-adapter](https://github.com/precice/fenicsx-adapter) needs to be installed.)
- numpy  >1, <2
- gmsh
- petsc4py
- mpi4py

Apart from dolfinx, all other dependencies can be installed using the [requirements.txt](requirements.txt).

Depending on the system, these dependencies can be installed using pip, conda, or as system packages.

## How to run

In the [examples](examples/) folder, there are scripts for a [pseudo problem](examples/pseudo-coupling/), a [3D bar](examples/coupled-bar/), and a [3D notch](examples/coupled-notch/). To run these problems, for example, the 3D bar, run:

```bash
python macro.py examples/coupled-bar/config-coupled-bar.json
```

Alternatively, the [run.sh](run.sh) creates a virtual environment, installs the dependencies in it, and runs the problem. The problem to be run can be passed in the following way:

```bash
./run.sh -problem bar
```

**Note**: The [requirements.txt](requirements.txt) does not have all the dependencies because of dependency version conflicts that pip cannot resolve. It is recommended to install petsc4py, mpi4py as system packages.

For more information about implementation details or how to create a config file see the [docs](doc/Overview.md).

## Known Issues

When providing a mesh in .msh format, physical groups need to be defined such that gmsh will create a `facet_tags` object upon loading.
