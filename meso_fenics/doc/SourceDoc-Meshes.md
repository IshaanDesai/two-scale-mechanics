---
title: Source Documentation - meshes.py
permalink: src-doc-meshes.html
keywords: source, meshes, doc
summary: Source Documentation of Meso FenicsX
---

1. [Configuration](configuration.html)
2. [Data Representation](data-format.html)
3. [Simulation Types](simulations.html)
4. [Source Doc](src-doc.html)
   - [config.py](src-doc-config.html)
   - [coupling.py](src-doc-coupling.html)
   - [fnx.py](src-doc-fnx.html)
   - [main.py](src-doc-main.html)
   - [mesh_utils.py](src-doc-mesh-utils.html)
   - [__meshes.py__](src-doc-meshes.html)
   - [simulation.py](src-doc-simulation.html)
   - [util.py](src-doc-util.html)

### Global Classes:

- [Mesh](src-doc-meshes.html#Mesh)
- [Bar](src-doc-meshes.html#Bar)
- [Notch](src-doc-meshes.html#Notch)

#### Meshes
Loads mesh and normalizes it to fit within unit cube. Constructs boundary condition objects
based on mesh and config parameters.

#### Bar
Customizes Mesh using predefined Bar mesh.

#### Notch
Customizes Mesh using predefined Notch mesh.