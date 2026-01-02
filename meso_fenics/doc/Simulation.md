---
title: Simulation Types
permalink: simulations.html
keywords: configuration, doc
summary: Simulation types within Meso FenicsX
---

1. [Configuration](configuration.html)
2. [Data Representation](data-format.html)
3. [__Simulation Types__](simulations.html)
4. [Source Doc](src-doc.html)
   - [config.py](src-doc-config.html)
   - [coupling.py](src-doc-coupling.html)
   - [fnx.py](src-doc-fnx.html)
   - [main.py](src-doc-main.html)
   - [mesh_utils.py](src-doc-mesh-utils.html)
   - [meshes.py](src-doc-meshes.html)
   - [simulation.py](src-doc-simulation.html)
   - [util.py](src-doc-util.html)

Meso FenicsX is capable of running in different simulation modes to fulfill different purposes.

1. MesoSim: only runs the pure mesoscopic simulation, for example to verify the mesh or boundary conditions
2. PseudoCoupledSim: Mocks the CoupledSim by loading in the final stresses and optionally tangent data
3. CoupledSim: Couples the mesoscopic simulation to micro simulations at each DOF using preCICE and micro-manager