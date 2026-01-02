---
title: Source Documentation - simulation.py
permalink: src-doc-simulation.html
keywords: source, simulation, doc
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
   - [meshes.py](src-doc-meshes.html)
   - [__simulation.py__](src-doc-simulation.html)
   - [util.py](src-doc-util.html)

### Global Classes:

- [MesoSim](src-doc-simulation.html#MesoSim)
- [PseudoCoupledSim](src-doc-simulation.html#PseudoCoupledSim)
- [CoupledSim](src-doc-simulation.html#CoupledSim)

#### MesoSim
Computes the pure mesoscopic problem.

#### PseudoCoupledSim
Computes the Multiscale Problem using input data to populate functions.

#### CoupledSim
Computes the Multiscale Problem by coupling micro simulations. 