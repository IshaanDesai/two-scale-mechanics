---
title: Source Documentation - config.py
permalink: src-doc-config.html
keywords: source, config, doc
summary: Source Documentation of Meso FenicsX
---

1. [Configuration](configuration.html)
2. [Data Representation](data-format.html)
3. [Simulation Types](simulations.html)
4. [Source Doc](src-doc.html)
   - [__config.py__](src-doc-config.html)
   - [coupling.py](src-doc-coupling.html)
   - [fnx.py](src-doc-fnx.html)
   - [main.py](src-doc-main.html)
   - [mesh_utils.py](src-doc-mesh-utils.html)
   - [meshes.py](src-doc-meshes.html)
   - [simulation.py](src-doc-simulation.html)
   - [util.py](src-doc-util.html)

Simple JSON based config format. 
To extend, add backing field in `__init__`, populate in `load` and mark accessor as property.