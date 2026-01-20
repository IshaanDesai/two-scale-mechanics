1. [Configuration](Configuration.md)
2. [Data Representation](Data_Format.md)
3. [__Simulation Types__](Simulation.md)
4. [Source Doc](SourceDoc.md)
   - [config.py](SourceDoc-Config.md)
   - [coupling.py](SourceDoc-Coupling.md)
   - [fnx.py](SourceDoc-Fnx.md)
   - [main.py](SourceDoc-Main.md)
   - [mesh_utils.py](SourceDoc-MeshUtils.md)
   - [meshes.py](SourceDoc-Meshes.md)
   - [simulation.py](SourceDoc-Simulation.md)
   - [util.py](SourceDoc-Util.md)

Meso FenicsX is capable of running in different simulation modes to fulfill different purposes.

1. MesoSim: only runs the pure mesoscopic simulation, for example to verify the mesh or boundary conditions
2. PseudoCoupledSim: Mocks the CoupledSim by loading in the final stresses and optionally tangent data
3. CoupledSim: Couples the mesoscopic simulation to micro simulations at each DOF using preCICE and micro-manager