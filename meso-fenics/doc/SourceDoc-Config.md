1. [Configuration](Configuration.md)
2. [Data Representation](Data_Format.md)
3. [Simulation Types](Simulation.md)
4. [Source Doc](SourceDoc.md)
   - [__config.py__](SourceDoc-Config.md)
   - [coupling.py](SourceDoc-Coupling.md)
   - [fnx.py](SourceDoc-Fnx.md)
   - [main.py](SourceDoc-Main.md)
   - [mesh_utils.py](SourceDoc-MeshUtils.md)
   - [meshes.py](SourceDoc-Meshes.md)
   - [simulation.py](SourceDoc-Simulation.md)
   - [util.py](SourceDoc-Util.md)

Simple JSON based config format. 
To extend, add backing field in `__init__`, populate in `load` and mark accessor as property.