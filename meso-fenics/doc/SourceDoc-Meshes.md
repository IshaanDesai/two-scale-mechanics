1. [Configuration](Configuration.md)
2. [Data Representation](Data_Format.md)
3. [Simulation Types](Simulation.md)
4. [Source Doc](SourceDoc.md)
   - [config.py](SourceDoc-Config.md)
   - [coupling.py](SourceDoc-Coupling.md)
   - [fnx.py](SourceDoc-Fnx.md)
   - [main.py](SourceDoc-Main.md)
   - [mesh_utils.py](SourceDoc-MeshUtils.md)
   - [__meshes.py__](SourceDoc-Meshes.md)
   - [simulation.py](SourceDoc-Simulation.md)
   - [util.py](SourceDoc-Util.md)

### Global Classes:

- [Mesh](#Meshes)
- [Bar](#Bar)
- [Notch](#Notch)

#### Meshes

Loads mesh and normalizes it to fit within unit cube. Constructs boundary condition objects
based on mesh and config parameters.

#### Bar

Customizes Mesh using predefined Bar mesh.

#### Notch

Customizes Mesh using predefined Notch mesh.
