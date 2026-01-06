1. [Configuration](Configuration.md)
2. [Data Representation](Data_Format.md)
3. [Simulation Types](Simulation.md)
4. [Source Doc](SourceDoc.md)
   - [config.py](SourceDoc-Config.md)
   - [__coupling.py__](SourceDoc-Coupling.md)
   - [fnx.py](SourceDoc-Fnx.md)
   - [main.py](SourceDoc-Main.md)
   - [mesh_utils.py](SourceDoc-MeshUtils.md)
   - [meshes.py](SourceDoc-Meshes.md)
   - [simulation.py](SourceDoc-Simulation.md)
   - [util.py](SourceDoc-Util.md)

### Global Classes:

- [CouplingBuffer](#CouplingBuffer)
- [Projectors](#Projectors)
- [Mergers](#Mergers)
- [DataTransformer](#DataTransformers)

#### CouplingBuffer
When coupling with preCICE the size of vector data per vertex is limited to the dimension of the mesh.
Therefore, when a fem.Function should be coupled with value_size > dim(mesh) multiple transfers must be
performed via the preCICE-fenicsx adapter. CouplingBuffer provides an automatic way to link the reference
fem.Function to the required buffers for coupling.

CouplingBuffer utilizes a Projector: Its purpose is to map the original fem.Function
to a np.ndarray of shape (num_dofs x num_buffers x buffer_size).

CouplingBuffer utilizes a Merger: Its purpose is to write the data within the buffers
back into the original fem.Function

#### Projectors
Are used by CouplingBuffer to map the original fem.Function
to a np.ndarray of shape (num_dofs x num_buffers x buffer_size).

Available Implementations are:

- InplaceSplitter : only reshapes the underlying np.ndarray, thus not requiring further transformation buffers
- SelectionSplitter : if not all data needs to be transferred, a selection can be provided which will then be applied to the result of the InplaceSplitter

#### Mergers
Are used by CouplingBuffer to write the data within the buffers
back into the original fem.Function.

Available Implementations are:

- InplaceMerger : writes data from buffers directly back into original, thus not requiring further transformation buffers
- SelectionMerger : First merges data into internal buffer using InplaceMerge, then applies selection to reduce or replicate data to revert effect of SelectionSplitter

#### DataTransformers
As the data format within Meso FenicsX may not align with the format within the respective micro solver, transformations are required.
This is handled by DataTransformer. It applies the changes within the buffers of the CouplingBuffer.
Therefore, DataTransformer must be called in the end when writing to a CouplingBuffer and called first when reading from a CouplingBuffer.