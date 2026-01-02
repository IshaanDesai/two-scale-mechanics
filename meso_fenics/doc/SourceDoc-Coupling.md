---
title: Source Documentation - coupling.py
permalink: src-doc-coupling.html
keywords: source, coupling, doc
summary: Source Documentation of Meso FenicsX
---

1. [Configuration](configuration.html)
2. [Data Representation](data-format.html)
3. [Simulation Types](simulations.html)
4. [Source Doc](src-doc.html)
   - [config.py](src-doc-config.html)
   - [__coupling.py__](src-doc-coupling.html)
   - [fnx.py](src-doc-fnx.html)
   - [main.py](src-doc-main.html)
   - [mesh_utils.py](src-doc-mesh-utils.html)
   - [meshes.py](src-doc-meshes.html)
   - [simulation.py](src-doc-simulation.html)
   - [util.py](src-doc-util.html)

### Global Classes:

- [CouplingBuffer](src-doc-coupling.html#CouplingBuffer)
- [Projectors](src-doc-coupling.html#Projectors)
- [Mergers](src-doc-coupling.html#Mergers)
- [DataTransformer](src-doc-coupling.html#DataTransformers)

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