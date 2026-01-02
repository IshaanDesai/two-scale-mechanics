---
title: Data Format
permalink: data-format.html
keywords: data format, doc
summary: Data format within Meso FenicsX
---

1. [Configuration](configuration.html)
2. [__Data Representation__](data-format.html)
3. [Simulation Types](simulations.html)
4. [Source Doc](src-doc.html)
   - [config.py](src-doc-config.html)
   - [coupling.py](src-doc-coupling.html)
   - [fnx.py](src-doc-fnx.html)
   - [main.py](src-doc-main.html)
   - [mesh_utils.py](src-doc-mesh-utils.html)
   - [meshes.py](src-doc-meshes.html)
   - [simulation.py](src-doc-simulation.html)
   - [util.py](src-doc-util.html)

In small strain setting, we use standard mandel notation for optimization purposes.
Therefore, the stress tensor is represented as a vector:

$$
\begin{bmatrix}
\sigma_{11} & \sigma_{22} & \sigma_{33} & \sqrt2\sigma_{23} & \sqrt2\sigma_{13} & \sqrt2\sigma_{12}
\end{bmatrix}
$$

The strain tensor is also represented as a vector:

$$
\begin{bmatrix}
\varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & \sqrt2\varepsilon_{23} & \sqrt2\varepsilon_{13} & \sqrt2\varepsilon_{12}
\end{bmatrix}
$$

The stiffness tensor can then be represented as a second order tensor (6x6):

$$
\begin{bmatrix}
C_{(11)(11)} & C_{(11)(22)} & C_{(11)(33)} & \sqrt2C_{(11)(23)} & \sqrt2C_{(11)(13)} & \sqrt2C_{(11)(12)}\\
C_{(22)(11)} & C_{(22)(22)} & C_{(22)(33)} & \sqrt2C_{(22)(23)} & \sqrt2C_{(22)(13)} & \sqrt2C_{(22)(12)}\\
C_{(33)(11)} & C_{(33)(22)} & C_{(33)(33)} & \sqrt2C_{(33)(23)} & \sqrt2C_{(33)(13)} & \sqrt2C_{(33)(12)}\\
\sqrt2C_{(23)(11)} & \sqrt2C_{(23)(22)} & \sqrt2C_{(23)(33)} & 2C_{(23)(23)} & 2C_{(23)(13)} & 2C_{(23)(12)}\\
\sqrt2C_{(13)(11)} & \sqrt2C_{(13)(22)} & \sqrt2C_{(13)(33)} & 2C_{(13)(23)} & 2C_{(13)(13)} & 2C_{(13)(12)}\\
\sqrt2C_{(12)(11)} & \sqrt2C_{(12)(22)} & \sqrt2C_{(12)(33)} & 2C_{(12)(23)} & 2C_{(12)(13)} & 2C_{(12)(12)}\\
\end{bmatrix}
$$

In large strain setting, no augmentation is performed.