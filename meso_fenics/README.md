# Meso FenicsX model

Current implementation supports a bar or notch mesh.

### How to run

Run with:
`python macro.py [--case <bar|notch>] [--micro <ADA|NASMAT|pyFANS>]`

and select the corresponding option for the used micro-solver

case: Argument | Description | Default
--- | --- | ---
bar | Bar Mesh | true
notch | Notch Mesh | false


micro: Argument | Description | Default
--- | --- | ---
ADA | pyFANS-Adapter or Surrogate micro-model | true
NASMAT | NASMAT micro-model | false
pyFANS | pyFANS micro-model | false



### Internal Data Representation

We use the the standard mandel notation with:

Stress tensor represented as a vector:

$$
\begin{bmatrix}
\sigma_{11} & \sigma_{22} & \sigma_{33} & \sqrt2\sigma_{23} & \sqrt2\sigma_{13} & \sqrt2\sigma_{12}
\end{bmatrix}
$$

Strain tensor represented as a vector

$$
\begin{bmatrix}
\varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & \sqrt2\varepsilon_{23} & \sqrt2\varepsilon_{13} & \sqrt2\varepsilon_{12}
\end{bmatrix}
$$

Stiffness tensor represented as a matrix

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
