# Meso-scale CalculiX simulation

One C3D8 element with a corner boundary condition.

## Run the participant

TODO

## Notation

CalculiX uses the **Voigt notation**. Tensors are ordered in the following way

Stress tensor is represented as a vector

$$
\begin{bmatrix}
\sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
\end{bmatrix}
$$

Strain tensor is represented as a vector

$$
\begin{bmatrix}
\varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & 2\varepsilon_{12} & 2\varepsilon_{13} & 2\varepsilon_{23}
\end{bmatrix}
$$

Stiffness tensor is represented as a matrix

$$
\begin{bmatrix}
C_{(11)(11)} & C_{(11)(22)} & C_{(11)(33)} & C_{(11)(12)} & C_{(11)(13)} & C_{(11)(23)}\\
C_{(22)(11)} & C_{(22)(22)} & C_{(22)(33)} & C_{(22)(12)} & C_{(22)(13)} & C_{(22)(23)}\\
C_{(33)(11)} & C_{(33)(22)} & C_{(33)(33)} & C_{(33)(12)} & C_{(33)(13)} & C_{(33)(23)}\\
C_{(12)(11)} & C_{(12)(22)} & C_{(12)(33)} & C_{(12)(12)} & C_{(12)(13)} & C_{(12)(23)}\\
C_{(13)(11)} & C_{(13)(22)} & C_{(13)(33)} & C_{(13)(12)} & C_{(13)(13)} & C_{(13)(23)}\\
C_{(23)(11)} & C_{(23)(22)} & C_{(23)(33)} & C_{(23)(12)} & C_{(23)(13)} & C_{(23)(23)}\\
\end{bmatrix}
$$
