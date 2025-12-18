# Fourier Accelerated Nodal Solvers (FANS)

## Install FANS

Clone [FANS](https://github.com/DataAnalyticsEngineering/FANS) and build [pyFANS](https://github.com/DataAnalyticsEngineering/FANS/tree/develop/pyfans).

## Configure the FANS simulation for the Micro Manager

FANS is configured in a JSON configuration file. The name `input.json` is fixed. More information about this file is in the [FANS configuration documentation](https://github.com/DataAnalyticsEngineering/FANS?tab=readme-ov-file#input-file-format).

## Run the FANS simulation

Point to the `pyFANS.so` file in the `micro_file_name` entry in the `micro-manager-config.json` configuration.

Run the Micro Manager

```bash
micro-manager-precice micro-manager-config.json
```

## Notation

FANS uses the **Mandel notation**. Tensors are ordered in the following way

Stress tensor represented as a vector

$$
\begin{bmatrix}
\sigma_{11} & \sigma_{22} & \sigma_{33} & \sqrt2\sigma_{12} & \sqrt2\sigma_{13} & \sqrt2\sigma_{23}
\end{bmatrix}
$$

Strain tensor represented as a vector

$$
\begin{bmatrix}
\varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & \sqrt2\varepsilon_{12} & \sqrt2\varepsilon_{13} & \sqrt2\varepsilon_{23}
\end{bmatrix}
$$

Stiffness tensor represented as a matrix

$$
\begin{bmatrix}
C_{(11)(11)} & C_{(11)(22)} & C_{(11)(33)} & \sqrt2C_{(11)(12)} & \sqrt2C_{(11)(13)} & \sqrt2C_{(11)(23)}\\
C_{(22)(11)} & C_{(22)(22)} & C_{(22)(33)} & \sqrt2C_{(22)(12)} & \sqrt2C_{(22)(13)} & \sqrt2C_{(22)(23)}\\
C_{(33)(11)} & C_{(33)(22)} & C_{(33)(33)} & \sqrt2C_{(33)(12)} & \sqrt2C_{(33)(13)} & \sqrt2C_{(33)(23)}\\
\sqrt2C_{(12)(11)} & \sqrt2C_{(12)(22)} & \sqrt2C_{(12)(33)} & 2C_{(12)(12)} & 2C_{(12)(13)} & 2C_{(12)(23)}\\
\sqrt2C_{(13)(11)} & \sqrt2C_{(13)(22)} & \sqrt2C_{(13)(33)} & 2C_{(13)(12)} & 2C_{(13)(13)} & 2C_{(13)(23)}\\
\sqrt2C_{(23)(11)} & \sqrt2C_{(23)(22)} & \sqrt2C_{(23)(33)} & 2C_{(23)(12)} & 2C_{(23)(13)} & 2C_{(23)(23)}\\
\end{bmatrix}
$$
