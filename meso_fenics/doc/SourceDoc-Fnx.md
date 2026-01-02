---
title: Source Documentation - fnx.py
permalink: src-doc-fnx.html
keywords: source, fnx, doc
summary: Source Documentation of Meso FenicsX
---

1. [Configuration](configuration.html)
2. [Data Representation](data-format.html)
3. [Simulation Types](simulations.html)
4. [Source Doc](src-doc.html)
   - [config.py](src-doc-config.html)
   - [coupling.py](src-doc-coupling.html)
   - [__fnx.py__](src-doc-fnx.html)
   - [main.py](src-doc-main.html)
   - [mesh_utils.py](src-doc-mesh-utils.html)
   - [meshes.py](src-doc-meshes.html)
   - [simulation.py](src-doc-simulation.html)
   - [util.py](src-doc-util.html)

### Global Classes:

- [Evaluator](src-doc-fnx.html#Evaluator)
- [MesoProblem](src-doc-fnx.html#MesoProblem)
- [MultiscaleProblem](src-doc-fnx.html#MultiscaleProblem)

#### Evaluator
Can be used to evaluate an arbitrary fem expression. Call `interpolate` to evaluate based on current solution.
Data can be accessed via `var_val` field as a fem.Function.

#### MesoProblem
Constructs the pure mesoscopic system of equations using FenicsX for small and large strain.
In small strain, the governing equations are:
$$
\nabla \cdot \sigma = 0
$$
$$
\epsilon = \nabla_s u
$$
$$
\psi = \frac{1}{2}\lambda \cdot (1+\frac{\alpha}{2}tr(\epsilon)^2) \cdot tr(\epsilon)^2 + \mu (1 + \frac{\alpha}{2}<\epsilon,\epsilon>) \cdot <\epsilon,\epsilon>
$$
$$
\sigma = \frac{\partial\psi}{\partial\epsilon}
$$
$$
\tau = \frac{\partial\sigma}{\partial\epsilon}
$$

Weak Form:
$$
\int_\Omega (\nabla\cdot\sigma)\cdot v dx = 0 \Rightarrow \int_\Omega (\tau : \nabla_s u) : \nabla_s v dx = \int_{\Gamma_{N}} t \cdot v ds
$$

In large strain, the governing equations are:
$$
\nabla \cdot P = 0
$$
$$
F = Id + \nabla u
$$
$$
C = F^TF
$$
$$
E = \frac{1}{2}(C-Id)
$$
$$
\psi = \frac{1}{2}\lambda \cdot (1+\frac{\alpha}{2}tr(\epsilon)^2) \cdot tr(\epsilon)^2 + \mu (1 + \frac{\alpha}{2}<\epsilon,\epsilon>) \cdot <\epsilon,\epsilon>
$$
$$
S = \frac{\partial\psi}{\partial E}
$$
$$
P = FS
$$
P is the first Piola Kirchhoff Stress and is comparable to $\sigma$ in small strain.

#### MultiscaleProblem
The MultiscaleProblem inherits from MesoProblem, however it removes $\psi$ and replaces
in small strain $\sigma$ as well as $\tau$ and in large strain $S$ with fem.Functions.
Those need to be populated either via coupling or other means.
