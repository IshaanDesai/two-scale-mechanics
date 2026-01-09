1. [Configuration](Configuration.md)
2. [Data Representation](Data_Format.md)
3. [Simulation Types](Simulation.md)
4. [Source Doc](SourceDoc.md)
   - [config.py](SourceDoc-Config.md)
   - [coupling.py](SourceDoc-Coupling.md)
   - [__fnx.py__](SourceDoc-Fnx.md)
   - [main.py](SourceDoc-Main.md)
   - [mesh_utils.py](SourceDoc-MeshUtils.md)
   - [meshes.py](SourceDoc-Meshes.md)
   - [simulation.py](SourceDoc-Simulation.md)
   - [util.py](SourceDoc-Util.md)

### Global Classes:

- [Evaluator](#Evaluator)
- [MesoProblem](#MesoProblem)
- [MultiscaleProblem](#MultiscaleProblem)

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
