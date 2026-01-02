from typing import Optional

import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io
from dolfinx.fem.petsc import NonlinearProblem, LinearProblem
import gmsh

from .config import Config
from .meshes import Mesh

class Evaluator:
    """
    Evaluates the given fenicsx variable/expression for a function space.
    """
    def __init__(self, var, space):
        """
        Constructs evaluator for

        Params
        ------
        var: object
            fenicsx variable/expression
        space: fem.Functionspace
            fenicsx space
        """
        self.var_exp = fem.Expression(var, space.element.interpolation_points)
        self.var_val = fem.Function(space)

    def interpolate(self):
        """
        Interpolates the expression and populates self.var_val.
        """
        self.var_val.interpolate(self.var_exp)
        return self

class MesoProblem:
    def __init__(self, config: Config, mesh: Mesh):
        self.mesh = mesh
        self.petsc_options = {
            "snes_type": "newtonls",
            "snes_linesearch_type": "none",
            "snes_atol": 1e-8,
            "snes_rtol": 1e-8,
            "snes_monitor": None,
            "ksp_error_if_not_converged": True,
            "ksp_type": "gmres",
            "ksp_rtol": 1e-10,
            "ksp_max_it": 1000,
            "ksp_monitor": None,
            "pc_type": "hypre",
            "pc_hypre_type": "boomeramg",
            "pc_hypre_boomeramg_max_iter": 1,
            "pc_hypre_boomeramg_cycle_type": "v",
        }

        # material properties
        self.lam   = fem.Constant(self.mesh.domain, config.problem_lambda)
        self.mu    = fem.Constant(self.mesh.domain, config.problem_mu    )
        self.alpha = fem.Constant(self.mesh.domain, config.problem_alpha )

        # FEM
        element_degree = config.problem_element_degree
        self.V = fem.functionspace(self.mesh.domain, ('CG', element_degree, (3,)))
        self.W3 = fem.functionspace(self.mesh.domain, ("CG", element_degree, (3,)))
        self.uh = fem.Function(self.V)
        self.uh.name = "displacement"
        self.u = ufl.TrialFunction(self.V)
        self.v = ufl.TestFunction(self.V)
        self.vm_stress_fun = fem.Function(self.V)
        self.vm_stress_fun.name = 'von_mises_stress'

        self.bc_dc = self.mesh.get_bc_diriclet(self.V)
        self.bc_nm = self.mesh.get_bc_neumann(self.v)

        self.is_small_strain = None
        if config.problem_strain_type == 'small_strain':
            self.is_small_strain = True
            self.W  = fem.functionspace(self.mesh.domain, ("CG", element_degree, (6,)))
            self.WT = fem.functionspace(self.mesh.domain, ("CG", element_degree, (6, 6)))

            self.eps_var = ufl.variable(self.symgrad_mandel(self.uh))
            self.sig_fun = fem.Function(self.W)
            self.tan_fun = fem.Function(self.WT)

            self.psi_exp = self.calc_psi()
            self.sigma_exp = ufl.diff(self.psi_exp, self.eps_var)
            self.tangent_exp = ufl.diff(self.sigma_exp, self.eps_var)

            self.meso_res = ufl.inner(self.sigma_exp, self.symgrad_mandel(self.v)) * ufl.dx - self.bc_nm(self.v)
            self.meso_jac = ufl.inner(ufl.dot(self.tangent_exp, self.symgrad_mandel(self.u)), self.symgrad_mandel(self.v)) * ufl.dx
        elif config.problem_strain_type == 'large_strain':
            self.is_small_strain = False
            self.W  = fem.functionspace(self.mesh.domain, ("CG", element_degree, (3, 3)))
            self.WT = None          # only use tangent implicitly

            Id = ufl.Identity(3)
            self.F = Id + ufl.grad(self.uh)
            self.eps_exp = fem.Constant(self.mesh.domain, 0.5) * (self.F.T * self.F - Id)
            self.eps_var = ufl.variable(self.eps_exp)
            self.sig_fun = fem.Function(self.W)
            self.tan_fun = None     # only use tangent implicitly

            self.psi_exp = self.calc_psi()
            self.sigma_exp = ufl.diff(self.psi_exp, self.eps_var)
            self.tangent_exp = None

            self.meso_res = ufl.inner(self.F * self.sigma_exp, ufl.grad(self.v)) * ufl.dx - self.bc_nm(self.v)
            self.meso_jac = ufl.derivative(self.meso_res, self.uh, self.u)
        else:
            raise ValueError("Unknown strain type")

        self.meso_problem = NonlinearProblem(
            self.meso_res,
            self.uh,
            bcs=self.bc_dc,
            J=self.meso_jac,
            petsc_options=self.petsc_options,
            petsc_options_prefix="nonlinpoisson"
        )

    def solve(self):
        self.meso_problem.solve()

        sig_eval = Evaluator(ufl.variable(self.sigma_exp), self.W).interpolate()
        self.sig_fun.x.array[:] = sig_eval.var_val.x.array[:]

        if self.is_small_strain:
            tan_eval = Evaluator(ufl.variable(self.tangent_exp), self.WT).interpolate()
            self.tan_fun.x.array[:] = tan_eval.var_val.x.array[:]

    def calc_psi(self):
        tr_e = None
        if self.is_small_strain: tr_e = self.eps_var[0] + self.eps_var[1] + self.eps_var[2]
        else: tr_e = ufl.tr(self.eps_var)
        e2 = ufl.inner(self.eps_var, self.eps_var)
        tr_e2 = tr_e ** 2
        half_alpha = 0.5 * self.alpha
        return 0.5 * self.lam * (1.0 + half_alpha * tr_e2) * tr_e2 + self.mu * (1 + half_alpha * e2) * e2

    def calc_von_mises_stress(self):
        if self.is_small_strain:
            sig = self.sig_fun.x.array.reshape(-1, 6)
            s1 = 0.5 * (np.power(sig[:, 0] - sig[:, 1], 2) +
                        np.power(sig[:, 1] - sig[:, 2], 2) +
                        np.power(sig[:, 2] - sig[:, 0], 2))
            # we div by 2 to remove sqrt 2 from sig
            # accounted by factor 3 -> 1.5
            s2 = 1.5 * (np.power(sig[:, 3], 2) + np.power(sig[:, 4], 2) + np.power(sig[:, 5], 2))
            self.vm_stress_fun.x.array.reshape(-1, 3)[:, 0] = np.sqrt(s1 + s2)
        else:
            sig = self.sig_fun.x.array.reshape(-1, 3, 3)
            s1 = 0.5 * (np.power(sig[:, 0, 0] - sig[:, 1, 1], 2) +
                        np.power(sig[:, 2, 1] - sig[:, 2, 2], 2) +
                        np.power(sig[:, 2, 2] - sig[:, 0, 0], 2))
            s2 = 3.0 * (np.power(sig[:, 0, 1], 2) +
                        np.power(sig[:, 1, 2], 2) +
                        np.power(sig[:, 2, 0], 2))
            self.vm_stress_fun.x.array.reshape(-1, 3)[:, 0] = np.sqrt(s1 + s2)

    def symgrad_mandel(self, vec):
        halfsqrt2 = 0.5 * np.sqrt(2)
        return ufl.as_vector([vec[0].dx(0), vec[1].dx(1), vec[2].dx(2),
                              halfsqrt2 * (vec[1].dx(2) + vec[2].dx(1)),
                              halfsqrt2 * (vec[0].dx(2) + vec[2].dx(0)),
                              halfsqrt2 * (vec[0].dx(1) + vec[1].dx(0))])

class MultiscaleProblem(MesoProblem):
    def __init__(self, config: Config, mesh: Mesh):
        super().__init__(config, mesh)

        if self.is_small_strain:
            a = ufl.inner(ufl.dot(self.tan_fun, self.symgrad_mandel(self.u)), self.symgrad_mandel(self.v)) * ufl.dx
            L = self.bc_nm(self.v)
            self.ms_problem = LinearProblem(
                a=a,
                L=L,
                bcs=self.bc_dc,
                petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
                petsc_options_prefix="Poisson",
            )
        else:
            self.res = ufl.inner(self.F * self.sig_fun, ufl.grad(self.v)) * ufl.dx - self.bc_nm(self.v)
            self.jac = ufl.derivative(self.res, self.uh, self.u)
            self.ms_problem = NonlinearProblem(
                self.meso_res,
                self.uh,
                bcs=self.bc_dc,
                J=self.meso_jac,
                petsc_options=self.petsc_options,
                petsc_options_prefix="nonlinpoisson"
            )

    def solve_meso(self):
        super().solve()

    def solve(self):
        if self.is_small_strain:
            result = self.ms_problem.solve()
            self.uh.x.array[:] = result.x.array[:]
        else:
            self.ms_problem.solve()

