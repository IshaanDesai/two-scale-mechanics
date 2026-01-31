from functools import partial
from typing import Optional

import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io
from dolfinx.fem.petsc import (NonlinearProblem, LinearProblem, assign,
                               assemble_vector, assemble_matrix, create_vector,
                               create_matrix, apply_lifting, set_bc)
import gmsh

from .config import Config
from .meshes import Mesh


class Anderson_AdaGuard:
    def __init__(self, r=1e-4):
        self.w_hist = list()
        self.u_hist = list()
        self.r = r

    def insert(self, w: np.ndarray, u: np.ndarray):
        self.w_hist.append(w)
        self.u_hist.append(u)
        if len(self.w_hist) > 2:
            self.w_hist.pop(0)
            self.u_hist.pop(0)

    def gamma_safeguard(self, gamma: float):
        nu   = np.linalg.norm(self.w_hist[1]) / np.linalg.norm(self.w_hist[0])
        r    = np.minimum(nu, self.r)
        beta = r * nu

        if gamma == 0 or gamma >= 1:
            return 0
        elif np.abs(gamma) / np.abs(1 - gamma) > beta:
            return beta / (gamma * (beta * np.sign(gamma)))
        else:
            return 0

    def accelerate(self, du: fem.Function, uh_old: fem.Function):
        self.insert(du.x.array.copy(), uh_old.x.array.copy())
        if len(self.w_hist) == 1:
            return self.w_hist[0]

        delta_w = self.w_hist[1] - self.w_hist[0]
        gamma = np.inner(delta_w, self.w_hist[1]) / np.sum(delta_w * delta_w)
        lam = self.gamma_safeguard(gamma)

        delta_u = self.u_hist[1] - self.u_hist[0]
        return self.w_hist[1] - lam * gamma * (delta_u + delta_w)


class LineSearchStep:
    def __init__(self, uh: fem.Function, calc_res: callable, max_iter=10, alpha_start=1.0, tau=0.8, alpha_min=1e-4, r=1e-4):
        self.uh = uh
        self.calc_res = calc_res
        self.max_iter = max_iter
        self.alpha_start = alpha_start
        self.alpha_min = alpha_min
        self.tau = tau
        self.reduction_factor = 0.5

        self.is_active = False
        self.curr_alpha = alpha_start
        self.uh0 = np.zeros_like(uh.x.array)
        self.res0 = 0
        self.last_res = self.res0
        self.curr_du = None
        self.iter = 0
        self.suf_decrease = 1e-3

        self.aa = Anderson_AdaGuard(r)

    def is_searching(self):
        return self.is_active

    def initialize(self, du: fem.Function):
        # call after ksp solve
        self.is_active = True
        self.iter = 0
        self.curr_alpha = min(self.alpha_start, 1.2*self.curr_alpha)
        self.uh0[:] = self.uh.x.array
        self.res0 = self.calc_res().norm(0)
        self.last_res = self.res0

        du_aa = self.aa.accelerate(du, self.uh)
        self.curr_du = du_aa
        self.update_uh(self.curr_alpha, du_aa)

    def step(self):
        new_res = self.calc_res().norm(0)
        ignore_res_check = False
        if self.iter == 0: # just computed R(u + 1 * du) = F
            if new_res <= (1 - self.suf_decrease * self.curr_alpha) * self.res0:
                # no need to perform LS
                self.is_active = False
                return True
            else:
                print(f"Solving LS - iter: {self.iter}, alpha: 0, res: {self.res0}")
                ignore_res_check = True

        if not ignore_res_check and (new_res > self.last_res):
            self.is_active = False
            print(f"Solving LS - iter: {self.iter}, bad alpha: {self.curr_alpha}, res: {new_res}")
            self.curr_alpha /= self.tau
            return True

        if new_res < self.res0 * self.reduction_factor:
            self.is_active = False
            print(f"Final LS - iter: {self.iter}, accepted alpha: {self.curr_alpha}, res: {new_res}")
            return True # last attempt was accepted

        self.last_res = new_res
        self.curr_alpha *= self.tau
        self.iter += 1

        if self.curr_alpha < self.alpha_min or self.iter >= self.max_iter:
            self.is_active = False
            return True

        print(f"Solving LS - iter: {self.iter}, alpha: {self.curr_alpha}, res: {new_res}")
        self.update_uh(self.curr_alpha, self.curr_du)
        return False

    def update_uh(self, alpha, du):
        self.uh.x.array[:] = self.uh0 + alpha * du
        self.uh.x.scatter_forward()


class GradApprox:
    def __init__(self, task, steps, uh):
        self.task = task
        self.max_steps = steps
        self.steps = 0
        self.uh0 = np.zeros_like(uh.x.array)
        self.uh = uh
        self.du = None

    def initialize(self, du):
        self.steps = 0
        self.du = du
        self.uh0[:] = self.uh.x.array[:]

    def solve(self):
        self.task(self.du, self.uh0, self.uh, self.steps)
        self.steps += 1
        return self.is_running()

    def is_running(self):
        return self.steps < self.max_steps

    def finalize(self):
        self.task(self.du, self.uh0, self.uh, self.max_steps)
        self.steps = self.max_steps
        self.uh.x.array[:] = self.uh0
        self.uh.x.scatter_forward()



class NonlinearProblemStep:
    def __init__(self, comm, uh: fem.Function, F, J, bcs: list, options: dict = None, grad_approx: GradApprox = None):
        self.res = fem.form(F)
        self.jac = fem.form(J)
        self.du = fem.Function(uh.function_space)
        self.uh = uh
        self.A = create_matrix(self.jac)
        self.L = create_vector(fem.extract_function_spaces(self.res))

        self.bcs = bcs
        self.iter = 0
        self.max_iter = 100
        self.threshold = 1e-10
        self.solver = PETSc.KSP().create(comm)
        self.solver.setOperators(self.A)
        self.grad_approx = grad_approx
        opts = PETSc.Options()
        if options is not None:
            if "ksp_max_it" in options:
                self.max_iter = options["ksp_max_it"]
            if "ksp_rtol" in options:
                self.threshold = options["ksp_rtol"]
            if "ksp_type" in options:
                opts["ksp_type"] = options["ksp_type"]
            if "pc_type" in options:
                opts["pc_type"] = options["pc_type"]
            if "pc_hypre_type" in options:
                opts["pc_hypre_type"] = options["pc_hypre_type"]
            if "pc_hypre_boomeramg_max_iter" in options:
                opts["pc_hypre_boomeramg_max_iter"] = options["pc_hypre_boomeramg_max_iter"]
            if "pc_hypre_boomeramg_cycle_type" in options:
                opts["pc_hypre_boomeramg_cycle_type"] = options["pc_hypre_boomeramg_cycle_type"]
        #opts.prefixPush(self.solver.getOptionsPrefix())
        self.solver.setFromOptions()
        #opts.prefixPop()
        self.pc = self.solver.getPC()

        self.last_correction = self.threshold * 10
        self.construct_L_res = partial(
            NonlinearProblemStep._compute_residual_vec,
            self.L,
            self.res,
            self.jac,
            self.bcs,
            self.uh,
        )

        self.ls = LineSearchStep(self.uh, self.construct_L_res)

    def solve(self):
        if self.iter >= self.max_iter: return True
        if self.last_correction < self.threshold: return True

        if self.grad_approx is not None and self.grad_approx.is_running():
            ongoing = self.grad_approx.solve()
            if ongoing: return False
            else: self.grad_approx.finalize()

        if self.ls.is_searching():
            success = self.ls.step()
            if not success: return False
            self.last_correction = self.du.x.petsc_vec.norm(0) * self.ls.curr_alpha
            print(f"Solving KSP - iter {self.iter} du-norm: {self.last_correction} res-norm: {self.L.norm(0)}")
            if self.grad_approx is not None:
                self.grad_approx.initialize(self.du)
        
        self._solve_ksp()
        self.ls.initialize(self.du)
        return False


    def _solve_ksp(self):
        self.construct_L_res()
        self.A.zeroEntries()
        assemble_matrix(self.A, self.jac, bcs=self.bcs)
        self.A.assemble()

        self.solver.solve(self.L, self.du.x.petsc_vec)
        self.du.x.scatter_forward()

        self.iter += 1

    @staticmethod
    def _compute_residual_vec(L, res, jac, bcs, uh):
        with L.localForm() as loc_L:
            loc_L.set(0)
        assemble_vector(L, res)
        L.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        L.scale(-1)

        apply_lifting(L, [jac], [bcs], x0=[uh.x.petsc_vec], alpha=1)
        set_bc(L, bcs, uh.x.petsc_vec, 1.0)
        L.ghostUpdate(addv=PETSc.InsertMode.INSERT_VALUES, mode=PETSc.ScatterMode.FORWARD)

        return L


class Evaluator:
    """
    Evaluates the given fenicsx variable/expression for a function space.

    Attributes
    ----------
    var_exp : fem.Expression
        The FEniCSx expression to evaluate.
    var_val : fem.Function
        The function containing interpolated values.
    """

    def __init__(self, var, space):
        """
        Construct evaluator for a variable/expression.

        Parameters
        ----------
        var : object
            FEniCSx variable/expression to evaluate.
        space : fem.FunctionSpace
            FEniCSx function space for interpolation.
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
    """
    Meso-scale nonlinear elasticity problem solver.

    This class sets up and solves a nonlinear elasticity problem using FEniCSx,
    supporting both small and large strain formulations.

    Parameters
    ----------
    config : Config
        Configuration object containing problem parameters.
    mesh : Mesh
        Mesh object containing the domain and boundary conditions.

    Attributes
    ----------
    mesh : Mesh
        The computational mesh.
    petsc_options : dict
        PETSc solver options.
    lam : fem.Constant
        Lamé parameter lambda.
    mu : fem.Constant
        Lamé parameter mu.
    alpha : fem.Constant
        Material nonlinearity parameter.
    V : fem.FunctionSpace
        Vector function space for displacement.
    uh : fem.Function
        Displacement solution.
    is_small_strain : bool
        True if using small strain formulation.
    """

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
            "ksp_rtol": 1e-8,
            "ksp_max_it": 1000,
            "ksp_monitor": None,
            "pc_type": "hypre",
            "pc_hypre_type": "boomeramg",
            "pc_hypre_boomeramg_max_iter": 1,
            "pc_hypre_boomeramg_cycle_type": "v",
        }

        # material properties
        self.lam = fem.Constant(self.mesh.domain, config.problem_lambda)
        self.mu = fem.Constant(self.mesh.domain, config.problem_mu)
        self.alpha = fem.Constant(self.mesh.domain, config.problem_alpha)

        # FEM
        element_degree = config.problem_element_degree
        self.V = fem.functionspace(self.mesh.domain, ("CG", element_degree, (3,)))
        self.W3 = fem.functionspace(self.mesh.domain, ("DG", 0, (3,)))
        self.uh = fem.Function(self.V)
        self.uh.name = "displacement"
        self.u = ufl.TrialFunction(self.V)
        self.v = ufl.TestFunction(self.V)
        self.vm_stress_fun = fem.Function(self.W3)
        self.vm_stress_fun.name = "von_mises_stress"

        self.bc_dc = self.mesh.get_bc_diriclet(self.V)
        self.bc_nm = self.mesh.get_bc_neumann(self.v)

        self.is_small_strain = None
        if config.problem_strain_type == "small_strain":
            self.is_small_strain = True
            self.W = fem.functionspace(self.mesh.domain, ("DG", 0, (6,)))
            self.WT = fem.functionspace(
                self.mesh.domain, ("DG", 0, (6, 6))
            )

            self.eps_var = ufl.variable(self.symgrad_mandel(self.uh))
            self.sig_fun = fem.Function(self.W)
            self.tan_fun = fem.Function(self.WT)

            self.psi_exp = self.calc_psi()
            self.sigma_exp = ufl.diff(self.psi_exp, self.eps_var)
            self.tangent_exp = ufl.diff(self.sigma_exp, self.eps_var)

            self.meso_res = ufl.inner(
                self.sigma_exp, self.symgrad_mandel(self.v)
            ) * ufl.dx - self.bc_nm(self.v)
            self.meso_jac = (
                ufl.inner(
                    ufl.dot(self.tangent_exp, self.symgrad_mandel(self.u)),
                    self.symgrad_mandel(self.v),
                )
                * ufl.dx
            )
        elif config.problem_strain_type == "large_strain":
            self.is_small_strain = False
            self.W = fem.functionspace(self.mesh.domain, ("DG", 0, (3, 3)))
            self.WT = None  # only use tangent implicitly

            Id = ufl.Identity(3)
            self.F = Id + ufl.grad(self.uh)
            self.eps_exp = fem.Constant(self.mesh.domain, 0.5) * (
                self.F.T * self.F - Id
            )
            self.eps_var = ufl.variable(self.eps_exp)
            self.sig_fun = fem.Function(self.W)
            self.tan_fun = None  # only use tangent implicitly

            self.psi_exp = self.calc_psi()
            self.sigma_exp = ufl.diff(self.psi_exp, self.eps_var)
            self.tangent_exp = None

            self.meso_res = ufl.inner(
                self.F * self.sigma_exp, ufl.grad(self.v)
            ) * ufl.dx - self.bc_nm(self.v)
            self.meso_jac = ufl.derivative(self.meso_res, self.uh, self.u)
        else:
            raise ValueError("Unknown strain type")

        self.meso_problem = NonlinearProblemStep(self.mesh.domain.comm, self.uh, self.meso_res, self.meso_jac, self.bc_dc, self.petsc_options)

    def solve(self):
        """
        Solve the meso-scale problem and update stress and tangent fields.

        This method solves the nonlinear problem and interpolates the
        stress tensor and (for small strain) tangent modulus.
        """
        while not self.meso_problem.solve(): pass

        sig_eval = Evaluator(ufl.variable(self.sigma_exp), self.W).interpolate()
        self.sig_fun.x.array[:] = sig_eval.var_val.x.array[:]
        self.sig_fun.x.scatter_forward()

        if self.is_small_strain:
            tan_eval = Evaluator(ufl.variable(self.tangent_exp), self.WT).interpolate()
            self.tan_fun.x.array[:] = tan_eval.var_val.x.array[:]
            self.tan_fun.x.scatter_forward()

    def calc_psi(self):
        """
        Calculate the strain energy density function.

        Returns
        -------
        ufl.core.expr.Expr
            The strain energy density expression.
        """
        tr_e = None
        if self.is_small_strain:
            tr_e = self.eps_var[0] + self.eps_var[1] + self.eps_var[2]
        else:
            tr_e = ufl.tr(self.eps_var)
        e2 = ufl.inner(self.eps_var, self.eps_var)
        tr_e2 = tr_e**2
        half_alpha = 0.5 * self.alpha
        return (
            0.5 * self.lam * (1.0 + half_alpha * tr_e2) * tr_e2
            + self.mu * (1 + half_alpha * e2) * e2
        )

    def calc_von_mises_stress(self):
        """
        Calculate the von Mises stress from the stress tensor.

        Updates the vm_stress_fun attribute with the computed von Mises stress.
        Handles both small strain (Mandel notation) and large strain formulations.
        """
        if self.is_small_strain:
            sig = self.sig_fun.x.array.reshape(-1, 6)
            s1 = 0.5 * (
                np.power(sig[:, 0] - sig[:, 1], 2)
                + np.power(sig[:, 1] - sig[:, 2], 2)
                + np.power(sig[:, 2] - sig[:, 0], 2)
            )
            # we div by 2 to remove sqrt 2 from sig
            # accounted by factor 3 -> 1.5
            s2 = 1.5 * (
                np.power(sig[:, 3], 2) + np.power(sig[:, 4], 2) + np.power(sig[:, 5], 2)
            )
            self.vm_stress_fun.x.array.reshape(-1, 3)[:, 0] = np.sqrt(s1 + s2)
        else:
            sig = self.sig_fun.x.array.reshape(-1, 3, 3)
            s1 = 0.5 * (
                np.power(sig[:, 0, 0] - sig[:, 1, 1], 2)
                + np.power(sig[:, 2, 1] - sig[:, 2, 2], 2)
                + np.power(sig[:, 2, 2] - sig[:, 0, 0], 2)
            )
            s2 = 3.0 * (
                np.power(sig[:, 0, 1], 2)
                + np.power(sig[:, 1, 2], 2)
                + np.power(sig[:, 2, 0], 2)
            )
            self.vm_stress_fun.x.array.reshape(-1, 3)[:, 0] = np.sqrt(s1 + s2)
        self.vm_stress_fun.x.scatter_forward()

    def symgrad_mandel(self, vec):
        """
        Compute the symmetric gradient in Mandel notation.

        Parameters
        ----------
        vec : ufl.core.expr.Expr
            Vector function for which to compute the symmetric gradient.

        Returns
        -------
        ufl.core.expr.Expr
            Symmetric gradient as a 6-component vector in Mandel notation.
        """
        halfsqrt2 = 0.5 * np.sqrt(2)
        return ufl.as_vector(
            [
                vec[0].dx(0),
                vec[1].dx(1),
                vec[2].dx(2),
                halfsqrt2 * (vec[1].dx(2) + vec[2].dx(1)),
                halfsqrt2 * (vec[0].dx(2) + vec[2].dx(0)),
                halfsqrt2 * (vec[0].dx(1) + vec[1].dx(0)),
            ]
        )


class MultiscaleProblem(MesoProblem):
    """
    Multiscale problem solver using pre-computed constitutive response.

    This class extends MesoProblem to solve problems where the stress and
    tangent modulus are provided from external micro-scale simulations.

    Parameters
    ----------
    config : Config
        Configuration object containing problem parameters.
    mesh : Mesh
        Mesh object containing the domain and boundary conditions.
    """

    class StressGradApproximator:
        def __init__(self, sig_fun: fem.Function, sig_grad_fun: fem.Function):
            self.sig_fun = sig_fun
            self.sig_grad_fun = sig_grad_fun
            self.sig_step_data = np.zeros((2, sig_fun.x.array.shape[0]))
            self.eps = 1e-6

        def __call__(self, du: fem.Function, uh0: np.ndarray, uh: fem.Function, step: int):
            # gather prev compute sig data
            if step > 0: self.sig_step_data[step-1][:] = self.sig_fun.x.array[:]
            # create new uh data
            if step == 0: uh.x.array[:] = uh0 + (self.eps * du.x.array / du.x.petsc_vec.norm(0))
            if step == 1: uh.x.array[:] = uh0 - (self.eps * du.x.array / du.x.petsc_vec.norm(0))
            # compute grad
            if step == 2: self.sig_grad_fun.x.array[:] = (self.sig_step_data[0] - self.sig_step_data[1]) / (2.0 * self.eps)


    def __init__(self, config: Config, mesh: Mesh):
        super().__init__(config, mesh)

        if self.is_small_strain:
            self.res = ufl.inner(self.sig_fun, self.symgrad_mandel(self.v)) * ufl.dx - self.bc_nm(self.v)
            self.jac = (
                ufl.inner(
                    ufl.dot(self.tan_fun, self.symgrad_mandel(self.u)),
                    self.symgrad_mandel(self.v),
                )
                * ufl.dx
            )
            self.ms_problem = NonlinearProblemStep(
                self.mesh.domain.comm,
                self.uh,
                self.res,
                self.jac,
                self.bc_dc,
                self.petsc_options
            )
        else:
            self.res = ufl.inner(
                self.F * self.sig_fun, ufl.grad(self.v)
            ) * ufl.dx - self.bc_nm(self.v)

            self.sig_grad_fun = fem.Function(self.W)
            self.sig_grad_approximator = MultiscaleProblem.StressGradApproximator(self.sig_fun, self.sig_grad_fun)
            self.grad_approx = GradApprox(self.sig_grad_approximator, 2, self.uh)
            self.jac = ufl.inner(
                ufl.derivative(self.F, self.uh, self.u) * self.sig_fun +
                self.F * self.sig_grad_fun * ufl.grad(self.u),
                ufl.grad(self.v)
            ) * ufl.dx
            self.ms_problem = NonlinearProblemStep(
                self.mesh.domain.comm,
                self.uh,
                self.res,
                self.jac,
                self.bc_dc,
                self.petsc_options,
                self.grad_approx
            )

    def solve_meso(self):
        """
        Solve using the meso-scale constitutive law.

        Calls the parent class solve method to compute stress and tangent
        from the built-in material model.
        """
        super().solve()

    def init_with_ms(self):
        """
        Solve problem approximately only using micro-scale tangent.
        Use as first guess for subsequent iterations.
        """
        sig = ufl.dot(self.tan_fun, self.symgrad_mandel(self.uh))
        res = ufl.inner(sig, self.symgrad_mandel(self.v)) * ufl.dx - self.bc_nm(self.v)
        jac = ufl.inner(ufl.dot(self.tan_fun, self.symgrad_mandel(self.u)), self.symgrad_mandel(self.v)) * ufl.dx
        problem = NonlinearProblem(
            res,
            self.uh,
            bcs=self.bc_dc,
            J=jac,
            petsc_options=self.petsc_options,
            petsc_options_prefix="nonlinpoisson",
        )
        problem.solve()

        sig_eval = Evaluator(ufl.variable(sig), self.W).interpolate()
        self.sig_fun.x.array[:] = sig_eval.var_val.x.array[:]
        self.sig_fun.x.scatter_forward()

    def solve(self):
        self.ms_problem.solve()
