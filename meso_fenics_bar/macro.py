from typing import Callable
import argparse

import dolfinx.io.gmsh
import numpy as np
import ufl

from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io
from dolfinx.fem.petsc import NonlinearProblem, LinearProblem
import gmsh

from fenicsxprecice import Adapter, CouplingMesh

class Evaluator:
    def __init__(self, var, space):
        self.var_exp = fem.Expression(var, space.element.interpolation_points)
        self.var_val = fem.Function(space)

    def interpolate(self):
        self.var_val.interpolate(self.var_exp)

class MeshParams:
    def __init__(self, x0, y0, z0, lx, ly, lz, nx=0, ny=0, nz=0):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.lx = lx
        self.ly = ly
        self.lz = lz
        self.nx = nx
        self.ny = ny
        self.nz = nz

    def get(self):
        return self.x0, self.y0, self.z0, self.lx, self.ly, self.lz, self.nx, self.ny, self.nz

class MeshOptions:
    class Option:
        def initialize(self): pass

        def __init__(self):
            self.mesh = None
            self.dbc_locator = None
            self.nbc_value = None
            self.cell_markers = None
            self.facet_markers = None
            self.mesh_params = None

    class _BarOption(Option):
        def initialize(self):
            super().__init__()
            mesh_data = io.gmsh.read_from_msh('bar3d.msh', MPI.COMM_WORLD, gdim=3)
            self.mesh = mesh_data.mesh
            self.cell_markers = mesh_data.cell_tags
            self.facet_markers = mesh_data.facet_tags
            self.dbc_locator = lambda x: np.isclose(x[0], 0)
            self.nbc_value = -0.01
            self.mesh_params = MeshParams(0.0, 0.0, 0.0, 1.0, 0.25, 0.25, 5, 2, 2)

    class _NotchOption(Option):
        def initialize(self):
            super().__init__()
            mesh_data = io.gmsh.read_from_msh('notch.msh', MPI.COMM_WORLD, gdim=3)

            # get bounds and shift lower corner to origin
            low = np.min(mesh_data.mesh.geometry.x, axis=0)
            high = np.max(mesh_data.mesh.geometry.x, axis=0)
            mesh_data.mesh.geometry.x[:] += np.abs(low)

            # fix facet markers for nm bc using facet center of mass (com) locations
            facet_markers = mesh_data.facet_tags
            facet_indices = facet_markers.indices.copy()
            facet_values = facet_markers.values.copy()
            com = mesh.compute_midpoints(mesh_data.mesh, 2, facet_markers.indices)
            com_y = com[:, 1]
            facet_values[np.isclose(com_y, 110, atol=1e-8)] = 3
            updated_facet_markers = mesh.meshtags(mesh_data.mesh, 2, facet_indices, facet_values)

            self.mesh = mesh_data.mesh
            self.cell_markers = mesh_data.cell_tags
            self.facet_markers = updated_facet_markers
            self.dbc_locator = lambda x: np.isclose(x[1], 0)
            self.nbc_value = 0.01
            dims = np.abs(high - low)
            self.mesh_params = MeshParams(0.0, 0.0, 0.0, dims[0], dims[1], dims[2])

    BAR: Option = _BarOption()
    NOTCH: Option = _NotchOption()

class FenicsXWrapper:
    def __init__(self, mesh: MeshOptions.Option):
        mesh.initialize()
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

        # domain
        self.domain = mesh.mesh
        self.cell_markers = mesh.cell_markers
        self.facet_markers = mesh.facet_markers
        self.ds = ufl.Measure('ds', domain=self.domain, subdomain_data=self.facet_markers)

        # material properties
        self.lam   = fem.Constant(self.domain, 10.0)
        self.mu    = fem.Constant(self.domain, 5.0)
        self.alpha = fem.Constant(self.domain, 100.0)

        # FEM
        self.V = fem.functionspace(self.domain, ("CG", 2, (3,)))
        self.uh = fem.Function(self.V)
        self.uh.name = "displacement"
        self.u = ufl.TrialFunction(self.V)
        self.v = ufl.TestFunction(self.V)

        self.W = fem.functionspace(self.domain, ("CG", 2, (6,)))
        self.WT = fem.functionspace(self.domain, ("CG", 2, (6, 6)))
        self.W3 = fem.functionspace(self.domain, ("CG", 2, (3,)))

        self.eps = ufl.variable(FenicsXWrapper.symgrad_mandel(self.uh))
        self.sig = fem.Function(self.W)
        self.tan = fem.Function(self.WT)
        self.vm_stress = fem.Function(self.V)
        self.vm_stress.name = "von_mises_stress"

        # BCs
        self.dc_bcs = [self.d_bc(fem.Constant(self.domain, 0.), self.V.sub(i), mesh.dbc_locator) for i in range(3)]
        self.bc_nm = ufl.inner(mesh.nbc_value, self.v[1]) * self.ds(3)

        # variational form
        self.res = ufl.inner(self.sig, FenicsXWrapper.symgrad_mandel(self.v)) * ufl.dx - self.bc_nm(self.v)
        self.jac = ufl.inner(ufl.dot(self.tan, FenicsXWrapper.symgrad_mandel(self.u)), FenicsXWrapper.symgrad_mandel(self.v)) * ufl.dx

        # coupling buffers
        self.eps_buffers = FenicsXWrapper.gen_coupling_buffers(self.W3, 2)
        self.sig_buffers = FenicsXWrapper.gen_coupling_buffers(self.W3, 2)
        self.tan_buffers = FenicsXWrapper.gen_coupling_buffers(self.W3, 7)

    def init_pure_meso(self):
        psi = self.calc_psi()
        sigma = ufl.diff(psi, self.eps)
        tangent = ufl.diff(sigma, self.eps)
        init_res = ufl.inner(sigma, FenicsXWrapper.symgrad_mandel(self.v)) * ufl.dx - self.bc_nm(self.v)
        init_jac = ufl.inner(ufl.dot(tangent, FenicsXWrapper.symgrad_mandel(self.u)), FenicsXWrapper.symgrad_mandel(self.v)) * ufl.dx
        problem = NonlinearProblem(init_res, self.uh, bcs=self.dc_bcs, J=init_jac, petsc_options=self.petsc_options,
                                   petsc_options_prefix="nonlinpoisson")
        problem.solve()

        sig_eval = self.eval_var(ufl.variable(sigma), self.W)
        self.sig.x.array[:] = sig_eval.var_val.x.array[:]

        tan_eval = self.eval_var(ufl.variable(tangent), self.WT)
        self.tan.x.array[:] = tan_eval.var_val.x.array[:]

    def solve(self):
        a = ufl.inner(ufl.dot(self.tan, FenicsXWrapper.symgrad_mandel(self.u)), FenicsXWrapper.symgrad_mandel(self.v)) * ufl.dx
        L = self.bc_nm
        problem = LinearProblem(
            a=a,
            L=L,
            bcs=self.dc_bcs,
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
            petsc_options_prefix="Poisson",
        )
        result = problem.solve()
        self.uh.x.array[:] = result.x.array[:]

    def calc_psi(self):
        tr_e = self.eps[0] + self.eps[1] + self.eps[2]
        e2 = ufl.inner(self.eps, self.eps)
        tr_e2 = tr_e ** 2
        half_alpha = 0.5 * self.alpha
        return 0.5 * self.lam * (1.0 + half_alpha * tr_e2) * tr_e2 + self.mu * (1 + half_alpha * e2) * e2

    def d_bc(self, bc_val, V: fem.FunctionSpace, loc_fun: Callable = lambda x: np.isclose(x[0], 0)):
        fdim = self.domain.topology.dim - 1
        b_facets = dolfinx.mesh.locate_entities_boundary(self.domain, fdim, loc_fun)
        b_dofs = fem.locate_dofs_topological(V, fdim, b_facets)
        return fem.dirichletbc(bc_val, b_dofs, V)

    def calc_von_mises_stress(self):
        sig = self.sig.x.array.reshape(-1, 6)
        s1 = 0.5 * (np.power(sig[:, 0] - sig[:, 1], 2) +
                    np.power(sig[:, 1] - sig[:, 2], 2) +
                    np.power(sig[:, 2] - sig[:, 0], 2))
        # we div by 2 to remove sqrt 2 from sig
        # accounted by factor 3 -> 1.5
        s2 = 1.5 * (np.power(sig[:, 3], 2) + np.power(sig[:, 4], 2) + np.power(sig[:, 5], 2))
        self.vm_stress.x.array.reshape(-1, 3)[:, 0] = np.sqrt(s1 + s2)

    @staticmethod
    def get_buffer_funcs(buffers):
        buffer_list, _ = buffers
        return buffer_list

    @staticmethod
    def symgrad_mandel(vec):
        halfsqrt2 = 0.5 * np.sqrt(2)
        return ufl.as_vector([vec[0].dx(0), vec[1].dx(1), vec[2].dx(2),
                              halfsqrt2 * (vec[1].dx(2) + vec[2].dx(1)),
                              halfsqrt2 * (vec[0].dx(2) + vec[2].dx(0)),
                              halfsqrt2 * (vec[0].dx(1) + vec[1].dx(0))])

    @staticmethod
    def eval_var(var, space):
        e = Evaluator(var, space)
        e.interpolate()
        return e

    @staticmethod
    def gen_coupling_buffers(space: fem.FunctionSpace, count: int):
        buffers = [fem.Function(space) for _ in range(count)]
        return buffers, space

    @staticmethod
    def write_to_buffers(buffers, values):
        buffer_list, space = buffers
        f_dim = space.value_size

        assert len(buffer_list) == len(values)
        for func, v_list in zip(buffer_list, values):
            func.x.array.reshape(-1, f_dim)[:, :] = np.array(v_list)

    @staticmethod
    def copy_from_buffers(targets, buffers):
        buffer_list, space = buffers
        f_dim = space.value_size

        assert len(buffer_list) == len(targets)
        for func, tgt in zip(buffer_list, targets):
            tgt[:, :] = func.x.array.reshape(-1, f_dim)[:, :]

    @staticmethod
    def copy_to_buffers(sources, buffers):
        buffer_list, space = buffers
        f_dim = space.value_size

        assert len(buffer_list) == len(sources)
        for func, src in zip(buffer_list, sources):
            func.x.array.reshape(-1, f_dim)[:, :] = src[:, :]

class DataTransformer:
    ADAPTER_OR_SURROGATE = 0
    PYFANS = 1
    NASMAT = 2

    def __init__(self, type):
        """
        type: either ADAPTER_OR_SURROGATE, PYFANS or NASMAT
        """
        self.sig_impl = DataTransformer._noop
        self.eps_impl = DataTransformer._noop
        self.tan_impl = DataTransformer._noop

        if type == self.ADAPTER_OR_SURROGATE:
            self.sig_impl = DataTransformer._noop
            self.eps_impl = DataTransformer._noop
            self.tan_impl = DataTransformer._noop
        if type == self.PYFANS:
            self.sig_impl = DataTransformer._handle_sig_eps_fans
            self.eps_impl = DataTransformer._handle_sig_eps_fans
            self.tan_impl = DataTransformer._handle_tan_fans
        if type == self.NASMAT:
            self.sig_impl = DataTransformer._handle_sig_nasmat
            self.eps_impl = DataTransformer._handle_eps_nasmat
            self.tan_impl = DataTransformer._handle_tan_nasmat

    def handle_sig(self, sig_buffer): pass
    def handle_eps(self, eps_buffer): pass
    def handle_tan(self, tan_buffer): pass

    @staticmethod
    def _swap_arr(a, b):
        tmp = a[:].copy()
        a[:] = b[:]
        b[:] = tmp[:]

    @staticmethod
    def _noop(dummy_arg): pass

    @staticmethod
    def _handle_sig_eps_fans(coupling_buffer):
        buffers, _ = coupling_buffer
        sig_high = buffers[1].x.array.reshape(-1, 3)
        DataTransformer._swap_arr(sig_high[:, 0], sig_high[:, 2])

    @staticmethod
    def _handle_tan_fans(tan_buffer):
        buffers, _ = tan_buffer
        #tan1 = buffers[0].x.array.reshape(-1, 3)
        tan2 = buffers[1].x.array.reshape(-1, 3)
        tan3 = buffers[2].x.array.reshape(-1, 3)
        tan4 = buffers[3].x.array.reshape(-1, 3)
        tan5 = buffers[4].x.array.reshape(-1, 3)
        tan6 = buffers[5].x.array.reshape(-1, 3)
        tan7 = buffers[6].x.array.reshape(-1, 3)

        DataTransformer._swap_arr(tan2[:, 0], tan2[:, 2])
        DataTransformer._swap_arr(tan3[:, 2], tan4[:, 1])
        DataTransformer._swap_arr(tan5[:, 0], tan5[:, 2])
        DataTransformer._swap_arr(tan6[:, 0], tan7[:, 2])
        DataTransformer._swap_arr(tan6[:, 1], tan7[:, 1])
        DataTransformer._swap_arr(tan6[:, 2], tan7[:, 0])

    @staticmethod
    def _handle_sig_nasmat(sig_buffer): pass
    @staticmethod
    def _handle_eps_nasmat(eps_buffer): pass
    @staticmethod
    def _handle_tan_nasmat(tan_buffer): pass

class MesoProblem:
    def __init__(self, data_transformer: DataTransformer, mesh_type: MeshOptions.Option = MeshOptions.NOTCH):
        # domain_path = None, mesh_params = Default_Mesh
        #if domain_path is None:
        #    domain_path = "bar3d.msh"
        #    MesoProblem.generate_bar_mesh(domain_path, mesh_params)
        self.mesh = mesh_type
        self.fnx = FenicsXWrapper(mesh_type)
        self.eps_eval = FenicsXWrapper.eval_var(self.fnx.eps, self.fnx.W)
        self.tan_conv_buffer = np.zeros((self.fnx.sig.x.array.reshape(-1, 6).shape[0], 21))
        self.tensor_mapping = np.array([
            [0,  1,  2,  3,  4,  5],
            [1,  6,  7,  8,  9, 10],
            [2,  7, 11, 12, 13, 14],
            [3,  8, 12, 15, 16, 17],
            [4,  9, 13, 16, 18, 19],
            [5, 10, 14, 17, 19, 20],
        ])
        self.data_transformer = data_transformer

    def split_eps_data(self):
        self.eps_eval.interpolate()
        self.eps_eval.var_val.x.petsc_vec.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        tmp_view = self.eps_eval.var_val.x.array.reshape(-1, 6)
        FenicsXWrapper.copy_to_buffers([tmp_view[:, 0:3], tmp_view[:, 3:6]], self.fnx.eps_buffers)
        self.data_transformer.handle_eps(self.fnx.eps_buffers)

    def merge_sig_data(self):
        self.data_transformer.handle_sig(self.fnx.sig_buffers)
        tmp_view = self.fnx.sig.x.array.reshape(-1, 6)
        FenicsXWrapper.copy_from_buffers([tmp_view[:, 0:3], tmp_view[:, 3:6]], self.fnx.sig_buffers)
        self.fnx.sig.x.petsc_vec.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

    def merge_and_conv_tan_data(self):
        self.data_transformer.handle_tan(self.fnx.tan_buffers)
        FenicsXWrapper.copy_from_buffers([self.tan_conv_buffer[:, 3*i:3*(i+1)] for i in range(7)], self.fnx.tan_buffers)
        self.fnx.tan.x.array.reshape(-1, 6, 6)[:, :, :] = self.as_sym_tensor_6x6(self.tan_conv_buffer)
        self.fnx.tan.x.petsc_vec.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

    def calc_von_mises_stress(self):
        self.fnx.calc_von_mises_stress()
        self.fnx.vm_stress.x.petsc_vec.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

    def as_sym_tensor_6x6(self, a):
        return a[:, self.tensor_mapping]

    @staticmethod
    def coupling_bc(x):
        return np.logical_or.reduce(np.equal(x, x))

    @staticmethod
    def generate_bar_mesh(msh_file, params: MeshParams, arrangement='AlternateRight'):
        x0, y0, z0, lx, ly, lz, nx, ny, nz = params.get()
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)  # to disable meshing info
        geom = gmsh.model.geo

        p = []
        p.append(geom.add_point(x0, y0, z0))
        p.append(geom.add_point(x0 + lx, y0, z0))
        p.append(geom.add_point(x0 + lx, y0 + ly, z0))
        p.append(geom.add_point(x0, y0 + ly, z0))
        p.append(geom.add_point(x0, y0, z0 + lz))
        p.append(geom.add_point(x0 + lx, y0, z0 + lz))
        p.append(geom.add_point(x0 + lx, y0 + ly, z0 + lz))
        p.append(geom.add_point(x0, y0 + ly, z0 + lz))

        l = []
        l.append(geom.add_line(p[0], p[1]))
        l.append(geom.add_line(p[1], p[2]))
        l.append(geom.add_line(p[2], p[3]))
        l.append(geom.add_line(p[3], p[0]))
        l.append(geom.add_line(p[0 + 4], p[1 + 4]))
        l.append(geom.add_line(p[1 + 4], p[2 + 4]))
        l.append(geom.add_line(p[2 + 4], p[3 + 4]))
        l.append(geom.add_line(p[3 + 4], p[0 + 4]))
        l.append(geom.add_line(p[0], p[0 + 4]))
        l.append(geom.add_line(p[1], p[1 + 4]))
        l.append(geom.add_line(p[2], p[2 + 4]))
        l.append(geom.add_line(p[3], p[3 + 4]))

        ll0 = geom.add_curve_loop([l[0], l[1], l[2], l[3]])  # front
        ll1 = geom.add_curve_loop([-l[4], -l[7], -l[6], -l[5]])  # back
        ll2 = geom.add_curve_loop([-l[1], l[9], l[5], -l[10]])  # right
        ll3 = geom.add_curve_loop([-l[3], l[11], l[7], -l[8]])  # left
        ll4 = geom.add_curve_loop([-l[0], l[8], l[4], -l[9]])  # bot
        ll5 = geom.add_curve_loop([-l[2], l[10], l[6], -l[11]])  # top
        ll = [ll0, ll1, ll2, ll3, ll4, ll5]
        s = [geom.add_plane_surface([lli]) for lli in ll]
        geom.synchronize()
        sl = geom.add_surface_loop(s)
        v = geom.add_volume([sl])

        geom.synchronize()

        for li, ni in zip(l, [nx, ny, nx, ny, nx, ny, nx, ny, nz, nz, nz, nz]):
            gmsh.model.mesh.set_transfinite_curve(li, ni + 1)
        for si in s:
            gmsh.model.mesh.set_transfinite_surface(si, arrangement=arrangement)
        gmsh.model.mesh.set_transfinite_volume(v)

        gmsh.model.add_physical_group(3, [v], 0)  # tag 0, full domain
        for i in range(6):
            gmsh.model.add_physical_group(2, [ll[i]], i + 1)  # tag 1 to 6, faces

        gmsh.model.mesh.generate(dim=3)
        gmsh.write(msh_file)
        gmsh.finalize()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Runs the macro simulation using fenicsx for the notch or bar geometry.',
        usage='%(prog)s [options]'
    )
    parser.add_argument(
        '--case',
        default='bar',
        choices=['bar', 'notch'],
        type=str,
        help='The case type.'
    )
    parser.add_argument(
        '--micro',
        default='ADA',
        choices=['ADA', 'pyFANS', 'NASMAT'],
        type=str,
        help='The micro solver type.'
    )
    args = parser.parse_args()


    transform_type = None
    if args.micro == 'ADA': transform_type = DataTransformer.ADAPTER_OR_SURROGATE
    if args.micro == 'pyFANS': transform_type = DataTransformer.PYFANS
    if args.micro == 'NASMAT': transform_type = DataTransformer.NASMAT

    mesh_type = None
    if args.case == 'bar': mesh_type = MeshOptions.BAR
    if args.case == 'notch': mesh_type = MeshOptions.NOTCH

    data_transformer = DataTransformer(transform_type)
    mp = MesoProblem(data_transformer, mesh_type)
    mesh_params = mp.mesh.mesh_params

    mp.fnx.init_pure_meso()         # "populate" all fields with values
    mp.split_eps_data()             # write data to coupling buffers

    # set up coupling
    precice = Adapter(MPI.COMM_WORLD, 'precice-adapter-config.json')
    coupling_boundary = MesoProblem.coupling_bc # we are coupling on full domain
    access_region = [(mesh_params.x0, mesh_params.y0, mesh_params.z0),
                     (mesh_params.lx + 1e-14, mesh_params.ly + 1e-14, mesh_params.lz + 1e-14)]
    read_fields = {
        'stresses1to3': mp.fnx.W3,
        'stresses4to6': mp.fnx.W3,
        'cmat1'       : mp.fnx.W3,
        'cmat2'       : mp.fnx.W3,
        'cmat3'       : mp.fnx.W3,
        'cmat4'       : mp.fnx.W3,
        'cmat5'       : mp.fnx.W3,
        'cmat6'       : mp.fnx.W3,
        'cmat7'       : mp.fnx.W3
    }
    write_fields = {
        'strains1to3': FenicsXWrapper.get_buffer_funcs(mp.fnx.eps_buffers)[0],
        'strains4to6': FenicsXWrapper.get_buffer_funcs(mp.fnx.eps_buffers)[1]
    }
    coupling_mesh = CouplingMesh('meso-mesh', coupling_boundary, read_fields, write_fields)
    precice.initialize([coupling_mesh])
    
    is_coupling_ongoing, t, n = precice.is_coupling_ongoing(), 0.0, 0
    while is_coupling_ongoing:
        if precice.requires_writing_checkpoint():
            state = (mp.fnx.sig.x.array, mp.eps_eval.var_val.x.array, mp.fnx.uh.x.array)
            precice.store_checkpoint(state, t, n)
        dt = fem.Constant(mp.fnx.domain, precice.get_max_time_step_size())

        values = [list(precice.read_data("meso-mesh", name, dt).values()) for name in ['stresses1to3', 'stresses4to6']]
        FenicsXWrapper.write_to_buffers(mp.fnx.sig_buffers, values)
        mp.merge_sig_data()

        values = [list(precice.read_data("meso-mesh", name, dt).values()) for name in ['cmat1', 'cmat2', 'cmat3', 'cmat4', 'cmat5', 'cmat6', 'cmat7']]
        FenicsXWrapper.write_to_buffers(mp.fnx.tan_buffers, values)
        mp.merge_and_conv_tan_data()

        mp.fnx.solve()

        is_coupling_ongoing = precice.is_coupling_ongoing()
        if is_coupling_ongoing:
            mp.split_eps_data()
            for eps_, name in zip(FenicsXWrapper.get_buffer_funcs(mp.fnx.eps_buffers), ['strains1to3', 'strains4to6']):
                precice.write_data("meso-mesh", name, eps_)
            precice.advance(dt(0))
        is_coupling_ongoing = precice.is_coupling_ongoing()

        if precice.requires_reading_checkpoint():
            state, t_, n_ = precice.retrieve_checkpoint()
            sig_arr, eps_arr, uh_arr = state
            mp.fnx.sig.x.array[:] = sig_arr
            mp.eps_eval.var_val.x.array[:] = eps_arr
            mp.fnx.uh.x.array[:] = uh_arr
            t = t_
            n = n_
        else:
            t += float(dt)
            n += 1

        mp.calc_von_mises_stress()
        with dolfinx.io.VTXWriter(MPI.COMM_WORLD, f"multi_scale_{n}.bp", [mp.fnx.uh, mp.fnx.vm_stress], engine="BP4") as vtx:
            vtx.write(t)

    precice.finalize()
    print(np.linalg.norm(mp.fnx.uh.x.array))