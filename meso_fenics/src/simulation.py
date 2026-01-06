# pure meso, pseudo coupled, coupled

from .config import Config
from .fnx import MesoProblem, MultiscaleProblem, Evaluator
from .meshes import Mesh
from .util import Registry
from .coupling import CouplingBuffer, DataTransformer, Projectors, Mergers

from mpi4py import MPI
from fenicsxprecice import Adapter, CouplingMesh
from dolfinx import io, fem

import numpy as np
import h5py

class Simulation:
    TYPES = Registry()
    def __init__(self, config: Config):
        self.output_path = config.output_path
        self.mesh = Mesh()
        self.mesh.load(config)

        self.problem = None
        self._write_state = config.simulation_write_state
        self._write_state_type = config.simulation_write_state_type

    def run(self): pass

    def write_output(self, t, n):
        self.problem.calc_von_mises_stress()
        with io.VTXWriter(MPI.COMM_WORLD, f"{self.output_path}_{n}.bp", [self.problem.uh, self.problem.vm_stress_fun],
                          engine="BP4") as vtx:
            vtx.write(t)

        if self._write_state is not None and type(self._write_state_type) == list:
            with h5py.File(self._write_state, "w") as f:
                if 'E' in self._write_state_type:
                    eval = Evaluator(self.problem.eps_var, self.problem.W)
                    eval.interpolate()
                    f.create_dataset("strain_data", data=eval.var_val.x.array)

                if 'S' in self._write_state_type:
                    f.create_dataset("stress_data", data=self.problem.sig_fun.x.array)

                if 'U' in self._write_state_type:
                    f.create_dataset("displacement_data", data=self.problem.uh.x.array)

    @staticmethod
    def generate(config: Config):
        type = config.simulation_type
        if not Simulation.TYPES.is_name_valid(type):
            options = "["
            names = Simulation.TYPES.get_names()
            for name in names[:-1]: options += name + ","
            options = options[-1] + "]"
            raise RuntimeError(f"Invalid simulation type. Select from: {options}")
        return Simulation.TYPES.get_by_name(type)(config)

@Simulation.TYPES.register
class MesoSim(Simulation):
    def __init__(self, config: Config):
        super().__init__(config)
        self.problem = MesoProblem(config, self.mesh)

    def run(self):
        self.problem.solve()
        self.write_output(0.0, 0)

@Simulation.TYPES.register
class PseudoCoupledSim(Simulation):
    def __init__(self, config: Config):
        super().__init__(config)
        self.input_path = config.simulation_input
        if self.input_path is None: raise RuntimeError("PseudoCoupledSim requires input path")

        self.problem = MultiscaleProblem(config, self.mesh)

    def run(self):
        with h5py.File(self.input_path, "r") as f:
            stress_data = f["stress_data"][:]
            self.problem.sig_fun.x.array[:] = stress_data[:]
            self.problem.sig_fun.x.scatter_forward()

            if self.problem.is_small_strain:
                tan_data = f["tan_data"][:]
                self.problem.tan_fun.x.array[:] = tan_data[:]
                self.problem.tan_fun.x.scatter_forward()

        self.problem.solve()
        self.write_output(0.0, 0)

@Simulation.TYPES.register
class CoupledSim(Simulation):
    def __init__(self, config: Config):
        super().__init__(config)
        self.problem = MultiscaleProblem(config, self.mesh)
        self.transform = DataTransformer(config.simulation_micro_type)

        num_sig_buffers = int(np.ceil(self.problem.W.value_size / 3.))
        self.sig_buffer = CouplingBuffer(
            self.problem.sig_fun,
            self.problem.W3,
            num_sig_buffers,
            Projectors.InplaceSplitter(self.problem.W3),
            Mergers.InplaceMerger()
        )
        self.eps_eval   = Evaluator(self.problem.eps_var, self.problem.W)
        self.eps_buffer = CouplingBuffer(
            self.eps_eval.var_val,
            self.problem.W3,
            num_sig_buffers,
            Projectors.InplaceSplitter(self.problem.W3),
            Mergers.InplaceMerger()
        )
        self.tan_buffer = None
        if self.problem.is_small_strain:
            self.tan_buffer = CouplingBuffer(
                self.problem.tan_fun,
                self.problem.W3,
                7,
                Projectors.SelectionSplitter(
                    self.problem.W3,
                    np.array([0, 1,  2,  3,  4,  5,
                                 7,  8,  9, 10, 11,
                                    14, 15, 16, 17,
                                        21, 22, 23,
                                            28, 29,
                                                35])),
                Mergers.SelectionMerger(
                    self.problem.W3,
                    7,
                    self.problem.tan_fun.x.array.reshape(-1, 6, 6).shape[0],
                    np.array([
                        [0,  1,  2,  3,  4,  5],
                        [1,  6,  7,  8,  9, 10],
                        [2,  7, 11, 12, 13, 14],
                        [3,  8, 12, 15, 16, 17],
                        [4,  9, 13, 16, 18, 19],
                        [5, 10, 14, 17, 19, 20],
                    ]))
            )
        else:
            # Need to remove transforms if using large strain
            self.transform.clear_transforms()

        self.problem.solve_meso() # init fields to set up coupling
        self.eps_eval.interpolate()
        self.eps_buffer.write_origin_to_buffer(self.transform.get_transform_eps())

        adapter_config_path = 'res/precice-adapter-config-small-strain.json'
        if not self.problem.is_small_strain:
            adapter_config_path = 'res/precice-adapter-config-large-strain.json'
        self.precice = Adapter(MPI.COMM_WORLD, adapter_config_path)
        coupling_boundary = CoupledSim.coupling_bc

        self.read_functions = dict()
        read_fields = dict()
        # attach sig functions
        for i, func in enumerate(self.sig_buffer.get_functions()):
            key = f"stresses{3*(i+0)+1}to{3*(i+0)+1+2}"
            self.read_functions[key] = func
            read_fields[key] = self.problem.W3
        # attach tan functions
        if self.tan_buffer is not None:
            for i, func in enumerate(self.tan_buffer.get_functions()):
                self.read_functions[f"cmat{i+1}"] = func

        self.write_fields = dict()
        # attach eps functions
        for i, func in enumerate(self.eps_buffer.get_functions()):
            key = f"strains{3*(i+0)+1}to{3*(i+0)+1+2}"
            self.write_fields[key] = func

        coupling_mesh = CouplingMesh('meso-mesh', coupling_boundary, read_fields, self.write_fields)
        self.precice.initialize([coupling_mesh])

    def write_checkpoint(self, t, n):
        if self.precice.requires_writing_checkpoint():
            self.precice.store_checkpoint(self.problem.uh, t, n)

    def read_checkpoint(self):
        if self.precice.requires_reading_checkpoint():
            state, t_, n_ = self.precice.retrieve_checkpoint()
            uh_arr = state
            self.problem.uh.x.array[:] = uh_arr.x.array[:]
            return True, t_, n_
        else:
            return False, None, None


    def run(self):
        is_coupling_ongoing, t, n = self.precice.is_coupling_ongoing(), 0.0, 0

        # handle first timestep: init MESO (was already computed)
        self.write_checkpoint(t, n)
        # no need to solve
        self.eps_eval.interpolate()
        self.eps_buffer.write_origin_to_buffer(self.transform.get_transform_eps())
        for name, func in self.write_fields.items(): self.precice.write_data("meso-mesh", name, func)
        dt = self.precice.get_max_time_step_size()
        self.precice.advance(dt)
        is_coupling_ongoing = self.precice.is_coupling_ongoing()
        loaded_checkpoint, t_, n_ = self.read_checkpoint()
        if loaded_checkpoint:
            t, n = t_, n_
        else:
            t += dt
            n += 1

        while is_coupling_ongoing:
            # Handle checkpointing
            self.write_checkpoint(t, n)
            dt = fem.Constant(self.mesh.domain, self.precice.get_max_time_step_size())

            # Read Data
            for name, func in self.read_functions.items(): self.precice.read_data("meso-mesh", name, dt, func)
            self.sig_buffer.write_buffer_to_origin(self.transform.get_transform_sig())
            self.sig_buffer.original.x.scatter_forward()
            if self.problem.is_small_strain:
                self.tan_buffer.write_buffer_to_origin(self.transform.get_transform_tan())
                self.tan_buffer.original.x.scatter_forward()

            # Solve EQ
            self.problem.solve()

            # Write Data
            is_coupling_ongoing = self.precice.is_coupling_ongoing()
            if is_coupling_ongoing:
                self.eps_eval.interpolate()
                self.eps_buffer.write_origin_to_buffer(self.transform.get_transform_eps())
                for name, func in self.write_fields.items(): self.precice.write_data("meso-mesh", name, func)
                self.precice.advance(dt(0))
            is_coupling_ongoing = self.precice.is_coupling_ongoing()

            # Load Checkpoint?
            loaded_checkpoint, t_, n_ = self.read_checkpoint()
            if loaded_checkpoint:
                t, n = t_, n_
            else:
                t += float(dt)
                n += 1

            # Output
            self.write_output(t, n)

    @staticmethod
    def coupling_bc(x):
        return np.logical_or.reduce(np.equal(x, x))