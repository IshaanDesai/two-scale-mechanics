# pure meso, pseudo coupled, coupled

from .config import Config
from .fnx import MesoProblem, MultiscaleProblem, Evaluator
from .meshes import Mesh
from .util import Registry
from .coupling import CouplingBuffer, DataTransformer, Projectors, Mergers

from mpi4py import MPI
from fenicsxprecice import Adapter, CouplingMesh
from fenicsxprecice.adapter_core import (
    convert_fenicsx_to_precice,
    get_fenicsx_interpolation_points,
    interpolate_boundary_function,
    FunctionType,
)
from dolfinx import io, fem

import numpy as np
import h5py


class Simulation:
    """
    Base class for different types of simulations.

    This class provides a common interface for running simulations and
    writing output. Specific simulation types are registered and created
    via the factory pattern.

    Parameters
    ----------
    config : Config
        Configuration object containing simulation parameters.

    Attributes
    ----------
    TYPES : Registry
        Registry of available simulation types.
    output_path : str
        Path for output files.
    mesh : Mesh
        The computational mesh.
    problem : MesoProblem or MultiscaleProblem
        The problem solver instance.
    """

    TYPES = Registry()

    def __init__(self, config: Config):
        self.output_path = config.output_path
        self.mesh = Mesh()
        self.mesh.load(config)

        self.problem = None
        self._write_state = config.simulation_write_state
        self._write_state_type = config.simulation_write_state_type

        self._coords_W = None
        self._cells_W = None
        self._values_to_send_W = None
        self._coords_WT = None
        self._cells_WT = None
        self._values_to_send_WT = None

    def load_quad_coords(self):
        """
        Loads quadrature point coordinates into buffers.
        Only call once self.problem is populated.
        """
        (
            self._coords_W,
            self._cells_W,
            self._values_to_send_W,
        ) = get_fenicsx_interpolation_points(
            self.problem.W, Simulation.always_true, MPI.COMM_WORLD, 25
        )

        if self.problem.WT is not None:
            (
                self._coords_WT,
                self._cells_WT,
                self._values_to_send_WT,
            ) = get_fenicsx_interpolation_points(
                self.problem.WT, Simulation.always_true, MPI.COMM_WORLD, 25
            )

    def run(self):
        """
        Run the simulation.

        This method should be implemented by subclasses.
        """
        pass

    @staticmethod
    def always_true(x):
        """
        Define function that returns True for any np array

        Parameters
        ----------
        x : numpy.ndarray
            Array of coordinates.

        Returns
        -------
        numpy.ndarray
            Boolean array that is True for all points.
        """
        return np.logical_or.reduce(np.equal(x, x))

    def write_output(self, t, n, i=None):
        """
        Write simulation output to file.

        Parameters
        ----------
        t : float
            Current simulation time.
        n : int
            Current time step number.
        i : int, optional
            Optional time step iteration
        """
        self.problem.calc_von_mises_stress()
        iter = ""
        if i is not None:
            iter = f"_{i}"
        rank = ""
        if MPI.COMM_WORLD.Get_size() > 1:
            rank = f"_rank{MPI.COMM_WORLD.Get_rank()}"
        with io.VTXWriter(
            MPI.COMM_WORLD,
            f"{self.output_path}_{n}{iter}.bp",
            [self.problem.uh, self.problem.vm_stress_fun],
            engine="BP4",
        ) as vtx:
            vtx.write(t)

        if self._write_state is not None and type(self._write_state_type) == list:
            with h5py.File(f"{self._write_state}_{n}{iter}{rank}.h5", "w") as f:
                if "E" in self._write_state_type and self._coords_W is not None:
                    eval = Evaluator(self.problem.eps_var, self.problem.W)
                    eval.interpolate()
                    eps = convert_fenicsx_to_precice(eval.var_val, self._coords_W, 25)
                    f.create_dataset("strain_data", data=eps)

                if "S" in self._write_state_type and self._coords_W is not None:
                    sig = convert_fenicsx_to_precice(
                        self.problem.sig_fun, self._coords_W, 25
                    )
                    f.create_dataset("stress_data", data=sig)

                if (
                    "T" in self._write_state_type
                    and self._coords_WT is not None
                    and self.problem.tan_fun is not None
                ):
                    tan = convert_fenicsx_to_precice(
                        self.problem.tan_fun, self._coords_WT, 25
                    )
                    f.create_dataset("tangent_data", data=tan)

                if "U" in self._write_state_type:
                    f.create_dataset("displacement_data", data=self.problem.uh.x.array)

    @staticmethod
    def generate(config: Config):
        """
        Factory method to create a simulation instance based on configuration.

        Parameters
        ----------
        config : Config
            Configuration object specifying the simulation type.

        Returns
        -------
        Simulation
            An instance of the appropriate Simulation subclass.

        Raises
        ------
        RuntimeError
            If the simulation type is invalid.
        """
        type = config.simulation_type
        if not Simulation.TYPES.is_name_valid(type):
            options = "["
            names = Simulation.TYPES.get_names()
            for name in names[:-1]:
                options += name + ","
            options = options[-1] + "]"
            raise RuntimeError(f"Invalid simulation type. Select from: {options}")
        return Simulation.TYPES.get_by_name(type)(config)


@Simulation.TYPES.register
class MesoSim(Simulation):
    """
    Meso-scale only simulation.

    Solves a single-scale problem using the built-in constitutive law.

    Parameters
    ----------
    config : Config
        Configuration object containing simulation parameters.
    """

    def __init__(self, config: Config):
        super().__init__(config)
        self.problem = MesoProblem(config, self.mesh)
        self.load_quad_coords()

    def run(self):
        """
        Run the meso-scale simulation.

        Solves the problem once and writes output.
        """
        self.problem.solve()
        self.write_output(0.0, 0)


@Simulation.TYPES.register
class PseudoCoupledSim(Simulation):
    """
    Pseudo-coupled simulation using pre-computed micro-scale data.

    Loads stress and tangent modulus from a file and solves the
    meso-scale problem using this external constitutive response.

    Parameters
    ----------
    config : Config
        Configuration object containing simulation parameters.

    Raises
    ------
    RuntimeError
        If input path is not provided in the configuration.
    """

    def __init__(self, config: Config):
        super().__init__(config)
        self.input_path = config.simulation_input
        if self.input_path is None:
            raise RuntimeError("PseudoCoupledSim requires input path")

        self.problem = MultiscaleProblem(config, self.mesh)
        self.load_quad_coords()

    def run(self):
        """
        Run the pseudo-coupled simulation.

        Loads stress and tangent data from file, solves the problem,
        and writes output.
        """
        with h5py.File(self.input_path, "r") as f:
            stress_data = f["stress_data"][:].reshape(
                -1, self.problem.W.value_size
            )  # n_quad x stress_size
            interpolate_boundary_function(
                {tuple(k): v for k, v in zip(self._coords_W, stress_data)},
                FunctionType.VECTOR,
                self._values_to_send_W,
                self.problem.sig_fun,
                self._cells_W,
                MPI.COMM_WORLD,
                False,
                25,
            )
            self.problem.sig_fun.x.scatter_forward()

            if self.problem.is_small_strain:
                tan_data = f["tan_data"][:].reshape(
                    -1, self.problem.WT.value_size
                )  # n_quad x tan_size
                interpolate_boundary_function(
                    {tuple(k): v for k, v in zip(self._coords_WT, tan_data)},
                    FunctionType.VECTOR,
                    self._values_to_send_WT,
                    self.problem.tan_fun,
                    self._cells_WT,
                    MPI.COMM_WORLD,
                    False,
                    25,
                )
                self.problem.tan_fun.x.scatter_forward()

        self.problem.solve()
        self.write_output(0.0, 0)


@Simulation.TYPES.register
class CoupledSim(Simulation):
    """
    Fully coupled multiscale simulation using preCICE.

    Couples the meso-scale problem with external micro-scale solvers
    via the preCICE coupling library.

    Parameters
    ----------
    config : Config
        Configuration object containing simulation parameters.

    Attributes
    ----------
    transform : DataTransformer
        Handles data format transformations for different micro-solvers.
    sig_buffer : CouplingBuffer
        Buffer for stress data exchange.
    eps_buffer : CouplingBuffer
        Buffer for strain data exchange.
    tan_buffer : CouplingBuffer or None
        Buffer for tangent modulus exchange (small strain only).
    precice : Adapter
        preCICE adapter for coupling.
    """

    def __init__(self, config: Config):
        super().__init__(config)
        self.problem = MultiscaleProblem(config, self.mesh)
        self.load_quad_coords()
        self.transform = DataTransformer(config.simulation_micro_type)

        (
            self.sig_buffer,
            self.eps_buffer,
            self.eps_eval,
            self.tan_buffer,
        ) = CoupledSim._construct_coupling_buffers(self.problem, self.transform)

        self.init_with_micro = config.simulation_init_with_micro
        if self.init_with_micro and not self.problem.is_small_strain:
            raise RuntimeError(
                "Large strain must initialize with meso values. Set init_with_micro to false!"
            )
        self.adapter_config_path = CoupledSim._get_adapter_path(
            self.problem.is_small_strain,
            config.simulation_precice_xml_path,
            config.simulation_slurm_id
        )
        self.precice = Adapter(MPI.COMM_WORLD, self.adapter_config_path)
        coupling_boundary = CoupledSim.coupling_bc

        (
            self.read_functions,
            read_fields,
            self.write_fields,
        ) = CoupledSim._construct_coupling_dicts(
            self.problem, self.sig_buffer, self.eps_buffer, self.tan_buffer
        )

        coupling_mesh = CouplingMesh(
            "meso-mesh", coupling_boundary, read_fields, self.write_fields
        )
        self.precice.initialize([coupling_mesh])

    @staticmethod
    def _construct_coupling_buffers(problem, transform):
        """
        Constructs coupling buffers for preCICE coupling.

        Parameters
        ----------
        problem : MultiscaleProblem
            Problem object containing multi-scale fenics problem

        transform : DataTransformer
            Data transformer object for coupling.
        """
        num_sig_buffers = int(np.ceil(problem.W.value_size / 3.0))
        sig_buffer = CouplingBuffer(
            problem.sig_fun,
            problem.W3,
            num_sig_buffers,
            Projectors.InplaceSplitter(problem.W3),
            Mergers.InplaceMerger(),
        )
        eps_eval = Evaluator(problem.eps_var, problem.W)
        eps_buffer = CouplingBuffer(
            eps_eval.var_val,
            problem.W3,
            num_sig_buffers,
            Projectors.InplaceSplitter(problem.W3),
            Mergers.InplaceMerger(),
        )
        tan_buffer = None
        if problem.is_small_strain:
            tan_buffer = CouplingBuffer(
                problem.tan_fun,
                problem.W3,
                7,
                Projectors.SelectionSplitter(
                    problem.W3,
                    np.array(
                        [
                            0,
                            1,
                            2,
                            3,
                            4,
                            5,
                            7,
                            8,
                            9,
                            10,
                            11,
                            14,
                            15,
                            16,
                            17,
                            21,
                            22,
                            23,
                            28,
                            29,
                            35,
                        ]
                    ),
                ),
                Mergers.SelectionMerger(
                    problem.W3,
                    7,
                    problem.tan_fun.x.array.reshape(-1, 6, 6).shape[0],
                    np.array(
                        [
                            [0, 1, 2, 3, 4, 5],
                            [1, 6, 7, 8, 9, 10],
                            [2, 7, 11, 12, 13, 14],
                            [3, 8, 12, 15, 16, 17],
                            [4, 9, 13, 16, 18, 19],
                            [5, 10, 14, 17, 19, 20],
                        ]
                    ),
                ),
            )
        else:
            # Need to remove transforms if using large strain
            transform.clear_transforms()

        return sig_buffer, eps_buffer, eps_eval, tan_buffer

    @staticmethod
    def _get_adapter_path(is_small_strain, prc_path, run_id):
        """
        Returns path to preCICE adapter config file based on simulation parameters.

        Parameters
        ----------
        is_small_strain : bool
            Is current simulation running in small strain mode
        prc_path : str
            Path to preCICE adapter config file
        run_id : str
            Unique run id to split potentially parallel jobs

        Returns
        -------
        config_path : str
            Path to preCICE adapter config file
        """
        strain_type = "small" if is_small_strain else "large"
        adapter_config_path_base = f"res/precice-adapter-config-{strain_type}-strain.json"
        adapter_config_path_run  = f"res/precice-adapter-config-{strain_type}-strain-{run_id}.json"

        # create job specific adapter config
        if MPI.COMM_WORLD.Get_rank() == 0:
            try:
                import json
                with open(adapter_config_path_base, "r") as infile:
                    base_conf = json.load(infile)

                if prc_path is not None:
                    base_conf["config_file_name"] = prc_path

                with open(adapter_config_path_run, "w") as outfile:
                    json.dump(base_conf, outfile, indent=4)
            except Exception as e:
                print(e)
                raise RuntimeError("Failed to modify adapter config file")
        MPI.COMM_WORLD.Barrier()

        return adapter_config_path_run

    @staticmethod
    def _delete_adapter_file(path):
        if MPI.COMM_WORLD.Get_rank() != 0: return

        import os
        if os.path.exists(path):
            os.remove(path)

    @staticmethod
    def _construct_coupling_dicts(problem, sig_buffer, eps_buffer, tan_buffer):
        """
        Constructs coupling dicts for preCICE coupling.

        Parameters
        ----------
        problem : MultiscaleProblem
            Problem object containing multi-scale fenics problem
        sig_buffer : CouplingBuffer
            stress buffer
        eps_buffer : CouplingBuffer
            strain buffer
        tan_buffer : Optional[CouplingBuffer]
            tangent buffer

        Returns
        -------
        read_functions : dict
            Dictionary of coupling read functions
        read_fields : dict
            Dictionary of coupling read fields. Same as read_functions["key"].function_space
        write_fields : dict
            Dictionary of coupling write fields
        """
        read_functions = dict()
        read_fields = dict()
        # attach sig functions
        for i, func in enumerate(sig_buffer.get_functions()):
            key = f"stresses{3 * (i + 0) + 1}to{3 * (i + 0) + 1 + 2}"
            read_functions[key] = func
            read_fields[key] = problem.W3
        # attach tan functions
        if tan_buffer is not None:
            for i, func in enumerate(tan_buffer.get_functions()):
                read_functions[f"cmat{i + 1}"] = func

        write_fields = dict()
        # attach eps functions
        for i, func in enumerate(eps_buffer.get_functions()):
            key = f"strains{3 * (i + 0) + 1}to{3 * (i + 0) + 1 + 2}"
            write_fields[key] = func

        return read_functions, read_fields, write_fields

    def write_checkpoint(self, t, n):
        """
        Write a checkpoint if required by preCICE.

        Parameters
        ----------
        t : float
            Current simulation time.
        n : int
            Current time step number.
        """
        if self.precice.requires_writing_checkpoint():
            self.precice.store_checkpoint(self.problem.uh, t, n)

    def read_checkpoint(self, t, n):
        """
        Read a checkpoint if required by preCICE.

        Returns
        -------
        loaded : bool
            True if checkpoint was loaded, False otherwise.
        t : float or None
            Restored simulation time if checkpoint loaded.
        n : int or None
            Restored time step number if checkpoint loaded.
        """
        if self.precice.requires_reading_checkpoint():
            state, t_, n_ = self.precice.retrieve_checkpoint()
            #uh_arr = state
            #self.problem.uh.x.array[:] = uh_arr.x.array[:]
            return True, t_, n_
        else:
            return False, t, n

    def init_meso(self, t, n, dt):
        self.problem.solve_meso()
        self.coupling_write()
        self.write_output(t, n, "init-meso-sol")

    def init_micro(self, t, n, dt):
        self.coupling_write()  # write 0, only use tan result
        self.precice.advance(dt)
        self.read_checkpoint(t, n)
        self.write_checkpoint(t, n)
        self.coupling_read(dt)
        self.problem.init_with_ms()
        self.coupling_write()  # write eps based on tan
        self.write_output(t, n, "init-micro-sol")

    def coupling_write(self):
        self.eps_eval.interpolate()
        self.eps_buffer.write_origin_to_buffer(self.transform.get_transform_eps())
        for name, func in self.write_fields.items():
            self.precice.write_data("meso-mesh", name, func)

    def coupling_read(self, dt):
        for name, func in self.read_functions.items():
            self.precice.read_data("meso-mesh", name, dt, func)
        self.sig_buffer.write_buffer_to_origin(self.transform.get_transform_sig())
        self.sig_buffer.original.x.scatter_forward()
        if self.problem.is_small_strain:
            self.tan_buffer.write_buffer_to_origin(self.transform.get_transform_tan())
            self.tan_buffer.original.x.scatter_forward()

    def run(self):
        """
        Run the coupled multiscale simulation.

        Executes the coupling loop, exchanging data with micro-solvers
        via preCICE and solving the meso-scale problem at each iteration.
        """
        is_coupling_ongoing, t, n, it, init_done = (
            self.precice.is_coupling_ongoing(),
            0.0,
            0,
            0,
            False,
        )

        def update_counters(reloaded, t_, n_, it, dt):
            if reloaded:
                return t_, n_, (it + 1 if init_done else it)
            else:
                return t_ + float(dt), n_ + 1, (0 if init_done else it)

        self.write_checkpoint(t, n)
        dt = self.precice.get_max_time_step_size()
        if self.init_with_micro:
            self.init_micro(t, n, dt)
        else:
            self.init_meso(t, n, dt)
        self.precice.advance(dt)
        is_coupling_ongoing = self.precice.is_coupling_ongoing()
        t, n, it = update_counters(*self.read_checkpoint(t, n), it, dt)
        init_done = True

        while is_coupling_ongoing:
            self.write_checkpoint(t, n)
            dt = fem.Constant(self.mesh.domain, self.precice.get_max_time_step_size())

            self.coupling_read(dt)
            # if it == 0 and n == 1: self.write_output(t, n, "micro-out")
            self.problem.solve()
            self.write_output(t, n, it)

            is_coupling_ongoing = self.precice.is_coupling_ongoing()
            if is_coupling_ongoing:
                self.coupling_write()
                self.precice.advance(dt(0))
            is_coupling_ongoing = self.precice.is_coupling_ongoing()

            t, n, it = update_counters(*self.read_checkpoint(t, n), it, dt)

        # clean up
        CoupledSim._delete_adapter_file(self.adapter_config_path)

    @staticmethod
    def coupling_bc(x):
        """
        Define the coupling boundary (all points).

        Parameters
        ----------
        x : numpy.ndarray
            Array of coordinates.

        Returns
        -------
        numpy.ndarray
            Boolean array that is True for all points.
        """
        return Simulation.always_true(x)
