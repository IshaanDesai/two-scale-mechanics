import json
from collections import defaultdict

from .mesh_utils import Locators


def get_else(val, default):
    """
    Return the value if not None, otherwise return the default.

    Parameters
    ----------
    val : any
        The value to check.
    default : any
        The default value to return if val is None.

    Returns
    -------
    any
        The value or default.
    """
    if val is None:
        return default
    return val


def get_raise(val, msg=""):
    """
    Return the value if not None, otherwise raise a RuntimeError.

    Parameters
    ----------
    val : any
        The value to check.
    msg : str, optional
        Error message to include in the exception (default is "").

    Returns
    -------
    any
        The value if not None.

    Raises
    ------
    RuntimeError
        If val is None.
    """
    if val is None:
        raise RuntimeError(msg)
    return val


def make_default_dict(d=None):
    """
    Create a defaultdict that returns None for missing keys.

    Parameters
    ----------
    d : dict, optional
        Initial dictionary to populate the defaultdict (default is None).

    Returns
    -------
    defaultdict
        A defaultdict with None as the default value.
    """
    if d is None:
        return defaultdict(lambda: None)
    else:
        return defaultdict(lambda: None, d)


class Config:
    """
    Configuration class for loading and storing simulation parameters.

    Attributes
    ----------
    output_path : str
        Path for output files.
    mesh_path : str
        Path to the mesh file.
    bc_dc_use_tag : bool
        Whether to use tags for Dirichlet boundary conditions.
    bc_dc_locator : str
        Locator function name for Dirichlet BC.
    bc_dc_value : float
        Value for Dirichlet BC.
    bc_nm_use_tag : bool
        Whether to use tags for Neumann boundary conditions.
    bc_nm_locator : str
        Locator function name for Neumann BC.
    bc_nm_value : float
        Value for Neumann BC.
    bc_nm_dim : int
        Dimension for Neumann BC.
    problem_lambda : float
        Lamé parameter lambda.
    problem_mu : float
        Lamé parameter mu.
    problem_alpha : float
        Material nonlinearity parameter.
    problem_strain_type : str
        Type of strain formulation ('small_strain' or 'large_strain').
    problem_element_degree : int
        Finite element polynomial degree.
    simulation_type : str
        Type of simulation to run.
    simulation_input : str
        Path to simulation input file.
    simulation_micro_type : str
        Type of micro-solver.
    simulation_write_state : str
        Path to write state data.
    simulation_write_state_type : list
        Types of state data to write.
    simulation_init_tan : bool
        Whether to initialize the simulation using tangent.
    simulation_prc_path : str
        Path to preCICE config file
    simulation_slurm_id : str
        Slurm ID for simulation. (job id)
    """

    def __init__(self):
        self._data = None

        self._mesh_path = None

        self._bc_dc_use_tag = None
        self._bc_dc_locator = None
        self._bc_dc_value = None
        self._bc_nm_use_tag = None
        self._bc_nm_locator = None
        self._bc_nm_value = None
        self._bc_nm_dim = None

        self._problem_lambda = None
        self._problem_mu = None
        self._problem_alpha = None
        self._problem_strain_type = None
        self._problem_elem_degree = None

        self._simulation_type = None
        self._simulation_input = None
        self._simulation_micro_t = None
        self._simulation_wstate = None
        self._simulation_wstate_t = None
        self._simulation_init_tan = None
        self._simulation_prc_path = None
        self._simulation_slurm_id = None

        self._output_path = None

    def load(self, path: str):
        """
        Load configuration from a JSON file.

        Parameters
        ----------
        path : str
            Path to the JSON configuration file.

        Raises
        ------
        RuntimeError
            If required configuration fields are missing or invalid.
        """
        with open(path, "r") as f:
            self._data = json.load(f)
        self._data = make_default_dict(self._data)

        # general
        self._output_path = get_else(self._data["output_path"], "mf_out")

        # MESH
        if self._data["mesh"] is None:
            raise RuntimeError("Must provide mesh descriptor in config.")
        self._data["mesh"] = make_default_dict(self._data["mesh"])
        self._mesh_path = self._data["mesh"]["path"]
        self._bc_dc_use_tag = get_raise(
            self._data["mesh"]["bc_dc_use_tag"], "Mesh requires bc_dc_use_tag field."
        )
        if not self._bc_dc_use_tag:
            self._bc_dc_locator = get_raise(
                self._data["mesh"]["bc_dc_locator"],
                "Mesh either requires tags or locator.",
            )
            if not Locators.is_name_valid(self._bc_dc_locator):
                raise RuntimeError("Invalid bc_dc_locator name.")
        self._bc_dc_value = get_else(self._data["mesh"]["bc_dc_value"], 0.0)
        self._bc_nm_use_tag = get_raise(
            self._data["mesh"]["bc_nm_use_tag"], "Mesh requires bc_nm_use_tag field."
        )
        if not self._bc_nm_use_tag:
            self._bc_nm_locator = get_raise(
                self._data["mesh"]["bc_nm_locator"],
                "Mesh either requires tags or locator.",
            )
            if not Locators.is_name_valid(self._bc_nm_locator):
                raise RuntimeError("Invalid bc_nm_locator name.")
        self._bc_nm_value = get_else(self._data["mesh"]["bc_nm_value"], 0.0)
        self._bc_nm_dim = get_else(self._data["mesh"]["bc_nm_dim"], 1)

        # PROBLEM
        if self._data["problem"] is not None:
            self._data["problem"] = make_default_dict(self._data["problem"])
            self._problem_lambda = get_else(self._data["problem"]["lambda"], 10.0)
            self._problem_mu = get_else(self._data["problem"]["mu"], 5.0)
            self._problem_alpha = get_else(self._data["problem"]["alpha"], 100.0)
            self._problem_strain_type = get_else(
                self._data["problem"]["strain_type"], "small_strain"
            )
            self._problem_elem_degree = get_else(
                self._data["problem"]["elem_degree"], 2
            )
        else:
            self._problem_lambda = 10.0
            self._problem_mu = 5.0
            self._problem_alpha = 100.0
            self._problem_strain_type = "small_strain"
            self._problem_elem_degree = 2

        # SIMULATION
        if self._data["simulation"] is None:
            raise RuntimeError("Must provide simulation descriptor.")
        self._data["simulation"] = make_default_dict(self._data["simulation"])
        self._simulation_type = get_raise(
            self._data["simulation"]["type"], "Simulation type required."
        )
        self._simulation_input = get_else(self._data["simulation"]["input"], None)
        self._simulation_micro_t = get_raise(
            self._data["simulation"]["micro_type"], "Simulation micro type required."
        )
        self._simulation_wstate = get_else(
            self._data["simulation"]["write_state"], None
        )
        self._simulation_wstate_t = get_else(
            self._data["simulation"]["write_state_type"], None
        )
        self._simulation_init_tan = get_else(
            self._data["simulation"]["init_with_micro"], False
        )
        self._simulation_prc_path = get_else(
            self._data["simulation"]["precice_xml_path"], None
        )
        self._simulation_slurm_id = get_else(
            self._data["simulation"]["slurm_id"], "default"
        )

    @property
    def output_path(self):
        return self._output_path

    @property
    def mesh_path(self):
        return self._mesh_path

    @property
    def bc_dc_use_tag(self):
        return self._bc_dc_use_tag

    @property
    def bc_dc_locator(self):
        return self._bc_dc_locator

    @property
    def bc_dc_value(self):
        return self._bc_dc_value

    @property
    def bc_nm_use_tag(self):
        return self._bc_nm_use_tag

    @property
    def bc_nm_locator(self):
        return self._bc_nm_locator

    @property
    def bc_nm_value(self):
        return self._bc_nm_value

    @property
    def bc_nm_dim(self):
        return self._bc_nm_dim

    @property
    def problem_lambda(self):
        return self._problem_lambda

    @property
    def problem_mu(self):
        return self._problem_mu

    @property
    def problem_alpha(self):
        return self._problem_alpha

    @property
    def problem_strain_type(self):
        return self._problem_strain_type

    @property
    def problem_element_degree(self):
        return self._problem_elem_degree

    @property
    def simulation_type(self):
        return self._simulation_type

    @property
    def simulation_input(self):
        return self._simulation_input

    @property
    def simulation_micro_type(self):
        return self._simulation_micro_t

    @property
    def simulation_write_state(self):
        return self._simulation_wstate

    @property
    def simulation_write_state_type(self):
        return self._simulation_wstate_t

    @property
    def simulation_init_with_micro(self):
        return self._simulation_init_tan

    @property
    def simulation_precice_xml_path(self):
        return self._simulation_prc_path

    @property
    def simulation_slurm_id(self):
        return self._simulation_slurm_id
