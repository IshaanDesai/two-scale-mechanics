import json
from collections import defaultdict

from .mesh_utils import Locators

def get_else(val, default):
    if val is None: return default
    return val

def get_raise(val, msg=""):
    if val is None: raise RuntimeError(msg)
    return val

def make_default_dict(d=None):
    if d is None:
        return defaultdict(lambda:None)
    else:
        return defaultdict(lambda:None, d)

class Config:
    def __init__(self):
        self._data                = None

        self._mesh_path           = None

        self._bc_dc_use_tag       = None
        self._bc_dc_locator       = None
        self._bc_dc_value         = None
        self._bc_nm_use_tag       = None
        self._bc_nm_locator       = None
        self._bc_nm_value         = None
        self._bc_nm_dim           = None

        self._problem_lambda      = None
        self._problem_mu          = None
        self._problem_alpha       = None
        self._problem_strain_type = None
        self._problem_elem_degree = None

        self._simulation_type     = None
        self._simulation_input    = None
        self._simulation_micro_t  = None
        self._simulation_wstate   = None
        self._simulation_wstate_t = None

        self._output_path         = None

    def load(self, path: str):
        with open(path, 'r') as f:
            self._data = json.load(f)
        self._data = make_default_dict(self._data)

        # general
        self._output_path             = get_else(self._data['output_path'], 'mf_out')

        # MESH
        if self._data['mesh'] is None:
            raise RuntimeError('Must provide mesh descriptor in config.')
        self._data['mesh'] = make_default_dict(self._data['mesh'])
        self._mesh_path               = self._data['mesh']['path']
        self._bc_dc_use_tag           = get_raise(self._data['mesh']['bc_dc_use_tag'], 'Mesh requires bc_dc_use_tag field.')
        if not self._bc_dc_use_tag:
            self._bc_dc_locator       = get_raise(self._data['mesh']['bc_dc_locator'], 'Mesh either requires tags or locator.')
            if not Locators.is_name_valid(self._bc_dc_locator):
                raise RuntimeError('Invalid bc_dc_locator name.')
        self._bc_dc_value             = get_else(self._data['mesh']['bc_dc_value'], 0.0)
        self._bc_nm_use_tag           = get_raise(self._data['mesh']['bc_nm_use_tag'], 'Mesh requires bc_nm_use_tag field.')
        if not self._bc_nm_use_tag:
            self._bc_nm_locator       = get_raise(self._data['mesh']['bc_nm_locator'], 'Mesh either requires tags or locator.')
            if not Locators.is_name_valid(self._bc_nm_locator):
                raise RuntimeError('Invalid bc_nm_locator name.')
        self._bc_nm_value             = get_else(self._data['mesh']['bc_nm_value'], 0.0)
        self._bc_nm_dim               = get_else(self._data['mesh']['bc_nm_dim'], 1)

        # PROBLEM
        if self._data['problem'] is not None:
            self._data['problem'] = make_default_dict(self._data['problem'])
            self._problem_lambda      = get_else(self._data['problem']['lambda'], 10.0)
            self._problem_mu          = get_else(self._data['problem']['mu'], 5.0)
            self._problem_alpha       = get_else(self._data['problem']['alpha'], 100.0)
            self._problem_strain_type = get_else(self._data['problem']['strain_type'], 'small_strain')
            self._problem_elem_degree = get_else(self._data['problem']['elem_degree'], 2)
        else:
            self._problem_lambda      = 10.0
            self._problem_mu          = 5.0
            self._problem_alpha       = 100.0
            self._problem_strain_type = 'small_strain'
            self._problem_elem_degree = 2

        # SIMULATION
        if self._data['simulation'] is None:
            raise RuntimeError('Must provide simulation descriptor.')
        self._data['simulation']  = make_default_dict(self._data['simulation'])
        self._simulation_type     = get_raise(self._data['simulation']['type'], 'Simulation type required.')
        self._simulation_input    = get_else(self._data['simulation']['input'], None)
        self._simulation_micro_t  = get_raise(self._data['simulation']['micro_type'], 'Simulation micro type required.')
        self._simulation_wstate   = get_else(self._data['simulation']['write_state'], None)
        self._simulation_wstate_t = get_else(self._data['simulation']['write_state_type'], None)


    @property
    def output_path(self): return self._output_path
    @property
    def mesh_path(self): return self._mesh_path
    @property
    def bc_dc_use_tag(self): return self._bc_dc_use_tag
    @property
    def bc_dc_locator(self): return self._bc_dc_locator
    @property
    def bc_dc_value(self): return self._bc_dc_value
    @property
    def bc_nm_use_tag(self): return self._bc_nm_use_tag
    @property
    def bc_nm_locator(self): return self._bc_nm_locator
    @property
    def bc_nm_value(self): return self._bc_nm_value
    @property
    def bc_nm_dim(self): return self._bc_nm_dim
    @property
    def problem_lambda(self): return self._problem_lambda
    @property
    def problem_mu(self): return self._problem_mu
    @property
    def problem_alpha(self): return self._problem_alpha
    @property
    def problem_strain_type(self): return self._problem_strain_type
    @property
    def problem_element_degree(self): return self._problem_elem_degree
    @property
    def simulation_type(self): return self._simulation_type
    @property
    def simulation_input(self): return self._simulation_input
    @property
    def simulation_micro_type(self): return self._simulation_micro_t
    @property
    def simulation_write_state(self): return self._simulation_wstate
    @property
    def simulation_write_state_type(self): return self._simulation_wstate_t
