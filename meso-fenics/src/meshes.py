from typing import Callable, Optional

import numpy as np

from mpi4py import MPI
from dolfinx import io, fem, mesh
import ufl

from .config import Config
from .mesh_utils import Locators, Dimensions

class Mesh:
    """
    Mesh class for loading and managing computational meshes.

    This class handles mesh loading from GMSH files, domain normalization,
    and boundary condition setup.

    Attributes
    ----------
    domain : dolfinx.mesh.Mesh
        The computational domain mesh.
    cell_markers : dolfinx.mesh.MeshTags
        Cell markers for subdomains.
    facet_markers : dolfinx.mesh.MeshTags
        Facet markers for boundaries.
    dimensions : Dimensions
        Original dimensions of the mesh before normalization.
    ds : ufl.Measure
        Surface measure for boundary integrals.
    bc_dc_use_tag : bool
        Whether to use tags for Dirichlet boundary conditions.
    bc_dc_locator : callable
        Locator function for Dirichlet BC.
    bc_dc_value : fem.Constant
        Value for Dirichlet BC.
    bc_nm_use_tag : bool
        Whether to use tags for Neumann boundary conditions.
    bc_nm_locator : callable
        Locator function for Neumann BC.
    bc_nm_value : fem.Constant
        Value for Neumann BC.
    bc_nm_dim : int
        Dimension component for Neumann BC.
    """
    def __init__(self):
        self.domain         = None
        self.cell_markers   = None
        self.facet_markers  = None
        self.dimensions     = None # original dimensions
        self.ds             = None

        self.bc_dc_use_tag  = False
        self.bc_dc_locator: Optional[Callable] = Locators.plane_yz_low
        self.bc_dc_value    = None
        self.bc_nm_use_tag  = False
        self.bc_nm_locator: Optional[Callable] = Locators.plane_yz_high
        self.bc_nm_value    = None
        self.bc_nm_dim      = 1

    def load(self, config: Config):
        """
        Load mesh from file and configure boundary conditions.

        Parameters
        ----------
        config : Config
            Configuration object containing mesh path and BC parameters.
        """
        mesh_data = io.gmsh.read_from_msh(config.mesh_path, MPI.COMM_WORLD, gdim=3)
        self.domain = mesh_data.mesh
        low, dims = self.normalize_domain()
        self.dimensions = Dimensions(low, dims, [0, 0, 0])

        self.bc_dc_use_tag = config.bc_dc_use_tag
        self.bc_dc_locator = Locators.from_name(config.bc_dc_locator)
        self.bc_nm_use_tag = config.bc_nm_use_tag
        self.bc_nm_locator = Locators.from_name(config.bc_nm_locator)
        self.bc_nm_dim     = config.bc_nm_dim

        # create nm tags if locator should be used
        if not self.bc_nm_use_tag:
            facet_markers = mesh_data.facet_tags
            facet_indices = facet_markers.indices.copy()
            facet_values = facet_markers.values.copy()
            com = mesh.compute_midpoints(mesh_data.mesh, 2, facet_markers.indices)
            facet_values[self.bc_nm_locator(np.swapaxes(com, 0, 1))] = Locators.NEUMANN_TAG
            self.facet_markers = mesh.meshtags(mesh_data.mesh, 2, facet_indices, facet_values)
        else:
            self.facet_markers = mesh_data.facet_tags
        self.cell_markers = mesh_data.cell_tags

        self.ds = ufl.Measure('ds', domain=self.domain, subdomain_data=self.facet_markers)

        self.bc_dc_value = fem.Constant(self.domain, config.bc_dc_value)
        self.bc_nm_value = fem.Constant(self.domain, config.bc_nm_value)

    def normalize_domain(self):
        """
        Normalize the domain to fit in a unit cube.

        Shifts the mesh so the lower corner is at the origin and scales
        it to have maximum dimension of 1.

        Returns
        -------
        low : numpy.ndarray
            Original lower corner coordinates before normalization.
        dims : float
            Original maximum dimension before normalization.
        """
        # find bounding box
        low = np.min(self.domain.geometry.x, axis=0)
        high = np.max(self.domain.geometry.x, axis=0)

        # shift lower corner to origin
        self.domain.geometry.x[:] -= low

        # rescale to range [0, 1]
        dims = np.max(np.abs(high - low))
        self.domain.geometry.x[:] /= dims

        return low, dims

    def get_bc_diriclet(self, V: fem.FunctionSpace):
        """
        Create Dirichlet boundary conditions for the problem.

        Parameters
        ----------
        V : fem.FunctionSpace
            Vector function space for displacement.

        Returns
        -------
        list of fem.DirichletBC
            List of Dirichlet boundary conditions for each component.
        """
        bcs = list()
        fdim = self.domain.topology.dim - 1

        if self.bc_dc_use_tag:
            for i in range(3):
                facets = self.facet_markers.find(Locators.DIRICLET_TAG)
                dofs   = fem.locate_dofs_topological(V.sub(i), fdim, facets)
                bcs.append(fem.dirichletbc(self.bc_dc_value, dofs, V.sub(i)))
        else:
            for i in range(3):
                facets = mesh.locate_entities_boundary(self.domain, fdim, self.bc_dc_locator)
                dofs = fem.locate_dofs_topological(V.sub(i), fdim, facets)
                bcs.append(fem.dirichletbc(self.bc_dc_value, dofs, V.sub(i)))

        return bcs

    def get_bc_neumann(self, v):
        """
        Create Neumann boundary condition term for the weak form.

        Parameters
        ----------
        v : ufl.Argument
            Test function.

        Returns
        -------
        ufl.core.expr.Expr
            Boundary integral term for Neumann BC.
        """
        return ufl.inner(self.bc_nm_value, v[self.bc_nm_dim]) * self.ds(Locators.NEUMANN_TAG)
