""" foam_internalwrapper module

Wrapper module for OpenFOAM

"""
from simphony.core.cuba import CUBA
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.meta import api

from .foam_mesh import FoamMesh
from .foam_dicts import (modifyNumerics, modifyFields)
from .foam_dicts import get_foam_solver
import simphonyfoaminterface as foamface


class Wrapper(ABCModelingEngine):
    """ Wrapper to OpenFOAM

    """

    def __init__(self, **kwargs):

        self._meshes = {}
        super(Wrapper, self).__init__(**kwargs)

    def _load_cuds(self):
        """Load CUDS data into  engine."""
        cuds = self.get_cuds()
        if not cuds:
            return

        for component in cuds.iter(ABCMesh):
            # while mesh can be from other mesh engine add_dataset creates
            # new FoamMesh. While there is no other way to know the
            # new mesh name on cuds the corresponding mesh component is
            # removed from cuds and replaced with the new mesh.
            # If the component is wanted to remain on cuds it must be
            # duplicated on user level code using some other name
            new_mesh = self.add_dataset(component)
            cuds.remove(component.name)
            cuds.add(new_mesh)

    def run(self):
        """ run OpenFoam based on CM, BC and SP data

        Returns
        -------
        lastTime : lastTime
            Last time step taken.

        Raises
        ------
        Exception when solver not supported.
        Exception when file IO occurs.

        """

        if not self._meshes:
            error_str = "Meshes not added to wrapper. Use add_mesh method"
            raise ValueError(error_str)

        if len(self._meshes) > 1:
            error_str = "Multiple meshes not supported"
            raise ValueError(error_str)

        cuds = self.get_cuds()
        solver = get_foam_solver(cuds)

        mesh = self._meshes.values()[0]

#       a) Modify fvSchemes and fvSolution
        modifyNumerics(mesh, cuds, solver)

#       b) Set boundary condition and Fields
        modifyFields(mesh, cuds, solver)

#       c) Call solver
        solver_parameters = cuds.iter(api.SolverParameter)
        ncores = 1
        if solver_parameters is not None:
            for sp in solver_parameters:
                if CUBA.NUMBER_OF_CORES in sp._data:
                    ncores = sp._data[CUBA.NUMBER_OF_CORES]
                    break
        mesh._time = foamface.run(mesh.name, ncores, solver)

    def add_dataset(self, mesh, name=None):
        """Add a mesh to the OpenFoam modeling engine.

        Parameters
        ----------
        mesh : ABCMesh
            mesh to be added.
        name : string
            name to give to mesh (optional)

        Returns
        -------
        proxy : FoamMesh
            A proxy mesh to be used to update/query the internal representation
            stored inside the modeling-engine. See get_mesh for more
            information.

        Raises
        ------
        Exception if mesh already exists

        """

        mesh_name = mesh.name
        if name:
            mesh_name = name
        if mesh_name in self._meshes:
            raise ValueError('Mesh \'{}\` already exists'.format(mesh_name))
        else:
            cuds = self.get_cuds()
            if cuds:
                solver = get_foam_solver(cuds)
            else:
                solver = 'pimpleFoam'
            self._meshes[mesh_name] = FoamMesh(mesh_name, cuds,
                                               solver, mesh)
            return self._meshes[mesh_name]

    def remove_dataset(self, name):
        """Delete mesh from the OpenFoam modeling engine.

        Parameters
        ----------
        name : str
            name of the mesh to be deleted.


        Raises
        ------
        Exception if mesh not found

        """

        if name not in self._meshes:
            raise ValueError('Mesh \'{}\` does not exists'.format(name))
        else:
            foamface.deleteMesh(name)
            del self._meshes[name]

    def get_dataset(self, name):
        """Get a mesh.

        The returned mesh can be used to query and update the state of the
        mesh inside the OpenFoam modeling engine.

        Parameters
        ----------
        name : str
            name of the mesh to be retrieved.

        Returns
        -------
        FoamMesh

        Raises
        ------
        Exception if mesh not found

        """

        if name in self._meshes:
            return self._meshes[name]
        else:
            raise ValueError(
                'Mesh \'{}\` does not exist'.format(name))

    def iter_datasets(self, names=None):
        """Returns an iterator over a subset or all of the meshes.

        Parameters
        ----------
        names : list of str
            names of specific meshes to be iterated over.
            If names is not given, then all meshes will
            be iterated over.

        Returns
        ----------
        Iterator over a subset or all of the meshes

        Raises
        ------
        Exception if some mesh fron mesh names list not found

        """

        if names is None:
            for name in self._meshes:
                yield self._meshes[name]
        else:
            for name in names:
                if name in self._meshes:
                    yield self._meshes[name]
                else:
                    raise ValueError(
                        'Mesh \'{}\` does not exist'.format(
                            name))

    def get_dataset_names(self):
        """ Returns the names of the meshes.

        """

        return self._meshes.keys()
