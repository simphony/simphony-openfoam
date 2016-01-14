""" foam_controlwrapper module

Wrapper module for OpenFOAM

"""
import os

from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_mesh import ABCMesh

from .cuba_extension import CUBAExt
from .foam_mesh import FoamMesh
from .foam_runner import FoamRunner
from .foam_files import modify_files, write_default_files, remove_parser_files
from .foam_templates import get_foam_solver
import simphonyfoaminterface


class Wrapper(ABCModelingEngine):
    """ Wrapper to OpenFOAM

    """

    def __init__(self):
        super(Wrapper, self).__init__()
        self._meshes = {}
        self.CM = DataContainer()
        self.BC = DataContainer()
        self.SP = DataContainer()
        #: to be able to use CUBAExt keywords, which are not in accepted
        #  CUBA keywords these extensions to CM and SP is used
        self.CM_extensions = {}
        self.SP_extensions = {}

    def run(self):
        """Run OpenFoam based on CM, BC and SP data

        Raises
        ------
        Exception when solver not supported.
        Exception when file IO occurs.

        """

        if not self.CM[CUBA.NAME]:
            error_str = "Mesh name must be defined in CM[CUBA.NAME]"
            raise ValueError(error_str)

        if not self._meshes:
            error_str = "Meshes not added to wrapper. Use add_mesh method"
            raise ValueError(error_str)

        name = self.CM[CUBA.NAME]

        if not self._meshes[name]:
            error_str = "Mesh {} does not exist"
            raise ValueError(error_str.format(name))

        mesh = self._meshes[name]
        case = mesh.path

        GE = self.CM_extensions[CUBAExt.GE]
        solver = get_foam_solver(self.CM_extensions)
        turbulent = 'Turbulent' if not (CUBAExt.LAMINAR_MODEL in GE) else ''

        # write default files based on solver
        # (not field data files if exists)
        templateName = solver + turbulent
        write_default_files(case, templateName, mesh._time, True)

        # write first mesh from foams objectRegistry to disk
        mesh.write()

        # modify control and boundary data files based on SP and BC
        modify_files(
            case, mesh._time, self.SP, self.BC,
            solver, self.SP_extensions, self.CM_extensions)

        # run case
        if CUBAExt.NUMBER_OF_CORES in self.CM_extensions:
            ncores = self.CM_extensions[CUBAExt.NUMBER_OF_CORES]
        else:
            ncores = 1
        runner = FoamRunner(solver, case, ncores)
        runner.run()

        # remove PyFoam parser files
        remove_parser_files(os.getcwd())

        # save timestep to mesh
        mesh._time = runner.get_last_time()

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
        Exception if mesh not instance of ABCMesh

        """

        if not isinstance(mesh, ABCMesh):
            raise TypeError('Mesh not instance of ABCMesh')

        mesh_name = mesh.name
        if name:
            mesh_name = name
        if mesh_name in self._meshes:
            raise ValueError('Mesh \'{}\` already exists'.format(mesh_name))
        else:
            if self.BC:
                self._meshes[mesh_name] = FoamMesh(mesh_name, self.BC, mesh)
            else:
                self._meshes[mesh_name] = FoamMesh(mesh_name, {}, mesh)
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
            simphonyfoaminterface.deleteMesh(name)
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
