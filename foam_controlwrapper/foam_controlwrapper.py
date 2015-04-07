""" foam_controlwrapper module

Wrapper module for OpenFOAM

"""
import os

from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.abc_modeling_engine import ABCModelingEngine

from .cuba_extension import CUBAExt
from .foam_mesh import FoamMesh
from .foam_runner import FoamRunner
from .foam_files import modify_files, write_default_files, remove_parser_files
import simphonyfoaminterface


class FoamControlWrapper(ABCModelingEngine):
    """ Wrapper to OpenFOAM

    """

    def __init__(self):

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

        Returns
        -------
        lastTime : lastTime
            Last time step taken.

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
        solver = "simpleFoam"
        if CUBAExt.LAMINAR_MODEL in GE:
            if CUBAExt.VOF in GE:
                solver = "interFoam"
            else:
                solver = "simpleFoam"
        else:
            error_str = "GE does not define supported solver: GE = {}"
            raise NotImplementedError(error_str.format(GE))

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
        return mesh._time

    def add_mesh(self, mesh):
        """Add a mesh to the OpenFoam modeling engine.

        Parameters
        ----------
        mesh : ABCMesh
            mesh to be added.

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

        if mesh.name in self._meshes:
            raise ValueError('Mesh \'{}\` already exists'.format(mesh.name))
        else:
            self._meshes[mesh.name] = FoamMesh(mesh.name, mesh)
            return self._meshes[mesh.name]

    def delete_mesh(self, name):
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

    def get_mesh(self, name):
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

    def iter_meshes(self, names=None):
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

    def add_particles(self, particle_container):
        message = 'FoamWrapper does not handle particle container'
        raise NotImplementedError(message)

    def get_particles(self, name):
        message = 'FoamWrapper does not handle particle container'
        raise NotImplementedError(message)

    def delete_particles(self, name):
        message = 'FoamWrapper does not handle particle container'
        raise NotImplementedError(message)

    def iter_particles(self, names=None):
        message = 'FoamWrapper does not handle particle container'
        raise NotImplementedError(message)

    def add_lattice(self, lattice):
        message = 'FoamWrapper does not handle lattice'
        raise NotImplementedError(message)

    def get_lattice(self, name):
        message = 'FoamWrapper does not handle lattice'
        raise NotImplementedError(message)

    def delete_lattice(self, name):
        message = 'FoamWrapper does not handle lattice'
        raise NotImplementedError(message)

    def iter_lattices(self, names=None):
        message = 'FoamWrapper does not handle lattice'
        raise NotImplementedError(message)

def read_foammesh(name, path):
    """Read mesh from OpenFoam case files.

    Parameters
    ----------
    name : str
    name to give to mesh
    path : str
    case directory

    Raises
    ------
    Exception if some mesh from mesh names list not found

    """

    simphonyfoaminterface.init(name, path)
    simphonyfoaminterface.readMesh(name)
    nPoints = simphonyfoaminterface.getPointCount(name)
    nCells = simphonyfoaminterface.getCellCount(name)
    nFaces = simphonyfoaminterface.getFaceCount(name)
    nEdges = 0

    foamMesh = FoamMesh(name)
    foamMesh.generate_uuidmapping(nPoints, nEdges, nFaces, nCells)
    return foamMesh
