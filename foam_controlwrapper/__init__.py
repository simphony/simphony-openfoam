# Functions, classes and constants exported here will be available
# when the `openfoam` module is imported.

from simphony.engine import ABCEngineExtension
from simphony.engine import EngineInterface
from simphony.engine.decorators import register

from foam_internalwrapper.foam_internalwrapper import (
    Wrapper as InternalWrapper)
from foam_controlwrapper import Wrapper as FileIOWrapper
from .io_utils import read_foammesh
from .blockmesh_utils import (create_quad_mesh, create_block_mesh)

__all__ = ['read_foammesh', 'create_quad_mesh',
           'create_block_mesh', 'EngineType']


@register
class OpenFOAMExtension(ABCEngineExtension):
    """Simphony-OpenFOAM extension.
    This extension provides support for OpenFoam engines.
    """

    def get_supported_engines(self):
        """Get metadata about supported engines.
        Returns
        -------
        list: a list of EngineMetadata objects
        """
        foam_features = None
        foam = self.create_engine_metadata('OpenFOAM',
                                           foam_features,
                                           [EngineInterface.Internal,
                                            EngineInterface.FileIO])

        return [foam]

    def create_wrapper(self, cuds, engine_name, engine_interface):
        """Creates a wrapper to the requested engine.
        Parameters
        ----------
        cuds: CUDS
          CUDS computational model data
        engine_name: str
          name of the engine, must be supported by this extension
        engine_interface: EngineInterface
          the interface to interact with engine
        Returns
        -------
        ABCEngineExtension: A wrapper configured with cuds and ready to run
        """
        if engine_name != 'OpenFOAM':
            raise Exception('Only OpenFOAM engine is supported. '
                            'Unsupported engine: %s', engine_name)

        if engine_interface == EngineInterface.Internal:
            return InternalWrapper(cuds=cuds)

        if engine_interface == EngineInterface.FileIO:
            return FileIOWrapper(cuds=cuds)
