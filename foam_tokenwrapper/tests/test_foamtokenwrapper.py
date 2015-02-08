""" test_foamtokenwrapper module

This module contains the unitary tests for the
foam_tokenwrapper module functionalities

"""

import unittest
from simphony.core.cuba import CUBA
from foam_tokenwrapper.foam_tokenwrapper import FoamTokenWrapper
import logging
import os
import shutil
from PyFoam.Execution.UtilityRunner import UtilityRunner
logger = logging.getLogger(__name__)


class FoamTokenWrapperTestCase(unittest.TestCase):
    """Test case for FoamTokenWrapper class"""
    def setUp(self):
        pass

    def test_internalScalarPointValues(self):
        """Test mesh and point data read from OpenFoam"""
        caseName = os.path.join('foam_tokenwrapper', 'tests', 'pitzDaily')
        try:
            src = os.path.expandvars('$FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily')
            shutil.copytree(src, caseName)
        except shutil.Error as e:
            raise shutil.Error('Case directory not copied. Error: %s' % e)
        except OSError as e:
            raise OSError('Case directory not copied. Error: %s' % e)

        blockMesh = UtilityRunner(argv=["blockMesh", "-case", caseName], silent=True, logname="blockMesh")
        blockMesh.start()

        foam_tokenwrapper = FoamTokenWrapper()
        timeStep = "0"
        (simPhoNyMesh, foamMeshMap) = foam_tokenwrapper.get_mesh(caseName, timeStep)
        foam_tokenwrapper.get_point_data(simPhoNyMesh, foamMeshMap, CUBA.PRESSURE)

        assert(sum(1 for _ in simPhoNyMesh.iter_points()) == 25012)
        assert(sum(1 for _ in simPhoNyMesh.iter_cells()) == 12225)
        assert(sum(1 for _ in simPhoNyMesh.iter_faces()) == 49180)

if __name__ == '__main__':
    unittest.main()
