""" test_foamtokenwrapper module

This module contains the unitary tests for the
foam_tokenwrapper module functionalities

"""

import unittest
from simphony.core.cuba import CUBA
from foam_tokenwrapper.foam_tokenwrapper import FoamTokenWrapper
import logging
logger = logging.getLogger(__name__)


class FoamTokenWrapperTestCase(unittest.TestCase):
    """Test case for FoamTokenWrapper class"""
    def setUp(self):
        pass

    def test_internalScalarPointValues(self):
        """Test mesh and point data read from OpenFoam"""
        foam_tokenwrapper = FoamTokenWrapper()
        caseName = os.path.join('foam_tokenwrapper','tests','pitzDaily')
        timeStep = "760"
        (simPhoNyMesh,foamMeshMap) = foam_tokenwrapper.get_mesh(caseName, timeStep)
        
        foam_tokenwrapper.get_point_data(simPhoNyMesh, foamMeshMap, CUBA.PRESSURE)

        minValue = 1.0e10
        maxValue = -1.0e10
        for point in simPhoNyMesh.iter_points():
            minValue = min(minValue,point.data.get(CUBA.PRESSURE))
            maxValue = max(maxValue,point.data.get(CUBA.PRESSURE))
 
        self.assertAlmostEqual(minValue, -8.15526920527) 
        self.assertAlmostEqual(maxValue, 14.3523315923) 


if __name__ == '__main__':
    unittest.main()

