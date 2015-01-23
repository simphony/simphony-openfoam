""" test_foamtokenwrapper module

This module contains the unitary tests for the
foam_tokenwrapper module functionalities

"""

import unittest
from foam_wrapper.model import Model
from foam_tokenwrapper.foam_tokenwrapper import FoamTokenWrapper
import logging
logger = logging.getLogger(__name__)


class FoamTokenWrapperTestCase(unittest.TestCase):
    """Test case for FoamTokenWrapper class"""
    def setUp(self):
        """Creates dummy model to perform tests"""
        self.model = Model(1)
        

    def test_internalVectorCellValues(self):
        """Test internalVectorCellValues read from OpenFoam"""
        foam_tokenwrapper = FoamTokenWrapper(self.model)
        argv=["simpleFoam","-case","foam_tokenwrapper/tests/pitzDaily"]
        foam_tokenwrapper.initArgs(argv)
        foam_tokenwrapper.readMesh()
        
        values = foam_tokenwrapper.getInternalVectorCellValues("U")

        minVal = 1.0e10
        maxVal = -1.0e10
        ir = len(values)-2
        for i in range(0,ir,3):
            minVal = min(minVal,values[i])
            maxVal = max(maxVal,values[i])

        self.assertAlmostEqual(min(values), -2.59356) 
        self.assertAlmostEqual(max(values), 10.1793) 
         

if __name__ == '__main__':
    unittest.main()

