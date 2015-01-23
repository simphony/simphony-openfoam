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


    def test_internalScalarPointValues(self):
        """Test internalScalarPointValues read from OpenFoam"""
        foam_tokenwrapper = FoamTokenWrapper(self.model)
        argv=["test","-case","foam_tokenwrapper/tests/pitzDaily"]
        foam_tokenwrapper.initArgs(argv)
        foam_tokenwrapper.readMesh()
        
        values = foam_tokenwrapper.getInternalScalarPointValues("p")
        self.assertEqual(len(values), 25012)
        self.assertAlmostEqual(min(values), -8.15520920527) 
        self.assertAlmostEqual(max(values), 14.3523269266) 


if __name__ == '__main__':
    unittest.main()

