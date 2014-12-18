""" test_foamwrapper module

This module contains the unitary tests for the
mesh module functionalities

"""

import unittest
from foam_wrapper.model import Model
from foam_wrapper.foam_wrapper import Foam_wrapper
import os

class FoamWrapperTestCase(unittest.TestCase):
    """Test case for Foam_wrapper class"""
    def setUp(self):
        """Creates dummy model to perform tests"""
        self.model = Model(1)


    def test_meshRead(self):
        """Test mesh read from OpenFoam to Simphony"""
        foam_wrapper = Foam_wrapper(self.model)
        curPath = os.getcwd()
        casePath = curPath
        simphonyMesh = foam_wrapper.get_mesh("foam_wrapper/tests/FOAM_Meshes/pitzDaily")
        print "Number of mesh points: ",sum(1 for _ in simphonyMesh.iter_points())
        print "Number of mesh cells : ",sum(1 for _ in simphonyMesh.iter_cells())
        print "Number of mesh faces : ",sum(1 for _ in simphonyMesh.iter_faces())

if __name__ == '__main__':
    unittest.main()
