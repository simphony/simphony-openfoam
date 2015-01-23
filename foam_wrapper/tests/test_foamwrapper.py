""" test_foamwrapper module

This module contains the unitary tests for the
mesh module functionalities

"""

import unittest
from foam_wrapper.model import Model
from foam_wrapper.foam_wrapper import FoamWrapper
import logging
logger = logging.getLogger(__name__)

class FoamWrapperTestCase(unittest.TestCase):
    """Test case for FoamWrapper class"""
    def setUp(self):
        """Creates dummy model to perform tests"""
        self.model = Model(1)


    def test_meshRead(self):
        """Test mesh read from OpenFoam to Simphony"""
        foam_wrapper = FoamWrapper(self.model)
        simphonyMesh = foam_wrapper.get_mesh("foam_wrapper/tests/FOAM_Meshes/pitzDaily")
        self.assertEqual(sum(1 for _ in simphonyMesh.iter_points()), 25012)

        self.assertEqual(sum(1 for _ in simphonyMesh.iter_cells()), 12225)

        self.assertEqual(sum(1 for _ in simphonyMesh.iter_faces()), 49180)


if __name__ == '__main__':
    unittest.main()

