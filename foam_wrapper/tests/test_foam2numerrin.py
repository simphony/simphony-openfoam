""" test_foamwrapper module

This module contains the unitary tests for the
mesh module functionalities

"""

import unittest
from foam_wrapper.model import Model
from foam_wrapper.foam_wrapper import FoamWrapper
from numerrin_wrapper.numerrin_wrapper import NumerrinWrapper
import logging
logger = logging.getLogger(__name__)

class FoamNumerrinWrapperTestCase(unittest.TestCase):
    """Test case for FoamWapper and NumerrinWrapper class"""
    def setUp(self):
        """Creates dummy model to perform tests"""
        self.model = Model(1)

        
    def test_foamMesh2numerrinMesh(self):
        """Test mesh read from OpenFoam to SimPhony, pass to Numerrin to save to CGNS -format and back to SimPhony mesh"""
        foam_wrapper = FoamWrapper(self.model)
        print "Read mesh from OpenFOAM"
        simphonyMesh = foam_wrapper.get_mesh("foam_wrapper/tests/FOAM_Meshes/pitzDaily")
        self.assertEqual(sum(1 for _ in simphonyMesh.iter_points()), 25012)

        self.assertEqual(sum(1 for _ in simphonyMesh.iter_cells()), 12225)

        self.assertEqual(sum(1 for _ in simphonyMesh.iter_faces()), 49180)

        numerrin_wrapper = NumerrinWrapper(self.model)
        numerrin_wrapper.importMesh(simphonyMesh)
        numerrin_wrapper.setProgramString("WriteCGNS(\"mesh.cgns\") mesh\n")
        numerrin_wrapper.run()
        mesh2=numerrin_wrapper.exportMesh()
        self.assertEqual(sum(1 for _ in mesh2.iter_points()), 25012)

        self.assertEqual(sum(1 for _ in mesh2.iter_cells()), 12225)

        self.assertEqual(sum(1 for _ in mesh2.iter_faces()), 49180)
        
if __name__ == '__main__':
    unittest.main()
