""" test_foamwrapper module

This module contains the unitary tests for the
mesh module functionalities

"""

import unittest
from foam_wrapper.model import Model
from foam_wrapper.foam_wrapper import Foam_wrapper
from numerrin_wrapper.numerrin_wrapper import NumerrinWrapper

class FoamNumerrinWrapperTestCase(unittest.TestCase):
    """Test case for Foam_wrapper and NumerrinWrapper class"""
    def setUp(self):
        """Creates dummy model to perform tests"""
        self.model = Model(1)

        
    def test_foamMesh2numerrinMesh(self):
        """Test mesh read from OpenFoam to SimPhony, pass to Numerrin to save to CGNS -format and back to SimPhony mesh"""
        foam_wrapper = Foam_wrapper(self.model)
        print "Read mesh from OpenFOAM"
        simphonyMesh = foam_wrapper.get_mesh("FOAM_Meshes/pitzDaily")
        print "Number of mesh points: ",sum(1 for _ in simphonyMesh.iter_points())
        print "Number of mesh cells : ",sum(1 for _ in simphonyMesh.iter_cells())
        print "Number of mesh faces : ",sum(1 for _ in simphonyMesh.iter_faces())

        print "Mesh to Numerrin"
        numerrin_wrapper = NumerrinWrapper(self.model)
        numerrin_wrapper.importMesh(simphonyMesh)
        print "Save to CGNS -format"
        numerrin_wrapper.setProgramString("WriteCGNS(\"mesh.cgns\") mesh\n")
        numerrin_wrapper.run()
        print "Mesh from Numerrin"
        mesh2=numerrin_wrapper.exportMesh()
        print type(mesh2).__name__
        print "Number of mesh points: ",sum(1 for _ in mesh2.iter_points())
        print "Number of mesh cells : ",sum(1 for _ in mesh2.iter_cells())
        print "Number of mesh faces : ",sum(1 for _ in mesh2.iter_faces())
        
if __name__ == '__main__':
    unittest.main()
