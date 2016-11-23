""" test_io_utils module

This module contains the unitary tests for the
io_utils module functionalities

"""

import unittest
import os
import tempfile

from simphony.core.cuba import CUBA

from foam_controlwrapper.blockmesh_utils import create_quad_mesh
from foam_controlwrapper.io_utils import read_foammesh


class IOUtilsTestCase(unittest.TestCase):
    """Test case for io_utils"""
    def setUp(self):
        self.path = os.path.join(tempfile.mkdtemp())
        self.name = "test_mesh"
        corner_points = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
                         (5.0, 5.0, 0.0), (0.0, 5.0, 0.0),
                         (0.0, 0.0, 1.0), (5.0, 0.0, 1.0),
                         (5.0, 5.0, 1.0), (0.0, 5.0, 1.0)]
        self.mesh = create_quad_mesh(self.path, self.name, corner_points,
                                     5, 5, 5)


    def test_read_foammesh(self):
        """Test read_foammesh method

        """

        mesh_from_file = read_foammesh(self.name, self.path)

        item_type = CUBA.POINT
        self.assertEqual(mesh_from_file.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUBA.FACE
        self.assertEqual(mesh_from_file.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUBA.CELL
        self.assertEqual(mesh_from_file.count_of(item_type),
                         self.mesh.count_of(item_type))


if __name__ == '__main__':
    unittest.main()
