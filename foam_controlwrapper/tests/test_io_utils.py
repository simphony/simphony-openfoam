""" test_io_utils module

This module contains the unitary tests for the
io_utils module functionalities

"""

import unittest
import os
import tempfile

from foam_controlwrapper.foam_controlwrapper import Wrapper
from foam_controlwrapper.blockmesh_utils import create_quad_mesh
from foam_controlwrapper.io_utils import read_foammesh
from simphony.core.cuds_item import CUDSItem


class IOUtilsTestCase(unittest.TestCase):
    """Test case for io_utils"""
    def setUp(self):
        wrapper = Wrapper()
        self.path = os.path.join(tempfile.mkdtemp())
        self.name = "test_mesh"
        corner_points = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
                         (5.0, 5.0, 0.0), (0.0, 5.0, 0.0),
                         (0.0, 0.0, 1.0), (5.0, 0.0, 1.0),
                         (5.0, 5.0, 1.0), (0.0, 5.0, 1.0)]
        create_quad_mesh(self.path, self.name, wrapper, corner_points, 5, 5, 5)

        self.mesh = wrapper.get_dataset(self.name)

    def test_read_foammesh(self):
        """Test read_foammesh method

        """

        mesh_from_file = read_foammesh(self.name, self.path)

        item_type = CUDSItem.POINT
        self.assertEqual(mesh_from_file.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUDSItem.FACE
        self.assertEqual(mesh_from_file.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUDSItem.CELL
        self.assertEqual(mesh_from_file.count_of(item_type),
                         self.mesh.count_of(item_type))


if __name__ == '__main__':
    unittest.main()