""" test_foamcontrolwrapper module

This module contains the unitary tests for the
foam_controlwrapper module functionalities

"""

import unittest
from simphony.cuds.mesh import Mesh, Face, Point, Cell
from simphony.core.data_container import DataContainer
from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from simphony.core.cuba import CUBA
import logging

logger = logging.getLogger(__name__)


class FoamControlWrapperTestCase(unittest.TestCase):
    """Test case for FoamControlWrapper class"""
    def setUp(self):

        self.mesh = Mesh(name="mesh1")

        self.points = [
            Point(
                (0.0, 0.0, 0.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]})),
            Point(
                (1.0, 0.0, 0.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]})),
            Point(
                (1.0, 1.0, 0.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]})),
            Point(
                (0.0, 1.0, 0.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]})),
            Point(
                (0.0, 0.0, 1.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]})),
            Point(
                (1.0, 0.0, 1.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]})),
            Point(
                (1.0, 1.0, 1.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]})),
            Point(
                (0.0, 1.0, 1.0),
                data=DataContainer({CUBA.VELOCITY: [0, 0, 0]}))
        ]

        puids = [self.mesh.add_point(point) for point in self.points]

        self.faces = [
            Face([puids[0], puids[3], puids[7], puids[4]],
                 data=DataContainer({CUBA.LABEL: 0})),
            Face([puids[1], puids[2], puids[6], puids[5]],
                 data=DataContainer({CUBA.LABEL: 1})),
            Face([puids[0], puids[1], puids[5], puids[4]],
                 data=DataContainer({CUBA.LABEL: 2})),
            Face([puids[3], puids[2], puids[6], puids[7]],
                 data=DataContainer({CUBA.LABEL: 3})),
            Face([puids[0], puids[1], puids[2], puids[3]],
                 data=DataContainer({CUBA.LABEL: 4})),
            Face([puids[4], puids[5], puids[6], puids[7]],
                 data=DataContainer({CUBA.LABEL: 5}))

        ]

        [self.mesh.add_face(face) for face in self.faces]

        self.cells = [
            Cell(puids)
        ]

        [self.mesh.add_cell(cell) for cell in self.cells]

    def test_add_mesh(self):
        """Test add_mesh method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_mesh(self.mesh)

    def test_get_mesh(self):
        """Test get_mesh method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_mesh(self.mesh)
        mesh_inside_wrapper = wrapper.get_mesh(self.mesh.name)
        self.assertEquals(self.mesh.name, mesh_inside_wrapper.name)

    def test_multiple_meshes(self):
        """Test multiple meshes inside wrapper

        """

        wrapper = FoamControlWrapper()
        wrapper.add_mesh(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_mesh(mesh2)
        mesh_inside_wrapper1 = wrapper.get_mesh(self.mesh.name)
        mesh_inside_wrapper2 = wrapper.get_mesh(mesh2.name)

        assert(sum(1 for _ in mesh_inside_wrapper1.iter_points()) ==
               sum(1 for _ in mesh_inside_wrapper2.iter_points()))

    def test_iter_meshes(self):
        """Test iter_meshes method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_mesh(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_mesh(mesh2)

        assert(sum(1 for _ in wrapper.iter_meshes()) == 2)

    def test_delete_mesh(self):
        """Test delete_mesh method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_mesh(self.mesh)
        wrapper.delete_mesh(self.mesh.name)

if __name__ == '__main__':
    unittest.main()
