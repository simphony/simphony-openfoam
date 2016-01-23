""" test_foam_mesh module

This module contains the unitary tests for the
foam_mesh module functionalities

"""

import unittest

from simphony.cuds.mesh import Mesh, Face, Point, Cell, Edge
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.core.cuds_item import CUDSItem

from foam_internalwrapper.foam_mesh import FoamMesh


class FoamMeshTestCase(unittest.TestCase):
    """Test case for FoamMesh class"""
    def setUp(self):
        self.mesh = Mesh(name='test_mesh')
        self.solver = 'pimpleFoam'
        self.points = [
            Point(
                (0.0, 0.0, 0.0)),
            Point(
                (1.0, 0.0, 0.0)),
            Point(
                (1.0, 1.0, 0.0)),
            Point(
                (0.0, 1.0, 0.0)),
            Point(
                (0.0, 0.0, 1.0)),
            Point(
                (1.0, 0.0, 1.0)),
            Point(
                (1.0, 1.0, 1.0)),
            Point(
                (0.0, 1.0, 1.0))
        ]

        puids = self.mesh.add_points(self.points)

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

        self.edges = [Edge([puids[0], puids[3]])]

        self.mesh.add_faces(self.faces)

        self.cells = [
            Cell(puids,
                 data=DataContainer({CUBA.VELOCITY: [1, 0, 0],
                                     CUBA.PRESSURE: 4.0}))
        ]

        self.puids = puids

        self.mesh.add_cells(self.cells)

    def test_get_edge(self):
        """Test get_edge method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.get_edge(self.edges[0].uid)

    def test_get_face(self):
        """Test get_face method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        face = self.faces[4]
        face_f = foam_mesh.get_face(face.uid)
        self.assertEqual(face.points, face_f.points)

    def test_get_cell(self):
        """Test get_cell method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        cell = self.cells[0]
        cell_f = foam_mesh.get_cell(cell.uid)
        self.assertEqual(set(cell.points), set(cell_f.points))
        self.assertEqual(cell.data[CUBA.PRESSURE], cell_f.data[CUBA.PRESSURE])
        self.assertEqual(cell.data[CUBA.VELOCITY], cell_f.data[CUBA.VELOCITY])

    def test_add_points(self):
        """Test add_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add_points(self.points)

    def test_add_edges(self):
        """Test add_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add_edges(self.edges)

    def test_add_faces(self):
        """Test add_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add_faces(self.faces)

    def test_add_cells(self):
        """Test add_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add_cells(self.cells)

    def test_update_points(self):
        """Test update_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.update_points(self.points)

    def test_update_edges(self):
        """Test update_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.update_edges(self.edges)

    def test_update_faces(self):
        """Test update_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.update_faces(self.faces)

    def test_update_cells(self):
        """Test update_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        cell = self.cells[0]
        cell.data[CUBA.VELOCITY] = [2, 1, 3]
        foam_mesh.update_cells(self.cells)
        cell_f = foam_mesh.get_cell(cell.uid)
        self.assertIsInstance(cell_f.data, DataContainer)
        self.assertEqual(cell.data, cell_f.data)

        cell.points = [self.points[1].uid, self.points[2].uid,
                       self.points[3].uid]
        with self.assertRaises(Warning):
            foam_mesh.update_cells(self.cells)

    def test_iter_points(self):
        """Test iter_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        puids_f = [point.uid for point in foam_mesh.iter_points()]
        self.assertEqual(set(self.puids), set(puids_f))

    def test_iter_edges(self):
        """Test iter_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertEqual(foam_mesh.iter_edges(), [])

    def test_iter_faces(self):
        """Test iter_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for face_f in foam_mesh.iter_faces():
            face = self.mesh.get_face(face_f.uid)
            self.assertEqual(str(face.data[CUBA.LABEL]),
                             face_f.data[CUBA.LABEL])

    def test_iter_cells(self):
        """Test iter_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for cell_f in foam_mesh.iter_cells():
            cell = self.mesh.get_cell(cell_f.uid)
            self.assertEqual(cell.data[CUBA.VELOCITY],
                             cell_f.data[CUBA.VELOCITY])

    def test_has_cells(self):
        """Test has_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertTrue(foam_mesh.has_cells())

    def test_has_faces(self):
        """Test has_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertTrue(foam_mesh.has_faces())

    def test_count_of(self):
        """Test count_of method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)

        item_type = CUDSItem.POINT
        self.assertEqual(foam_mesh.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUDSItem.EDGE
        with self.assertRaises(ValueError):
            foam_mesh.count_of(item_type)

        item_type = CUDSItem.FACE
        self.assertEqual(foam_mesh.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUDSItem.CELL
        self.assertEqual(foam_mesh.count_of(item_type),
                         self.mesh.count_of(item_type))

    def test_generate_uuidmapping(self):
        """Test generate_uuidmapping method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        foam_mesh.generate_uuidmapping(8, 0, 6, 1)


if __name__ == '__main__':
    unittest.main()
