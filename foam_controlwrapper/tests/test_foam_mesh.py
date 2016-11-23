""" test_foam_mesh module

This module contains the unitary tests for the
foam_mesh module functionalities

"""

import unittest
import os

from simphony.cuds.mesh import Mesh, Face, Point, Cell, Edge
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from foam_controlwrapper.foam_mesh import FoamMesh


class TestMesh(Mesh):
    """ Test mesh class to containing definition of boundaries
    """
    def __init__(self, name):
        Mesh.__init__(self, name)
        self._boundaries = {}


class FoamMeshTestCase(unittest.TestCase):
    """Test case for FoamMesh class"""
    def setUp(self):

        self.mesh = TestMesh(name="mesh1")
        self.solver = 'pimpleFoam'

        self.points = [
            Point(
                (0.0, 0.0, 0.0)),
            Point(
                (1.0, 0.0, 0.0)),
            Point(
                (1.0, 0.0, 1.0)),
            Point(
                (0.0, 0.0, 1.0)),
            Point(
                (0.0, 1.0, 0.0)),
            Point(
                (1.0, 1.0, 0.0)),
            Point(
                (1.0, 1.0, 1.0)),
            Point(
                (0.0, 1.0, 1.0))
        ]

        puids = self.mesh._add_points(self.points)

        self.faces = [
            Face([puids[0], puids[3], puids[7], puids[4]]),
            Face([puids[1], puids[2], puids[6], puids[5]]),
            Face([puids[0], puids[1], puids[5], puids[4]]),
            Face([puids[3], puids[2], puids[6], puids[7]]),
            Face([puids[0], puids[1], puids[2], puids[3]]),
            Face([puids[4], puids[5], puids[6], puids[7]])
        ]

        self.edges = [Edge([puids[0], puids[3]])]

        self.fuids = self.mesh._add_faces(self.faces)

        self.cells = [
            Cell(puids,
                 data=DataContainer({CUBA.VELOCITY: [1, 0, 0],
                                     CUBA.PRESSURE: 4.0}))
        ]

        self.puids = puids

        self.mesh._add_cells(self.cells)
        self.boundaries = {"boundary"+str(i): [self.fuids[i]]
                           for i in range(6)}
        self.mesh._boundaries = self.boundaries

    def test_get_point(self):
        """Test _get_point method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        for point in self.mesh._iter_points():
            puid = foam_mesh._foamPointLabelToUuid[label]
            point_f = foam_mesh._get_point(puid)
            self.assertEqual(point.coordinates, point_f.coordinates)
            label += 1

    def test_get_edge(self):
        """Test _get_edge method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._get_edge(self.edges[0].uid)

    def test_get_face(self):
        """Test _get_face method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for face in foam_mesh._iter_faces():
            self.assertEqual(len(face.points), 4)

    def test_get_cell(self):
        """Test get_cell method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        for cell in self.mesh._iter_cells():
            cuid = foam_mesh._foamCellLabelToUuid[label]
            cell_f = foam_mesh._get_cell(cuid)
            self.assertEqual(cell.data[CUBA.PRESSURE],
                             cell_f.data[CUBA.PRESSURE])
            self.assertEqual(cell.data[CUBA.VELOCITY],
                             cell_f.data[CUBA.VELOCITY])
            label += 1

    def test_add_points(self):
        """Test _add_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._add_points(self.points)

    def test_add_edges(self):
        """Test _add_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._add_edges(self.edges)

    def test_add_faces(self):
        """Test _add_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._add_faces(self.faces)

    def test_add_cells(self):
        """Test _add_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._add_cells(self.cells)

    def test_update_points(self):
        """Test _update_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._update_points(self.points)

    def test_update_edges(self):
        """Test _update_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._update_edges(self.edges)

    def test_update_faces(self):
        """Test _update_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh._update_faces(self.faces)

    def test_update_cells(self):
        """Test _update_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        cuid = foam_mesh._foamCellLabelToUuid[label]
        cell_f = foam_mesh._get_cell(cuid)
        self.assertIsInstance(cell_f.data, DataContainer)
        cell = list(self.mesh._iter_cells())[label]
        self.assertEqual(cell.data, cell_f.data)

        updated_cells = []
        for cell in foam_mesh._iter_cells():
            cell.data[CUBA.VELOCITY] = [2, 1, 3]
            updated_cells.append(cell)
        foam_mesh._update_cells(updated_cells)

        label = 0
        cuid = foam_mesh._foamCellLabelToUuid[label]
        cell_f = foam_mesh._get_cell(cuid)
        self.assertIsInstance(cell_f.data, DataContainer)
        cell = list(self.mesh._iter_cells())[label]
        self.assertNotEqual(cell.data, cell_f.data)

        updated_cells = []
        for cell in foam_mesh._iter_cells():
            cell.points = [self.points[1].uid, self.points[2].uid,
                           self.points[3].uid]
            updated_cells.append(cell)
        with self.assertRaises(Warning):
            foam_mesh._update_cells(updated_cells)

    def test_iter_points(self):
        """Test _iter_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        for point in self.mesh._iter_points():
            puid = foam_mesh._foamPointLabelToUuid[label]
            point_f = foam_mesh._get_point(puid)
            self.assertEqual(point.coordinates, point_f.coordinates)
            label += 1

    def test_iter_edges(self):
        """Test _iter_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertEqual(foam_mesh._iter_edges(), [])

    def test_iter_faces(self):
        """Test _iter_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for face_f in foam_mesh._iter_faces():
            self.assertEqual(len(face_f.points), 4)

    def test_iter_cells(self):
        """Test _iter_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for cell_f in foam_mesh._iter_cells():
            label = foam_mesh._uuidToFoamLabel[cell_f.uid]
            cell = self.cells[label]
            self.assertEqual(cell.data[CUBA.VELOCITY],
                             cell_f.data[CUBA.VELOCITY])

    def test_has_faces(self):
        """Test _has_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertTrue(foam_mesh._has_faces())

    def test_has_cells(self):
        """Test _has_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertTrue(foam_mesh._has_cells())

    def test_count_of(self):
        """Test count_of method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)

        item_type = CUBA.POINT
        self.assertEqual(foam_mesh.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUBA.EDGE
        self.assertEqual(foam_mesh.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUBA.FACE
        self.assertEqual(foam_mesh.count_of(item_type),
                         self.mesh.count_of(item_type))

        item_type = CUBA.CELL
        self.assertEqual(foam_mesh.count_of(item_type),
                         self.mesh.count_of(item_type))

    def test_write(self):
        """Test write method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        foam_mesh.write()
        meshpath = os.path.join(foam_mesh.path, 'constant', 'polyMesh')
        self.assertTrue(os.path.exists(os.path.join(meshpath, 'points')))
        self.assertTrue(os.path.exists(os.path.join(meshpath, 'owner')))
        self.assertTrue(os.path.exists(os.path.join(meshpath, 'neighbour')))
        self.assertTrue(os.path.exists(os.path.join(meshpath, 'boundary')))
        self.assertTrue(os.path.exists(os.path.join(meshpath, 'faces')))

    def test_generate_uuidmapping(self):
        """Test generate_uuidmapping method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        foam_mesh.generate_uuidmapping(8, 0, 6, 1)


if __name__ == '__main__':
    unittest.main()
