""" test_foam_mesh module

This module contains the unitary tests for the
foam_mesh module functionalities

"""

import unittest

from simphony.cuds.mesh import Mesh, Face, Point, Cell, Edge
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from foam_internalwrapper.foam_mesh import FoamMesh


class TestMesh(Mesh):
    """ Test mesh class to containing definition of boundaries
    """
    def __init__(self, name):
        Mesh.__init__(self, name)
        self._boundaries = {}


class FoamMeshTestCase(unittest.TestCase):
    """Test case for FoamMesh class"""
    def setUp(self):
        self.mesh = TestMesh(name='test_mesh')
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

        puids = self.mesh.add(self.points)

        self.faces = [
            Face([puids[0], puids[3], puids[7], puids[4]]),
            Face([puids[1], puids[2], puids[6], puids[5]]),
            Face([puids[0], puids[1], puids[5], puids[4]]),
            Face([puids[3], puids[2], puids[6], puids[7]]),
            Face([puids[0], puids[1], puids[2], puids[3]]),
            Face([puids[4], puids[5], puids[6], puids[7]])
        ]

        self.edges = [Edge([puids[0], puids[3]])]

        self.fuids = self.mesh.add(self.faces)

        self.cells = [
            Cell(puids,
                 data=DataContainer({CUBA.VELOCITY: [1, 0, 0],
                                     CUBA.PRESSURE: 4.0}))
        ]

        self.puids = puids

        self.mesh.add(self.cells)

        self.boundaries = {"boundary"+str(i): [self.fuids[i]]
                           for i in range(6)}
        self.mesh._boundaries = self.boundaries

    def test_get_point(self):
        """Test get_point method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        for point in self.mesh.iter(item_type=CUBA.POINT):
            puid = foam_mesh._foamPointLabelToUuid[label]
            point_f = foam_mesh.get(puid)
            self.assertEqual(point.coordinates, point_f.coordinates)
            label += 1

    def test_get_edge(self):
        """Test get_edge method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.get(self.edges[0].uid)

    def test_get_face(self):
        """Test get_face method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for face in foam_mesh.iter(item_type=CUBA.FACE):
            self.assertEqual(len(face.points), 4)

    def test_get_cell(self):
        """Test get_cell method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        for cell in self.mesh.iter(item_type=CUBA.CELL):
            cuid = foam_mesh._foamCellLabelToUuid[label]
            cell_f = foam_mesh.get(cuid)
            self.assertEqual(cell.data[CUBA.PRESSURE],
                             cell_f.data[CUBA.PRESSURE])
            self.assertEqual(cell.data[CUBA.VELOCITY],
                             cell_f.data[CUBA.VELOCITY])
            label += 1

    def test_add_points(self):
        """Test add_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add(self.points)

    def test_add_edges(self):
        """Test add_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add(self.edges)

    def test_add_faces(self):
        """Test add_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add(self.faces)

    def test_add_cells(self):
        """Test add_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.add(self.cells)

    def test_update_points(self):
        """Test update_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.update(self.points)

    def test_update_edges(self):
        """Test update_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.update(self.edges)

    def test_update_faces(self):
        """Test update_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        with self.assertRaises(NotImplementedError):
            foam_mesh.update(self.faces)

    def test_update_cells(self):
        """Test update_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        cuid = foam_mesh._foamCellLabelToUuid[label]
        cell_f = foam_mesh.get(cuid)
        self.assertIsInstance(cell_f.data, DataContainer)
        cell = list(self.mesh.iter(item_type=CUBA.CELL))[label]
        self.assertEqual(cell.data, cell_f.data)

        updated_cells = []
        for cell in foam_mesh.iter(item_type=CUBA.CELL):
            cell.data[CUBA.VELOCITY] = [2, 1, 3]
            updated_cells.append(cell)
        foam_mesh.update(updated_cells)

        label = 0
        cuid = foam_mesh._foamCellLabelToUuid[label]
        cell_f = foam_mesh.get(cuid)
        self.assertIsInstance(cell_f.data, DataContainer)
        cell = list(self.mesh.iter(item_type=CUBA.CELL))[label]
        self.assertNotEqual(cell.data, cell_f.data)

        updated_cells = []
        for cell in foam_mesh.iter(item_type=CUBA.CELL):
            cell.points = [self.points[1].uid, self.points[2].uid,
                           self.points[3].uid]
            updated_cells.append(cell)
        with self.assertRaises(Warning):
            foam_mesh.update(updated_cells)

    def test_iter_points(self):
        """Test iter_points method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        label = 0
        for point in self.mesh.iter(item_type=CUBA.POINT):
            puid = foam_mesh._foamPointLabelToUuid[label]
            point_f = foam_mesh.get(puid)
            self.assertEqual(point.coordinates, point_f.coordinates)
            label += 1

    def test_iter_edges(self):
        """Test iter_edges method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertEqual(foam_mesh.iter(item_type=CUBA.EDGE), [])

    def test_iter_faces(self):
        """Test iter_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for face_f in foam_mesh.iter(item_type=CUBA.FACE):
            self.assertEqual(len(face_f.points), 4)

    def test_iter_cells(self):
        """Test iter_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        for cell_f in foam_mesh.iter(item_type=CUBA.CELL):
            label = foam_mesh._uuidToFoamLabel[cell_f.uid]
            cell = self.cells[label]
            self.assertEqual(cell.data[CUBA.VELOCITY],
                             cell_f.data[CUBA.VELOCITY])

    def test_has_cells(self):
        """Test has_cells method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertTrue(foam_mesh.has_type(item_type=CUBA.CELL))

    def test_has_faces(self):
        """Test has_faces method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        self.assertTrue(foam_mesh.has_type(item_type=CUBA.FACE))

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

    def test_generate_uuidmapping(self):
        """Test generate_uuidmapping method

        """

        foam_mesh = FoamMesh('test_mesh', {}, self.solver, self.mesh)
        foam_mesh.generate_uuidmapping(8, 0, 6, 1)


if __name__ == '__main__':
    unittest.main()
