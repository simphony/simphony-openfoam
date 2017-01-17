""" test_foaminternalwrapper module

This module contains the unitary tests for the
foam_controlwrapper module functionalities

"""

import unittest
import os
import shutil

from simphony.cuds.mesh import Mesh, Face, Point, Cell
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from foam_internalwrapper.foam_internalwrapper import Wrapper
from foam_controlwrapper.cuba_extension import CUBAExt

from foam_controlwrapper.blockmesh_utils import create_quad_mesh


class TestMesh(Mesh):
    """ Test mesh class to containing definition of boundaries
    """
    def __init__(self, name):
        Mesh.__init__(self, name)
        self._boundaries = {}


class WrapperTestCase(unittest.TestCase):
    """Test case for Wrapper class"""
    def setUp(self):
        self.mesh = TestMesh(name="mesh1")
        self.GE = (CUBAExt.INCOMPRESSIBLE, CUBAExt.LAMINAR_MODEL)
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

        self.fuids = self.mesh.add_faces(self.faces)

        self.cells = [
            Cell(puids,
                 data=DataContainer({CUBA.VELOCITY: [1, 0, 0],
                                     CUBA.PRESSURE: 4.0}))
        ]

        self.puids = puids

        self.mesh.add_cells(self.cells)

        self.boundaries = {"boundary"+str(i): [self.fuids[i]]
                           for i in range(6)}
        self.mesh._boundaries = self.boundaries

    def test_add_dataset(self):
        """Test add_dataset method

        """

        wrapper = Wrapper()
        wrapper.CM_extensions[CUBAExt.GE] = self.GE
        wrapper.add_dataset(self.mesh)
        self.assertEqual(sum(1 for _ in wrapper.iter_datasets()), 1)

    def test_remove_dataset(self):
        """Test remove_dataset method

        """

        wrapper = Wrapper()
        wrapper.CM_extensions[CUBAExt.GE] = self.GE
        wrapper.add_dataset(self.mesh)
        wrapper.remove_dataset(self.mesh.name)
        with self.assertRaises(ValueError):
            wrapper.get_dataset(self.mesh.name)

    def test_get_dataset(self):
        """Test get_dataset method

        """

        wrapper = Wrapper()
        wrapper.CM_extensions[CUBAExt.GE] = self.GE
        wrapper.add_dataset(self.mesh)
        mesh_inside_wrapper = wrapper.get_dataset(self.mesh.name)
        self.assertEqual(self.mesh.name, mesh_inside_wrapper.name)

        label = 0
        for point in self.mesh.iter(item_type=CUBA.POINT):
            puid = mesh_inside_wrapper._foamPointLabelToUuid[label]
            point_f = mesh_inside_wrapper.get(puid)
            self.assertEqual(point.coordinates, point_f.coordinates)
            label += 1

        label = 0
        for cell in self.mesh.iter(item_type=CUBA.CELL):
            cuid = mesh_inside_wrapper._foamCellLabelToUuid[label]
            cell_f = mesh_inside_wrapper.get(cuid)
            self.assertEqual(cell.data[CUBA.PRESSURE],
                             cell_f.data[CUBA.PRESSURE])
            self.assertEqual(cell.data[CUBA.VELOCITY],
                             cell_f.data[CUBA.VELOCITY])
            label += 1

    def test_iter_datasets(self):
        """Test iter_datasets method

        """

        wrapper = Wrapper()
        wrapper.CM_extensions[CUBAExt.GE] = self.GE
        wrapper.add_dataset(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_dataset(mesh2)

        self.assertEqual(sum(1 for _ in wrapper.iter_datasets()), 2)

    def test_multiple_meshes(self):
        """Test multiple meshes inside wrapper

        """

        wrapper = Wrapper()
        wrapper.CM_extensions[CUBAExt.GE] = self.GE
        wrapper.add_dataset(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_dataset(mesh2)
        mesh_inside_wrapper1 = wrapper.get_dataset(self.mesh.name)
        mesh_inside_wrapper2 = wrapper.get_dataset(mesh2.name)

        self.assertEqual(
            sum(1 for _ in mesh_inside_wrapper1.iter(item_type=CUBA.POINT)),
            sum(1 for _ in mesh_inside_wrapper2.iter(item_type=CUBA.POINT)))


class WrapperRunTestCase(unittest.TestCase):
    def setUp(self):
        wrapper = Wrapper()
        self.path = "test_path"
        name = "simplemesh"
        wrapper.CM[CUBA.NAME] = name
        wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                             CUBAExt.LAMINAR_MODEL)
        wrapper.SP[CUBA.TIME_STEP] = 1
        wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 3
        wrapper.SP[CUBA.DENSITY] = 1.0
        wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0
        wrapper.BC[CUBA.VELOCITY] = {'inlet': ('fixedValue', (0.1, 0, 0)),
                                     'outlet': 'zeroGradient',
                                     'walls': ('fixedValue', (0, 0, 0)),
                                     'frontAndBack': 'empty'}
        wrapper.BC[CUBA.PRESSURE] = {'inlet': 'zeroGradient',
                                     'outlet': ('fixedValue', 0),
                                     'walls': 'zeroGradient',
                                     'frontAndBack': 'empty'}
        self.wrapper = wrapper

        corner_points = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
                         (5.0, 5.0, 0.0), (0.0, 5.0, 0.0),
                         (0.0, 0.0, 1.0), (5.0, 0.0, 1.0),
                         (5.0, 5.0, 1.0), (0.0, 5.0, 1.0)]
        create_quad_mesh(self.path, name, self.wrapper, corner_points, 5, 5, 5)
        self.mesh_inside_wrapper = self.wrapper.get_dataset(name)

    def tearDown(self):
        if os.path.exists(self.path):
            shutil.rmtree(self.path)

    def test_run_time(self):
        """Test that field variable value is changed after
        consecutive calls of run method

        """
        self.wrapper.SP[CUBA.TIME_STEP] = 1

        self.wrapper.run()

        for cell in self.mesh_inside_wrapper.iter(item_type=CUBA.CELL):
            old_vel = cell.data[CUBA.VELOCITY]
            old_pres = cell.data[CUBA.PRESSURE]
            cell_uid = cell.uid

        self.wrapper.run()

        cell = self.mesh_inside_wrapper.get_cell(cell_uid)
        new_vel = cell.data[CUBA.VELOCITY]
        new_pres = cell.data[CUBA.PRESSURE]

        self.assertNotEqual(old_vel, new_vel)
        self.assertNotEqual(old_pres, new_pres)


if __name__ == '__main__':
    unittest.main()
