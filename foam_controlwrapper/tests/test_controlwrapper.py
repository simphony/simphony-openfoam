""" test_foamcontrolwrapper module

This module contains the unitary tests for the
foam_controlwrapper module functionalities

"""

import unittest
import os
import re
import shutil

from simphony.cuds.mesh import Mesh, Face, Point, Cell
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from foam_controlwrapper.cuba_extension import CUBAExt
from foam_controlwrapper.blockmesh_utils import create_quad_mesh
from foam_controlwrapper import foam_files


class FoamControlWrapperTestCase(unittest.TestCase):
    """Test case for FoamControlWrapper class"""
    def setUp(self):
        self.mesh = Mesh(name="mesh1")

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

        self.mesh.add_faces(self.faces)

        self.cells = [
            Cell(puids)
        ]

        self.puids = puids

        self.mesh.add_cells(self.cells)

    def tearDown(self):
        foam_files.remove_parser_files(os.getcwd())

    def test_add_dataset(self):
        """Test add_dataset method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_dataset(self.mesh)
        self.assertEqual(sum(1 for _ in wrapper.iter_datasets()), 1)

    def test_remove_dataset(self):
        """Test remove_dataset method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_dataset(self.mesh)
        wrapper.remove_dataset(self.mesh.name)
        with self.assertRaises(ValueError):
            wrapper.get_dataset(self.mesh.name)

    def test_get_dataset(self):
        """Test get_dataset method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_dataset(self.mesh)
        mesh_inside_wrapper = wrapper.get_dataset(self.mesh.name)
        self.assertEqual(self.mesh.name, mesh_inside_wrapper.name)

        for point in self.mesh.iter_points():
            point_w = mesh_inside_wrapper.get_point(point.uid)
            self.assertEqual(point.uid, point_w.uid)
            self.assertEqual(point.coordinates, point_w.coordinates)

        for face in self.mesh.iter_faces():
            face_w = mesh_inside_wrapper.get_face(face.uid)
            self.assertEqual(face.uid, face_w.uid)
            self.assertEqual(face.points, face_w.points)
            self.assertEqual(face.data, face_w.data)

        for cell in self.mesh.iter_cells():
            cell_w = mesh_inside_wrapper.get_cell(cell.uid)
            self.assertEqual(cell.uid, cell_w.uid)
            self.assertEqual(set(cell.points), set(cell_w.points))

    def test_get_dataset_names(self):
        """Test get_dataset_names method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_dataset(self.mesh)
        name1 = self.mesh.name
        mesh2 = self.mesh
        name2 = "mesh2"
        mesh2.name = name2
        wrapper.add_dataset(mesh2)

        self.assertEqual(list(wrapper.get_dataset_names())[0], name1)
        self.assertEqual(list(wrapper.get_dataset_names())[1], name2)

    def test_iter_datasets(self):
        """Test iter_datasets method

        """

        wrapper = FoamControlWrapper()
        wrapper.add_dataset(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_dataset(mesh2)

        self.assertEqual(sum(1 for _ in wrapper.iter_datasets()), 2)

    def test_multiple_meshes(self):
        """Test multiple meshes inside wrapper

        """

        wrapper = FoamControlWrapper()
        wrapper.add_dataset(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_dataset(mesh2)
        mesh_inside_wrapper1 = wrapper.get_dataset(self.mesh.name)
        mesh_inside_wrapper2 = wrapper.get_dataset(mesh2.name)

        self.assertEqual(
            sum(1 for _ in mesh_inside_wrapper1.iter_points()),
            sum(1 for _ in mesh_inside_wrapper2.iter_points()))


class FoamControlWrapperRunTestCase(unittest.TestCase):
    def setUp(self):
        wrapper = FoamControlWrapper()
        path = "test_path"
        name = "simplemesh"
        wrapper.CM[CUBA.NAME] = name
        wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                             CUBAExt.LAMINAR_MODEL)
        wrapper.SP[CUBA.TIME_STEP] = 1
        wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 3
        wrapper.SP[CUBA.DENSITY] = 1.0
        wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0
        wrapper.BC[CUBA.VELOCITY] = {'boundary0': (0.1, 0, 0),
                                     'boundary1': 'zeroGradient',
                                     'boundary2': (0, 0, 0),
                                     'boundary3': 'empty'}
        wrapper.BC[CUBA.PRESSURE] = {'boundary0': 'zeroGradient',
                                     'boundary1': 0,
                                     'boundary2': 'zeroGradient',
                                     'boundary3': 'empty'}
        self.wrapper = wrapper

        corner_points = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
                         (5.0, 5.0, 0.0), (0.0, 5.0, 0.0),
                         (0.0, 0.0, 1.0), (5.0, 0.0, 1.0),
                         (5.0, 5.0, 1.0), (0.0, 5.0, 1.0)]
        create_quad_mesh(path, name, self.wrapper, corner_points, 5, 5, 5)
        self.mesh_inside_wrapper = self.wrapper.get_dataset(name)

    def tearDown(self):
        if os.path.exists(self.mesh_inside_wrapper.path):
            shutil.rmtree(self.mesh_inside_wrapper.path)

    def test_run_time(self):
        """Test that field variable value is changed after
        consecutive calls of run method

        """
        self.wrapper.SP[CUBA.TIME_STEP] = 1

        self.wrapper.run()

        for cell in self.mesh_inside_wrapper.iter_cells():
            old_vel = cell.data[CUBA.VELOCITY]
            old_pres = cell.data[CUBA.PRESSURE]
            cell_uid = cell.uid

        self.wrapper.run()

        cell = self.mesh_inside_wrapper.get_cell(cell_uid)
        new_vel = cell.data[CUBA.VELOCITY]
        new_pres = cell.data[CUBA.PRESSURE]

        self.assertNotEqual(old_vel, new_vel)
        self.assertNotEqual(old_pres, new_pres)

    def test_parallel_run(self):
        """Test parallel running of OpenFoam

        """
        self.wrapper.SP[CUBA.TIME_STEP] = 3
        self.wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 2

        self.wrapper.run()

        os.listdir(self.mesh_inside_wrapper.path)
        self.assertEqual(self.wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES],
                         len([d for d in
                              os.listdir(self.mesh_inside_wrapper.path)
                              if re.match(r'processor*', d)]))


if __name__ == '__main__':
    unittest.main()
