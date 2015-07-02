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
from simphony.io.h5_cuds import H5CUDS

from foam_internalwrapper.foam_internalwrapper import FoamInternalWrapper
from foam_internalwrapper.cuba_extension import CUBAExt


class FoamInternalWrapperTestCase(unittest.TestCase):
    """Test case for FoamInternalWrapper class"""
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

        self.puids = puids

        [self.mesh.add_cell(cell) for cell in self.cells]

    def test_add_mesh(self):
        """Test add_mesh method

        """

        wrapper = FoamInternalWrapper()
        wrapper.add_mesh(self.mesh)
        self.assertEqual(sum(1 for _ in wrapper.iter_meshes()), 1)

    def test_delete_mesh(self):
        """Test delete_mesh method

        """

        wrapper = FoamInternalWrapper()
        wrapper.add_mesh(self.mesh)
        wrapper.delete_mesh(self.mesh.name)
        with self.assertRaises(ValueError):
            wrapper.get_mesh(self.mesh.name)

    def test_get_mesh(self):
        """Test get_mesh method

        """

        wrapper = FoamInternalWrapper()
        wrapper.add_mesh(self.mesh)
        mesh_inside_wrapper = wrapper.get_mesh(self.mesh.name)
        self.assertEqual(self.mesh.name, mesh_inside_wrapper.name)

        for point in self.mesh.iter_points():
            point_w = mesh_inside_wrapper.get_point(point.uid)
            self.assertEqual(point.uid, point_w.uid)
            self.assertEqual(point.coordinates, point_w.coordinates)

        for face in self.mesh.iter_faces():
            face_w = mesh_inside_wrapper.get_face(face.uid)
            self.assertEqual(face.uid, face_w.uid)

        for cell in self.mesh.iter_cells():
            cell_w = mesh_inside_wrapper.get_cell(cell.uid)
            self.assertEqual(cell.uid, cell_w.uid)
            self.assertEqual(set(cell.points), set(cell_w.points))

    def test_iter_meshes(self):
        """Test iter_meshes method

        """

        wrapper = FoamInternalWrapper()
        wrapper.add_mesh(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_mesh(mesh2)

        self.assertEqual(sum(1 for _ in wrapper.iter_meshes()), 2)

    def test_multiple_meshes(self):
        """Test multiple meshes inside wrapper

        """

        wrapper = FoamInternalWrapper()
        wrapper.add_mesh(self.mesh)
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_mesh(mesh2)
        mesh_inside_wrapper1 = wrapper.get_mesh(self.mesh.name)
        mesh_inside_wrapper2 = wrapper.get_mesh(mesh2.name)

        self.assertEqual(
            sum(1 for _ in mesh_inside_wrapper1.iter_points()),
            sum(1 for _ in mesh_inside_wrapper2.iter_points()))

    def test_add_particles(self):
        """Test add_particles method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.add_particles(DataContainer())

    def test_get_particles(self):
        """Test get_particles method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.get_particles('')

    def test_delete_particles(self):
        """Test delete_particles method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.delete_particles('')

    def test_iter_particles(self):
        """Test iter_particles method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.iter_particles()

    def test_add_lattice(self):
        """Test add_lattice method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.add_lattice('')

    def test_get_lattice(self):
        """Test get_lattice method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.get_lattice('')

    def test_delete_lattice(self):
        """Test delete_lattice method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.delete_lattice('')

    def test_iter_lattices(self):
        """Test iter_lattices method

        """

        wrapper = FoamInternalWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.iter_lattices()

    def test_run_time(self):
        """Test that field variable value is changed after
        consecutive calls of run method

        """

        wrapper = FoamInternalWrapper()
        name = 'simplemesh'
        wrapper.CM[CUBA.NAME] = name
        wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                             CUBAExt.LAMINAR_MODEL)
        wrapper.SP[CUBA.TIME_STEP] = 1
        wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1
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
        mesh_file = H5CUDS.open(os.path.join('foam_internalwrapper',
                                             'tests',
                                             'simplemesh.cuds'))
        mesh_from_file = mesh_file.get_mesh(name)

        mesh_inside_wrapper = wrapper.add_mesh(mesh_from_file)

        wrapper.run()

        for cell in mesh_inside_wrapper.iter_cells():
            old_vel = cell.data[CUBA.VELOCITY]
            old_pres = cell.data[CUBA.PRESSURE]
            cell_uid = cell.uid

        wrapper.run()

        cell = mesh_inside_wrapper.get_cell(cell_uid)
        new_vel = cell.data[CUBA.VELOCITY]
        new_pres = cell.data[CUBA.PRESSURE]

        self.assertNotEqual(old_vel, new_vel)
        self.assertNotEqual(old_pres, new_pres)

        if os.path.exists(mesh_inside_wrapper.path):
            shutil.rmtree(mesh_inside_wrapper.path)

        mesh_file.close()

if __name__ == '__main__':
    unittest.main()
