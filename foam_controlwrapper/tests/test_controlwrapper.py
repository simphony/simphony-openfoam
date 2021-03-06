""" test_foamcontrolwrapper module

This module contains the unitary tests for the
foam_controlwrapper module functionalities

"""

import unittest

import foam_controlwrapper as register

from simphony.cuds.mesh import Mesh, Face, Point, Cell
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.engine import EngineInterface
from simphony.api import CUDS, Simulation
from simphony.cuds.meta.api import Cfd


class TestMesh(Mesh):
    """ Test mesh class to containing definition of boundaries
    """
    def __init__(self, name):
        Mesh.__init__(self, name)
        self._boundaries = {}


class WrapperTestCase(unittest.TestCase):
    """Test case for Wrapper class"""
    def setUp(self):
        register.OpenFOAMExtension()

        self.engine = EngineInterface.FileIO

        self.mesh = TestMesh(name="mesh1")

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

        self.cuds = CUDS(name='cuds')
        self.cuds.add([Cfd(name='default cfd model'), self.mesh])

    def test_add_dataset(self):
        """Test add_dataset method

        """

        sim = Simulation(self.cuds, 'OpenFOAM',
                         engine_interface=self.engine)
        wrapper = sim._wrapper

        self.assertEqual(sum(1 for _ in wrapper.iter_datasets()), 1)

    def test_remove_dataset(self):
        """Test remove_dataset method

        """

        sim = Simulation(self.cuds, 'OpenFOAM',
                         engine_interface=self.engine)
        wrapper = sim._wrapper
        wrapper.remove_dataset(self.mesh.name)
        with self.assertRaises(ValueError):
            wrapper.get_dataset(self.mesh.name)

    def test_get_dataset(self):
        """Test get_dataset method

        """

        sim = Simulation(self.cuds, 'OpenFOAM',
                         engine_interface=self.engine)
        wrapper = sim._wrapper
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

    def test_get_dataset_names(self):
        """Test get_dataset_names method

        """

        sim = Simulation(self.cuds, 'OpenFOAM',
                         engine_interface=self.engine)
        wrapper = sim._wrapper
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

        sim = Simulation(self.cuds, 'OpenFOAM',
                         engine_interface=self.engine)
        wrapper = sim._wrapper
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_dataset(mesh2)

        self.assertEqual(sum(1 for _ in wrapper.iter_datasets()), 2)

    def test_multiple_meshes(self):
        """Test multiple meshes inside wrapper

        """

        sim = Simulation(self.cuds, 'OpenFOAM',
                         engine_interface=self.engine)
        wrapper = sim._wrapper
        mesh2 = self.mesh
        mesh2.name = "mesh2"
        wrapper.add_dataset(mesh2)
        mesh_inside_wrapper1 = wrapper.get_dataset(self.mesh.name)
        mesh_inside_wrapper2 = wrapper.get_dataset(mesh2.name)

        self.assertEqual(
            sum(1 for _ in mesh_inside_wrapper1.iter(item_type=CUBA.POINT)),
            sum(1 for _ in mesh_inside_wrapper2.iter(item_type=CUBA.POINT)))


if __name__ == '__main__':
    unittest.main()
