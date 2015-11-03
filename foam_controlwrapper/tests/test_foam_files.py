""" test_foam_files module

This module contains the unitary tests for the
foam_files module functionalities

"""

import unittest
import os
import re
import shutil

from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer


from foam_controlwrapper import foam_files
from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from foam_controlwrapper.blockmesh_utils import create_quad_mesh


class FoamFilesTestCase(unittest.TestCase):
    """Test case for FoamFiles class"""
    def setUp(self):
        self.path = 'test_path'
        self.solver = 'simpleFoam'
        self.time = '0'
        self.SP = DataContainer()
        self.BC = DataContainer()
        self.SPExt = {}
        self.CMExt = {}

        self.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1
        self.SP[CUBA.TIME_STEP] = 1
        self.SP[CUBA.DENSITY] = 1.0
        self.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0
        self.BC[CUBA.VELOCITY] = {'boundary0': (0.1, 0, 0),
                                  'boundary1': 'zeroGradient',
                                  'boundary2': (0, 0, 0),
                                  'boundary3': 'empty'}
        self.BC[CUBA.PRESSURE] = {'boundary0': 'zeroGradient',
                                  'boundary1': 0,
                                  'boundary2': 'zeroGradient',
                                  'boundary3': 'empty'}

    def tearDown(self):
        foam_files.remove_parser_files(os.getcwd())
        if os.path.exists(self.path):
            shutil.rmtree(self.path)

    def test_create_file_content(self):
        """Test create_file_content method

        """

        content = foam_files.create_file_content(self.path,
                                                 self.solver,
                                                 self.time,
                                                 True)
        file_name = os.path.join('system', 'controlDict')
        self.assertTrue('controlDict' in
                        re.findall(r'(\w+)', content[file_name]))

    def test_create_directories(self):
        """Test create_directories method

        """

        foam_files.create_directories(self.path)
        self.assertTrue(os.path.exists(os.path.join(self.path, 'system')))

    def test_write_default_files(self):
        """Test write_default_files method

        """

        foam_files.write_default_files(self.path,
                                       self.solver,
                                       self.time,
                                       True)
        self.assertTrue(os.path.exists(
            os.path.join(self.path, 'system', 'controlDict')))

    def test_modify_files(self):
        """Test modify_files method

        """

        wrapper = FoamControlWrapper()
        path = "test_path"
        name = "test_mesh"
        corner_points = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
                         (5.0, 5.0, 0.0), (0.0, 5.0, 0.0),
                         (0.0, 0.0, 1.0), (5.0, 0.0, 1.0),
                         (5.0, 5.0, 1.0), (0.0, 5.0, 1.0)]
        create_quad_mesh(path, name, wrapper, corner_points, 5, 5, 5)

        mesh_inside_wrapper = wrapper.get_dataset(name)

        foam_files.write_default_files(mesh_inside_wrapper.path,
                                       self.solver,
                                       self.time, True)

        file_name = os.path.join(mesh_inside_wrapper.path, self.time, 'U')
        with open(file_name, "r") as ufile:
            data = ufile.read()

        foam_files.modify_files(mesh_inside_wrapper.path, self.time, self.SP,
                                self.BC, self.solver,
                                self.SPExt, self.CMExt)
        with open(file_name, "r") as ufile:
            moddata = ufile.read()

        self.assertNotEqual(data, moddata)


if __name__ == '__main__':
    unittest.main()
