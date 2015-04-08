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
from simphony.io.h5_cuds import H5CUDS

from foam_controlwrapper import foam_files
from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper


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
        self.assertEqual('controlDict' in
                         re.findall(r'(\w+)',
                                    content[file_name]), True)

    def test_create_directories(self):
        """Test create_directories method

        """

        foam_files.create_directories(self.path)
        self.assertEqual(os.path.exists(
            os.path.join(self.path, 'system')), True)

    def test_write_default_files(self):
        """Test write_default_files method

        """

        foam_files.write_default_files(self.path,
                                       self.solver,
                                       self.time,
                                       True)
        self.assertEqual(os.path.exists(
            os.path.join(self.path, 'system', 'controlDict')), True)

    def test_modify_files(self):
        """Test modify_files method

        """

        wrapper = FoamControlWrapper()
        name = 'simplemesh'
        mesh_file = H5CUDS.open(os.path.join('foam_controlwrapper',
                                             'tests',
                                             'simplemesh.cuds'))
        mesh_from_file = mesh_file.get_mesh(name)

        mesh_inside_wrapper = wrapper.add_mesh(mesh_from_file)

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

        mesh_file.close()


if __name__ == '__main__':
    unittest.main()
