""" test_foam_files module

This module contains the unitary tests for the
foam_files module functionalities

"""

import unittest

from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from foam_controlwrapper.foam_files import FoamFiles
from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from simphony.io.h5_cuds import H5CUDS
import os
import shutil


class FoamFilesTestCase(unittest.TestCase):
    """Test case for FoamFiles class"""
    def setUp(self):
        self.files = FoamFiles()
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
        FoamFiles().remove_parser_files(os.getcwd())
        if os.path.exists(self.path):
            shutil.rmtree(self.path)

    def test_create_file_content(self):
        """Test create_file_content method

        """

        self.files.create_file_content(self.path,
                                       self.solver,
                                       self.time,
                                       True)

    def test_create_directories(self):
        """Test create_directories method

        """

        self.files.create_directories(self.path)

    def test_write_default_files(self):
        """Test write_default_files method

        """

        self.files.write_default_files(self.path,
                                       self.solver,
                                       self.time,
                                       True)

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

        self.files.write_default_files(mesh_inside_wrapper.path,
                                       self.solver,
                                       self.time, True)
        self.files.modify_files(mesh_inside_wrapper.path, self.time, self.SP,
                                self.BC, self.solver,
                                self.SPExt, self.CMExt)
        mesh_file.close()


if __name__ == '__main__':
    unittest.main()
