import unittest
import os
import shutil

from simphony.core.cuba import CUBA
from foam_controlwrapper.foam_controlwrapper import Wrapper
from foam_controlwrapper.cuba_extension import CUBAExt
from foam_controlwrapper.blockmesh_utils import create_quad_mesh


class WrapperRunTestCase(unittest.TestCase):
    def setUp(self):
        wrapper = Wrapper()
        self.path = "test_path_io"
        name = "simplemeshIO"
        wrapper.CM[CUBA.NAME] = name
        wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                             CUBAExt.LAMINAR_MODEL)
        wrapper.SP[CUBA.TIME_STEP] = 1
        wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1
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
        if os.path.exists(self.mesh_inside_wrapper.path):
            shutil.rmtree(self.mesh_inside_wrapper.path)
        if os.path.exists(self.path):
            shutil.rmtree(self.path)

    def test_run_time(self):
        """Test that field variable value is changed after
        consecutive calls of run method

        """

        self.wrapper.SP[CUBA.TIME_STEP] = 1
        self.wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 2

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


if __name__ == '__main__':
    unittest.main()
