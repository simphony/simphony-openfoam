import unittest
import os
import re
import shutil

from foam_controlwrapper.foam_controlwrapper import Wrapper
from simphony.core.cuba import CUBA
from foam_controlwrapper.cuba_extension import CUBAExt
from foam_controlwrapper.blockmesh_utils import create_quad_mesh


class WrapperRunTestCase(unittest.TestCase):
    def setUp(self):
        wrapper = Wrapper()
        name = "simplemesh_parallel"
        self.path = "test_parallel_path"
        wrapper.CM[CUBA.NAME] = name
        wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                             CUBAExt.LAMINAR_MODEL)

        wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 4

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

    def test_parallel_run(self):
        """Test parallel running of OpenFoam

        """

        self.wrapper.run()

        self.assertEqual(self.wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES],
                         len([d for d in
                              os.listdir(self.mesh_inside_wrapper.path)
                              if re.match(r'processor*', d)]))

if __name__ == '__main__':
    unittest.main()
