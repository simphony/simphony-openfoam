""" test_foamcontrolwrapper module

This module contains the unitary tests for the
foam_controlwrapper module functionalities

"""

import unittest
from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from simphony.core.cuba import CUBA
from foam_controlwrapper.cuba_extension import CUBAExt
import shutil
import logging
import os
from PyFoam.Execution.UtilityRunner import UtilityRunner

logger = logging.getLogger(__name__)


class FoamControlWrapperTestCase(unittest.TestCase):
    """Test case for FoamControlWrapper class"""
    def setUp(self):
        pass

    def test_pitzDailyRun(self):
        """Test to run pitzDaily laminar example in OpenFoam"
        """
        caseName = os.path.join('foam_controlwrapper', 'tests', 'pitzDaily')
        try:
            tutorials = os.path.expandvars('$FOAM_TUTORIALS')
            src = os.path.join(tutorials, 'incompressible',
                               'simpleFoam', 'pitzDaily')
            shutil.copytree(src, caseName)
        except shutil.Error as e:
            raise shutil.Error('Case directory not copied. Error: %s' % e)
        except OSError as e:
            raise OSError('Case directory not copied. Error: %s' % e)

        blockMesh = UtilityRunner(argv=["blockMesh", "-case", caseName],
                                  silent=True,
                                  logname="blockMesh")
        blockMesh.start()

        foam_controlwrapper = FoamControlWrapper()
        foam_controlwrapper.CM[CUBA.NAME] = caseName
        foam_controlwrapper.CM[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                              CUBAExt.LAMINAR_MODEL)
        foam_controlwrapper.SP[CUBA.TIME_STEP] = 1
        foam_controlwrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 100
        foam_controlwrapper.SP[CUBA.DENSITY] = 1.0
        foam_controlwrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0
# this is just an example. It is not enough for general setting of BC's
        foam_controlwrapper.BC[CUBA.VELOCITY] = {'inlet': (0.1, 0, 0),
                                                 'outlet': 'zeroGradient',
                                                 'upperWall': (0, 0, 0),
                                                 'lowerWall': (0, 0, 0),
                                                 'frontAndBack': 'empty'}
        foam_controlwrapper.BC[CUBA.PRESSURE] = {'inlet': 'zeroGradient',
                                                 'outlet': 0,
                                                 'upperWall': 'zeroGradient',
                                                 'lowerWall': 'zeroGradient',
                                                 'frontAndBack': 'empty'}

        timeStep = foam_controlwrapper.run()

        self.assertEqual(timeStep, '100')


if __name__ == '__main__':
    unittest.main()
