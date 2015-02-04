""" test_foamcontrolwrapper module

This module contains the unitary tests for the
foam_controlwrapper module functionalities

"""

import unittest
from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from simphony.core.cuba import CUBA
import logging
logger = logging.getLogger(__name__)


class FoamControlWrapperTestCase(unittest.TestCase):
    """Test case for FoamTokenWrapper class"""
    def setUp(self):
        pass

    def test_pitzDailyRun(self):
        """Test to run pitzDaily example in OpenFoam"""
        foam_controlwrapper = FoamControlWrapper()
        foam_controlwrapper.CM[CUBA.NAME]="foam_controlwrapper/tests/pitzDaily"
#        foam_controlwrapper.CM[CUBA.SOLVER]="simpleFoam"
        foam_controlwrapper.SP[CUBA.TIME_STEP]=1
        foam_controlwrapper.SP[CUBA.NUMBEROF_TIME_STEPS]=500
        foam_controlwrapper.SP[CUBA.DENSITY] = 1.0
        foam_controlwrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0e-5
        
        foam_controlwrapper.run()
        

if __name__ == '__main__':
    unittest.main()

