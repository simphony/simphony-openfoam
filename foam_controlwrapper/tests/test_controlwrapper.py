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
    """Test case for FoamControlWrapper class"""
    def setUp(self):
        pass

    def test_pitzDailyRun(self):
        """Test to run pitzDaily example in OpenFoam"
            -to get example working copy directory $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily
             to test -directory and run blockMesh on that directory
        """
        foam_controlwrapper = FoamControlWrapper()
        foam_controlwrapper.CM[CUBA.NAME]="foam_controlwrapper/tests/pitzDaily"
#        foam_controlwrapper.CM[CUBA.SOLVER]="simpleFoam"
        foam_controlwrapper.SP[CUBA.TIME_STEP]=1
        foam_controlwrapper.SP[CUBA.NUMBEROF_TIME_STEPS]=500
        foam_controlwrapper.SP[CUBA.DENSITY] = 1.0
        foam_controlwrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0e-5
# this is just an example. It is not enough for general setting of BC's
        foam_controlwrapper.BC[CUBA.VELOCITY]={ 'inlet'        : (10, 0, 0), 
                                                'outlet'       : 'zeroGradient', 
                                                'upperWall'    : (0, 0, 0),
                                                'lowerWall'    : (0, 0, 0),
                                                'frontAndBack' : 'empty'}
        foam_controlwrapper.BC[CUBA.PRESSURE]={ 'inlet'        : 'zeroGradient', 
                                                'outlet'       : 0, 
                                                'upperWall'    : 'zeroGradient',
                                                'lowerWall'    : 'zeroGradient',
                                                'frontAndBack' : 'empty'}
        
        foam_controlwrapper.run()
        

if __name__ == '__main__':
    unittest.main()

