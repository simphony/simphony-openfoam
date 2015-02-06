from simphony.core.cuba import CUBA
from foam_tokenwrapper.foam_tokenwrapper import FoamTokenWrapper
from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper

foam_controlwrapper = FoamControlWrapper()

# for OpenFoam the case -directory is given as CUBA.NAME
foam_controlwrapper.CM[CUBA.NAME]="foam_controlwrapper/examples/poiseuille"
# in the future the OpenFoam solver must be specified somewhere
#foam_controlwrapper.CM[CUBA.SOLVER]="simpleFoam"
foam_controlwrapper.SP[CUBA.TIME_STEP]=1
foam_controlwrapper.SP[CUBA.NUMBER_OF_TIME_STEPS]=1000
foam_controlwrapper.SP[CUBA.DENSITY] = 1.0
foam_controlwrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0

# this is just an example. It is not enough for general setting of BC's
foam_controlwrapper.BC[CUBA.VELOCITY]={'inlet'        : (0.1, 0, 0), 
                                       'outlet'       : 'zeroGradient', 
                                       'sides'        : (0, 0, 0),
                                       'frontAndBack' : 'empty'}
foam_controlwrapper.BC[CUBA.PRESSURE]={'inlet'        : 'zeroGradient', 
                                       'outlet'       : 0, 
                                       'sides'        : 'zeroGradient',
                                       'frontAndBack' : 'empty'}

# run returns the latest time (can be different from CUBA.TIME_STEP*CUBA.NUMBEROF_TIME_STEPS if convergence criteria is fullfilled)
lastTime = foam_controlwrapper.run()
              
foam_tokenwrapper = FoamTokenWrapper()
# SimPhoNy mesh and OpenFoam mesh uuid mapping is returned
(simPhoNyMesh,foamMeshMap) = foam_tokenwrapper.get_mesh(foam_controlwrapper.CM[CUBA.NAME],lastTime)
# uuid mapping is used to associate variable values to points
foam_tokenwrapper.get_point_data(simPhoNyMesh, foamMeshMap, CUBA.PRESSURE)


maxValue = -1.0e10
for point in simPhoNyMesh.iter_points():
    value = point.data.get(CUBA.PRESSURE)
    maxValue = max(maxValue, value)
    
print "Maximum pressure ",maxValue

