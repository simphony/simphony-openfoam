from simphony.core.cuba import CUBA
from simphony.engine import openfoam

wrapper = openfoam.FoamControlWrapper()
CUBAExt = openfoam.CUBAExt

# for OpenFoam the case -directory is given as CUBA.NAME
wrapper.CM[CUBA.NAME] = "foam_controlwrapper/examples/poiseuille"

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)
wrapper.SP[CUBA.TIME_STEP] = 1
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1000
wrapper.SP[CUBA.DENSITY] = 1.0
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'inlet': (0.1, 0, 0),
                             'outlet': 'zeroGradient',
                             'sides': (0, 0, 0),
                             'frontAndBack': 'empty'}
wrapper.BC[CUBA.PRESSURE] = {'inlet': 'zeroGradient',
                             'outlet': 0,
                             'sides': 'zeroGradient',
                             'frontAndBack': 'empty'}

# run returns the latest time
lastTime = wrapper.run()


print "Timesteps taken ", lastTime
