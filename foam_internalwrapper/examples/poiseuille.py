"""Example to solve 2D poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_internal
from simphony.engine import openfoam_file_io

wrapper = openfoam_internal.FoamInternalWrapper()
CUBAExt = openfoam_internal.CUBAExt

import os
path = os.path.abspath(os.curdir)
name = 'poiseuille'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)

wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 1

wrapper.SP[CUBA.TIME_STEP] = 0.001
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 100

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'boundary0': (0, 0, 0),
                             'boundary2': 'zeroGradient',
                             'boundary1': (0.1, 0, 0),
                             'boundary3': 'empty'}
wrapper.BC[CUBA.PRESSURE] = {'boundary0': 'zeroGradient',
                             'boundary2': 0,
                             'boundary1': 'zeroGradient',
                             'boundary3': 'empty'}

corner_points = [(0.0, 0.0, 0.0), (20.0e-3, 0.0, 0.0),
                 (20.0e-3, 1.0e-3, 0.0), (0.0, 1.0e-3, 0.0),
                 (0.0, 0.0, 0.1), (20.0e-3, 0.0, 0.1),
                 (20.0e-3, 1.0e-3, 0.1), (0.0, 1.0e-3, 0.1)]

# elements in x -direction
nex = 10
# elements in y -direction
ney = 4
openfoam_file_io.create_quad_mesh(path, name, wrapper, corner_points,
                                  nex, ney, 1)

mesh_inside_wrapper = wrapper.get_dataset(name)


# run returns the latest time
wrapper.run()
