"""Example to solve 2D poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_file_io
import tempfile

wrapper = openfoam_file_io.Wrapper()
CUBAExt = openfoam_file_io.CUBAExt

name = 'poiseuille'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)
wrapper.SP[CUBA.TIME_STEP] = 0.001
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1000
wrapper.SP[CUBA.DENSITY] = 1.0
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'inlet': ('fixedValue', (0.1, 0, 0)),
                             'outlet': 'zeroGradient',
                             'walls': ('fixedValue', (0, 0, 0)),
                             'frontAndBack': 'empty'}
wrapper.BC[CUBA.PRESSURE] = {'inlet': 'zeroGradient',
                             'outlet': ('fixedValue', 0),
                             'walls': 'zeroGradient',
                             'frontAndBack': 'empty'}

corner_points = [(0.0, 0.0, 0.0), (20.0e-3, 0.0, 0.0),
                 (20.0e-3, 1.0e-3, 0.0), (0.0, 1.0e-3, 0.0),
                 (0.0, 0.0, 0.1), (20.0e-3, 0.0, 0.1),
                 (20.0e-3, 1.0e-3, 0.1), (0.0, 1.0e-3, 0.1)]
# elements in x -direction
nex = 8
# elements in y -direction
ney = 4
openfoam_file_io.create_quad_mesh(tempfile.mkdtemp(), name, wrapper,
                                  corner_points, nex, ney, 1)

mesh_inside_wrapper = wrapper.get_dataset(name)

print "Case directory ", mesh_inside_wrapper.path

wrapper.run()
