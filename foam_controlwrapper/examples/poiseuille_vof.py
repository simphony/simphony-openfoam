"""Example to solve 2D 2 -phase poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_file_io
import tempfile

wrapper = openfoam_file_io.Wrapper()
CUBAExt = openfoam_file_io.CUBAExt


name = 'poiseuille_vof'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL,
                                     CUBAExt.VOF_MODEL)
wrapper.SP[CUBA.TIME_STEP] = 0.001
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 100

wrapper.SP_extensions[CUBAExt.PHASE_LIST] = ('water', 'air')
wrapper.SP[CUBA.DENSITY] = {'water': 1000.0, 'air': 1.0}
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = {'water': 0.001, 'air': 1.8e-5}
wrapper.SP_extensions[CUBAExt.SURFACE_TENSION] = {('water', 'air'): 72.86e-3}

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'inlet': ('fixedValue', (0, 0, 0)),
                             'outlet': ('fixedValue', (0.01e-3, 0, 0)),
                             'walls': 'zeroGradient',
                             'frontAndBack': 'empty'}
wrapper.BC[CUBA.DYNAMIC_PRESSURE] = {'inlet': 'zeroGradient',
                                     'outlet': 'zeroGradient',
                                     'walls': ('fixedValue', 0),
                                     'frontAndBack': 'empty'}
wrapper.BC[CUBA.VOLUME_FRACTION] = {'inlet': 'zeroGradient',
                                    'outlet': ('fixedValue', 1),
                                    'walls': 'zeroGradient',
                                    'frontAndBack': 'empty'}

corner_points = [(0.0, 0.0, 0.0), (20.0e-3, 0.0, 0.0),
                 (20.0e-3, 1.0e-3, 0.0), (0.0, 1.0e-3, 0.0),
                 (0.0, 0.0, 0.1e-3), (20.0e-3, 0.0, 0.1e-3),
                 (20.0e-3, 1.0e-3, 0.1e-3), (0.0, 1.0e-3, 0.1e-3)]

# elements in x -direction
nex = 40
# elements in y -direction
ney = 4
openfoam_file_io.create_quad_mesh(tempfile.mkdtemp(), name, wrapper,
                                  corner_points, nex, ney, 1)

mesh_inside_wrapper = wrapper.get_dataset(name)

# initial state. In VOF only one velocity and pressure field

print "Case directory ", mesh_inside_wrapper.path

updated_cells = []
for cell in mesh_inside_wrapper.iter_cells():
    xmid = sum(mesh_inside_wrapper.get_point(puid).coordinates[0]
               for puid in cell.points)
    xmid /= sum(1.0 for _ in cell.points)
    if xmid < 0.02/3.:
        cell.data[CUBA.VOLUME_FRACTION] = 1.0
    else:
        cell.data[CUBA.VOLUME_FRACTION] = 0.0

    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0, 0.0, 0.0]

    updated_cells.append(cell)

mesh_inside_wrapper.update_cells(updated_cells)

wrapper.run()
