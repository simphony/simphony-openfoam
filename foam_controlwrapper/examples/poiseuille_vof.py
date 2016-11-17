"""Example to solve 2D 2 -phase poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_file_io

from mayavi.scripts import mayavi2

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
wrapper.SP_extensions[CUBAExt.MAX_COURANT_NUMBER] = 0.5

wrapper.SP_extensions[CUBAExt.PHASE_LIST] = ('water', 'air')
wrapper.SP[CUBA.DENSITY] = {'water': 1000.0, 'air': 1.0}
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = {'water': 0.001, 'air': 1.8e-5}
wrapper.SP_extensions[CUBAExt.SURFACE_TENSION] = {('water', 'air'): 72.86e-3}

wrapper.BC[CUBA.VELOCITY] = {'inlet': 'zeroGradient',
                             'outlet': 'zeroGradient',
                             'walls': ('fixedValue', (0, 0, 0)),
                             'frontAndBack': 'empty'}
wrapper.BC[CUBA.DYNAMIC_PRESSURE] = {'inlet': ('fixedValue', 1),
                                     'outlet': ('fixedValue', 0),
                                     'walls': 'fixedFluxPressure',
                                     'frontAndBack': 'empty'}
wrapper.BC[CUBA.VOLUME_FRACTION] = {'inlet': ('fixedValue', 1),
                                    'outlet': 'zeroGradient',
                                    'walls': 'zeroGradient',
                                    'frontAndBack': 'empty'}

len_x = 20.0e-3
len_y = 1.0e-3
len_z = 0.1e-3

corner_points = [(0.0, 0.0, 0.0), (len_x, 0.0, 0.0),
                 (len_x, len_y, 0.0), (0.0, len_y, 0.0),
                 (0.0, 0.0, len_z), (len_x, 0.0, len_z),
                 (len_x, len_y, len_z), (0.0, len_y, len_z)]

# elements in x -direction
nex = 100
# elements in y -direction
ney = 10
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
    if xmid < len_x/3.:
        cell.data[CUBA.VOLUME_FRACTION] = 1.0
    else:
        cell.data[CUBA.VOLUME_FRACTION] = 0.0

    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0, 0.0, 0.0]

    updated_cells.append(cell)

mesh_inside_wrapper.update_cells(updated_cells)

wrapper.run()


@mayavi2.standalone
def view():
    from mayavi.modules.surface import Surface
    from simphony_mayavi.sources.api import CUDSSource

    mayavi.new_scene()  # noqa
    src = CUDSSource(cuds=mesh_inside_wrapper)
    mayavi.add_source(src)  # noqa
    s = Surface()
    mayavi.add_module(s)  # noqa


if __name__ == '__main__':
    view()
