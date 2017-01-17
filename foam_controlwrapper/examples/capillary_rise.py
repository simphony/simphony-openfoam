"""Example to solve 2D 2 -phase poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_file_io

from mayavi.scripts import mayavi2

import tube_mesh

import tempfile

wrapper = openfoam_file_io.Wrapper()
CUBAExt = openfoam_file_io.CUBAExt


name = 'capillary_rise'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL,
                                     CUBAExt.VOF_MODEL)
wrapper.SP[CUBA.TIME_STEP] = 1.0e-5
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 50000
wrapper.SP_extensions[CUBAExt.MAX_COURANT_NUMBER] = 0.2

wrapper.SP_extensions[CUBAExt.PHASE_LIST] = ('water', 'air')
wrapper.SP[CUBA.DENSITY] = {'water': 1000.0, 'air': 1.0}
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = {'water': 0.001, 'air': 1.48e-5}
wrapper.SP_extensions[CUBAExt.SURFACE_TENSION] = {('water', 'air'): 70.7106e-3}
wrapper.SP_extensions[CUBAExt.EXTERNAL_BODY_FORCE_MODEL] = 'gravitation'
wrapper.SP_extensions[CUBAExt.EXTERNAL_BODY_FORCE_MODEL_COEFFS] =\
    {'g': (0.0, -9.81, 0.0)}

wrapper.BC[CUBA.VELOCITY] = {'inlet': ('pressureIOVelocity', (0, 0, 0)),
                             'atmosphere': ('pressureIOVelocity', (0, 0, 0)),
                             'walls': ('fixedValue', (0, 0, 0)),
                             'frontAndBack': 'empty'}
wrapper.BC[CUBA.DYNAMIC_PRESSURE] = {'inlet': ('fixedValue', 0),
                                     'atmosphere': ('fixedValue', 0),
                                     'walls': 'fixedFluxPressure',
                                     'frontAndBack': 'empty'}
wrapper.BC[CUBA.VOLUME_FRACTION] = {'inlet': ('inletOutlet', 1),
                                    'atmosphere': 'zeroGradient',
                                    'walls': ('wettingAngle', 45),
                                    'frontAndBack': 'empty'}


# create mesh
openfoam_file_io.create_block_mesh(tempfile.mkdtemp(), name, wrapper,
                                   tube_mesh.blockMeshDict)

mesh_inside_wrapper = wrapper.get_dataset(name)

# initial state. In VOF only one velocity and pressure field

print "Case directory ", mesh_inside_wrapper.path

updated_cells = []
for cell in mesh_inside_wrapper.iter(item_type=CUBA.CELL):
    ymid = sum(mesh_inside_wrapper.get_point(puid).coordinates[1]
               for puid in cell.points)
    ymid /= sum(1.0 for _ in cell.points)
    if ymid < 8e-3:
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
