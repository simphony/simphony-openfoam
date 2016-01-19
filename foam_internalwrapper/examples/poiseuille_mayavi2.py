"""Example to solve 2D poiseuille flow

"""

import os

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_internal
from simphony.engine import openfoam_file_io
from simphony.core.cuds_item import CUDSItem

wrapper = openfoam_internal.Wrapper()
CUBAExt = openfoam_internal.CUBAExt

path = os.path.abspath(os.curdir)
name = 'poiseuille'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)

wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 1

wrapper.SP[CUBA.TIME_STEP] = 0.001
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1000
wrapper.SP[CUBA.DENSITY] = 1.0
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 0.001

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'boundary0': ('fixedValue', (0.1, 0, 0)),
                             'boundary1': 'zeroGradient',
                             'boundary2': ('fixedValue', (0, 0, 0)),
                             'boundary3': 'empty'}
wrapper.BC[CUBA.PRESSURE] = {'boundary0': 'zeroGradient',
                             'boundary1': ('fixedValue', 0),
                             'boundary2': 'zeroGradient',
                             'boundary3': 'empty'}

corner_points = [(0.0, 0.0, 0.0), (20.0e-3, 0.0, 0.0),
                 (20.0e-3, 1.0e-3, 0.0), (0.0, 1.0e-3, 0.0),
                 (0.0, 0.0, 0.1e-3), (20.0e-3, 0.0, 0.1e-3),
                 (20.0e-3, 1.0e-3, 0.1e-3), (0.0, 1.0e-3, 0.1e-3)]

# elements in x -direction
nex = 50
# elements in y -direction
ney = 6
openfoam_file_io.create_quad_mesh(path, name, wrapper, corner_points,
                                  nex, ney, 1)
mesh_inside_wrapper = wrapper.get_dataset(name)

#initial values
for cell in mesh_inside_wrapper.iter_cells():
    cell.data[CUBA.VELOCITY] = [0.1, 0, 0]

wrapper.run()

avg_velo = 0.0
for cell in mesh_inside_wrapper.iter_cells():
    avg_velo += cell.data[CUBA.VELOCITY][0]

print "Average velocity ",avg_velo/mesh_inside_wrapper.count_of(CUDSItem.CELL)
# Now view the data.
from mayavi.scripts import mayavi2
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

