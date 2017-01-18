"""Example to solve 2D vortex shedding flow

"""
from mayavi.scripts import mayavi2

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_file_io
import vortex_shedding_mesh

import tempfile
import math

wrapper = openfoam_file_io.Wrapper()
CUBAExt = openfoam_file_io.CUBAExt

name = 'vortex_shedding'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)
wrapper.SP[CUBA.TIME_STEP] = 0.025
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 6000
wrapper.SP[CUBA.DENSITY] = 1.0
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 2.0e-2

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'boundary0': ('fixedValue', (1.0, 0, 0)),
                             'boundary1': 'zeroGradient',
                             'boundary2': 'slip',
                             'boundary3': 'slip',
                             'boundary4': ('fixedValue', (0, 0, 0)),
                             'boundary5': 'empty'}
wrapper.BC[CUBA.PRESSURE] = {'boundary0': 'zeroGradient',
                             'boundary1': ('fixedValue', 0),
                             'boundary2': 'zeroGradient',
                             'boundary3': 'zeroGradient',
                             'boundary4': 'zeroGradient',
                             'boundary5': 'empty'}

# create mesh
openfoam_file_io.create_block_mesh(tempfile.mkdtemp(), name, wrapper,
                                   vortex_shedding_mesh.blockMeshDict)

mesh_inside_wrapper = wrapper.get_dataset(name)

print "Case directory ", mesh_inside_wrapper.path

wrapper.run()

# compute velocity magnitude (volume fraction here)
updated_cells = []
for cell in mesh_inside_wrapper.iter(item_type=CUBA.CELL):
    velo = cell.data[CUBA.VELOCITY]
    cell.data[CUBA.VOLUME_FRACTION] = math.sqrt(sum(velo[i]*velo[i]
                                                    for i in range(3)))
    updated_cells.append(cell)
mesh_inside_wrapper.update(updated_cells)


@mayavi2.standalone
def view(dataset):
    from mayavi.modules.surface import Surface
    from simphony_mayavi.sources.api import CUDSSource

    mayavi.new_scene()  # noqa
    src = CUDSSource(cuds=dataset)
    mayavi.add_source(src)  # noqa
    s = Surface()
    mayavi.add_module(s)  # noqa


if __name__ == '__main__':
    view(mesh_inside_wrapper)
