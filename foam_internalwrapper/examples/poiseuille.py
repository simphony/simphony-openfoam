"""Example to solve 2D poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_internal
from simphony.engine import openfoam_file_io

from mayavi.scripts import mayavi2

import tempfile

wrapper = openfoam_internal.Wrapper()
CUBAExt = openfoam_internal.CUBAExt

name = 'poiseuille'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)

wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 1

wrapper.SP[CUBA.TIME_STEP] = 0.001
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1000
wrapper.SP[CUBA.DENSITY] = 1.0
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 0.001

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
                 (0.0, 0.0, 0.1e-3), (20.0e-3, 0.0, 0.1e-3),
                 (20.0e-3, 1.0e-3, 0.1e-3), (0.0, 1.0e-3, 0.1e-3)]

# elements in x -direction
nex = 8
# elements in y -direction
ney = 4
openfoam_file_io.create_quad_mesh(tempfile.mkdtemp(), name, wrapper,
                                  corner_points, nex, ney, 1)
mesh_inside_wrapper = wrapper.get_dataset(name)

# run returns the latest time
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
