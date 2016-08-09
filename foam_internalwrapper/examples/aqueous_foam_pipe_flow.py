"""Example to solve 3D aqueous foam pipe flow using rheological
Herschel-Bulkley power law for bulk and wall shear stress dependent
slip velocity law for wall layer

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_file_io
from simphony.engine import openfoam_internal

from mayavi.scripts import mayavi2

import pipe_mesh
import tempfile

wrapper = openfoam_internal.Wrapper()
CUBAExt = openfoam_internal.CUBAExt

name = 'aqueous_foam'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)

wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 1

wrapper.SP[CUBA.TIME_STEP] = 0.0001
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 3000
wrapper.SP[CUBA.DENSITY] = 250.0
wrapper.SP_extensions[CUBAExt.VISCOSITY_MODEL] = 'HerschelBulkley'
wrapper.SP_extensions[CUBAExt.VISCOSITY_MODEL_COEFFS] =\
    {'HerschelBulkley': {'nu0': 0.01748,
                         'tau0': 0.0148,
                         'k': 0.00268,
                         'n': 0.5}}


wrapper.BC[CUBA.VELOCITY] = {'inlet': ('fixedValue', (0, 0, 0.53)),
                             'outlet': 'zeroGradient',
                             'walls': ('shearStressPowerLawSlipVelocity',
                                       {'rho': 250.0,
                                        'beta': 3.1e-3,
                                        'n': 1.16})}

wrapper.BC[CUBA.PRESSURE] = {'inlet': 'zeroGradient',
                             'outlet': ('fixedValue', 0),
                             'walls': 'zeroGradient'}


# create mesh
openfoam_file_io.create_block_mesh(tempfile.mkdtemp(), name, wrapper,
                                   pipe_mesh.blockMeshDict)

mesh_inside_wrapper = wrapper.get_dataset(name)


# run returns the latest time
wrapper.run()

average_pressure = 0.0
for cell in mesh_inside_wrapper.get_boundary_cells('inlet'):
    average_pressure += cell.data[CUBA.PRESSURE]

average_pressure /= len(mesh_inside_wrapper._boundaries['inlet'])

print "Average pressure on inlet: ", average_pressure


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
