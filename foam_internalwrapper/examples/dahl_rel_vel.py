"""Example to solve mixture model
"""

from simphony.core.cuba import CUBA

from simphony.engine import openfoam_internal
from simphony.engine import openfoam_file_io

import dahl_mesh
import tempfile


wrapper = openfoam_internal.Wrapper()
CUBAExt = openfoam_internal.CUBAExt

name = 'dahl'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL,
                                     CUBAExt.MIXTURE_MODEL)

wrapper.SP[CUBA.TIME_STEP] = 0.1
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1

wrapper.SP_extensions[CUBAExt.PHASE_LIST] = ('sludge', 'water')
wrapper.SP[CUBA.DENSITY] = {'sludge': 1900.0, 'water': 1000.0}
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = {'sludge': 0.01, 'water': 1e-3}
wrapper.SP_extensions[CUBAExt.VISCOSITY_MODEL] =\
    {'sludge': 'BinghamPlastic', 'water': 'Newtonian'}
wrapper.SP_extensions[CUBAExt.VISCOSITY_MODEL_COEFFS] =\
    {'BinghamPlastic': {'coeff': 0.00023143, 'exponent': 179.26,
                        'BinghamCoeff': 0.0005966, 'BinghamExponent': 1050.8,
                        'BinghamOffset': 0, 'muMax': 10}}
wrapper.SP_extensions[CUBAExt.STRESS_MODEL] = 'standard'
wrapper.SP_extensions[CUBAExt.RELATIVE_VELOCITY_MODEL] = 'fromMesoscale'
wrapper.SP_extensions[CUBAExt.EXTERNAL_BODY_FORCE_MODEL] = 'gravitation'
wrapper.SP_extensions[CUBAExt.EXTERNAL_BODY_FORCE_MODEL_COEFFS] =\
    {'g': (0.0, -9.81, 0.0)}
wrapper.BC[CUBA.VELOCITY] = {'boundary0': ('fixedValue', (0.0191, 0, 0)),
                             'boundary1': ('pressureIOVelocity', (0, 0, 0)),
                             'boundary2': ('fixedValue', (0, 0, 0)),
                             'boundary3': ('fixedValue', (0, 0, 0)),
                             'boundary4': 'slip',
                             'boundary5': 'empty'}
# CUBA.CONCENTRATION is used for dynamic pressure while not in CUBA keys
wrapper.BC[CUBA.CONCENTRATION] = {'boundary0': 'fixedFluxPressure',
                                  'boundary1': ('fixedValue', 0),
                                  'boundary2': 'fixedFluxPressure',
                                  'boundary3': 'fixedFluxPressure',
                                  'boundary4': 'fixedFluxPressure',
                                  'boundary5': 'empty'}

wrapper.BC[CUBA.VOLUME_FRACTION] = {'boundary0': ('fixedValue', 0.001),
                                    'boundary1': ('inletOutlet', 0.001),
                                    'boundary2': 'zeroGradient',
                                    'boundary3': 'zeroGradient',
                                    'boundary4': 'zeroGradient',
                                    'boundary5': 'empty'}

# create mesh
openfoam_file_io.create_block_mesh(tempfile.mkdtemp(), name, wrapper,
                                   dahl_mesh.blockMeshDict)

mesh_inside_wrapper = wrapper.get_dataset(name)

updated_cells = []
for cell in mesh_inside_wrapper.iter_cells():
    cell.data[CUBA.VOLUME_FRACTION] = 0.001
    cell.data[CUBA.CONCENTRATION] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0191, 0.0, 0.0]
    # CUBA.ANGULAR_VELOCITY is used for relative velocity
    cell.data[CUBA.ANGULAR_VELOCITY] = [0.0, 0.0, 0.0]
    updated_cells.append(cell)

mesh_inside_wrapper.update_cells(updated_cells)

V0 = [0.0, -0.002, 0.0]
a = 285.0

number_of_outer_timesteps = 6400

for time_i in range(number_of_outer_timesteps):
    # solve macroscopic scale
    print "Solve cfd"
    wrapper.run()
    print "Time: ", mesh_inside_wrapper._time
    print "Mesoscale as analytic coupling"
    print " Update relative velocity"

    updated_cells = []
    for cell in mesh_inside_wrapper.iter_cells():
        alphad = cell.data[CUBA.VOLUME_FRACTION]
        vdj = [V*pow(10.0, -a*max(alphad, 0.0)) for V in V0]
        cell.data[CUBA.ANGULAR_VELOCITY] = vdj
        updated_cells.append(cell)

    mesh_inside_wrapper.update_cells(updated_cells)
