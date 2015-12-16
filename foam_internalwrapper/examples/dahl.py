"""Example to solve 2D poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_internal
from simphony.engine import openfoam_file_io
import dahl_mesh

wrapper = openfoam_internal.FoamInternalWrapper()
CUBAExt = openfoam_internal.CUBAExt

import os
path = os.path.abspath(os.curdir)
name = 'dahl'

wrapper.CM[CUBA.NAME] = name
                                     
wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL,
                                     CUBAExt.MIXTURE)

wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 1

wrapper.SP[CUBA.TIME_STEP] = 0.1
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 5000

wrapper.SP_extensions[CUBAExt.PHASE_LIST] = ('sludge', 'water')
wrapper.SP[CUBA.DENSITY] = {'sludge': 1996.0, 'water': 996}
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = {'sludge': 'BinghamPlastic', 'water': 0.0017}

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'boundary0': (0, 0, 0),
                             'boundary1': 'zeroGradient', #pressureInletOutletVelocity uniform (0 0 0);
                             'boundary2': (0, 0, 0),
                             'boundary3': (0, 0, 0),
                             'boundary4': 'slip',
                             'boundary5': 'empty'}
wrapper.BC[CUBA.CONCENTRATION] = {'boundary0': 'fixedFluxPressure',
                                  'boundary1': 0,
                                  'boundary2': 'fixedFluxPressure',
                                  'boundary3': 'fixedFluxPressure',
                                  'boundary4': 'fixedFluxPressure',
                                  'boundary5': 'empty'}
wrapper.BC[CUBA.VOLUME_FRACTION] = {'boundary0': 0.001,
                                    'boundary1': ('inletOutlet', 0.001),
                                    'boundary2': 'zeroGradient',
                                    'boundary3': 'zeroGradient',
                                    'boundary4': 'zeroGradient',
                                    'boundary5': 'empty'}
                                    
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 10

#creating the mesh
file_name = 'blockMeshDict'
case = os.path.join(path, name)
templateName = 'simpleFoam'
openfoam_file_io.foam_files.write_default_files(case, templateName, '0', True)

full_name = os.path.join(os.path.join(
        os.path.join(case, 'constant'), 'polyMesh'), file_name)
with open(full_name, 'w') as f:
    f.write(dahl_mesh.blockMeshDict)

ncores = 1
solver = 'blockMesh'
runner = openfoam_file_io.foam_runner.FoamRunner(solver, case, ncores)
runner.run()

foam_mesh = openfoam_file_io.read_foammesh(name, path)

# add mesh to engine
wrapper.add_dataset(foam_mesh)


mesh_inside_wrapper = wrapper.get_dataset(name)

updated_cells = []
for cell in mesh_inside_wrapper.iter_cells():
    cell.data[CUBA.VOLUME_FRACTION] = 0.001
    cell.data[CUBA.CONCENTRATION] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0191, 0.0, 0.0]
    updated_cells.append(cell)

mesh_inside_wrapper.update_cells(updated_cells)

# run returns the latest time
wrapper.run()
