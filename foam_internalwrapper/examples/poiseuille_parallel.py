"""Example to solve 2D poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoamII
from simphony.io.h5_cuds import H5CUDS
import os

# Only for postprocessing purposes
import matplotlib.pyplot as plt

wrapper = openfoamII.FoamControlWrapper()
CUBAExt = openfoamII.CUBAExt

name = 'poiseuille'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)
                                     
wrapper.CM_extensions[CUBAExt.NUMBER_OF_CORES] = 4

wrapper.SP[CUBA.TIME_STEP] = 1
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1000
wrapper.SP[CUBA.DENSITY] = 1.0
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'boundary2': (0.1, 0, 0),
                             'boundary1': 'zeroGradient',
                             'boundary0': (0, 0, 0),
                             'boundary3': 'empty'}
wrapper.BC[CUBA.PRESSURE] = {'boundary2': 'zeroGradient',
                             'boundary1': 0,
                             'boundary0': 'zeroGradient',
                             'boundary3': 'empty'}

mesh_file = H5CUDS.open(os.path.join(name, 'poiseuille.cuds'))
mesh_from_file = mesh_file.get_mesh(name)

print "Mesh name ", mesh_from_file.name

mesh_inside_wrapper = wrapper.add_mesh(mesh_from_file)

print "Case directory ", mesh_inside_wrapper.path

for cell in mesh_inside_wrapper.iter_cells():
    cell.data[CUBA.PRESSURE] = 1.0
    cell.data[CUBA.VELOCITY] = [0.0, 0.0, 0.0]
    mesh_inside_wrapper.update_cell(cell)

# run returns the latest time
lastTime = wrapper.run()

print "post-processing"
XYZUVW = mesh_inside_wrapper.getXYZUVW()
plt.quiver(XYZUVW[:,0],XYZUVW[:,1],XYZUVW[:,3],XYZUVW[:,4])
plt.savefig("result.png")
plt.show()

