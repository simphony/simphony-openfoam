"""Example to solve 2D poiseuille flow

"""

from simphony.core.cuba import CUBA
from simphony.engine import openfoam_file_io
from simphony.io.h5_cuds import H5CUDS
import os

wrapper = openfoam_file_io.FoamControlWrapper()
CUBAExt = openfoam_file_io.CUBAExt

name = 'poiseuille'

wrapper.CM[CUBA.NAME] = name

wrapper.CM_extensions[CUBAExt.GE] = (CUBAExt.INCOMPRESSIBLE,
                                     CUBAExt.LAMINAR_MODEL)
wrapper.SP[CUBA.TIME_STEP] = 1
wrapper.SP[CUBA.NUMBER_OF_TIME_STEPS] = 1000
wrapper.SP[CUBA.DENSITY] = 1.0
wrapper.SP[CUBA.DYNAMIC_VISCOSITY] = 1.0

# this is just an example. It is not enough for general setting of BC's
wrapper.BC[CUBA.VELOCITY] = {'boundary0': (0.1, 0, 0),
                             'boundary1': 'zeroGradient',
                             'boundary2': (0, 0, 0),
                             'boundary3': 'empty'}
wrapper.BC[CUBA.PRESSURE] = {'boundary0': 'zeroGradient',
                             'boundary1': 0,
                             'boundary2': 'zeroGradient',
                             'boundary3': 'empty'}

mesh_file = H5CUDS.open(os.path.join(name, 'poiseuille.cuds'))
mesh_from_file = mesh_file.get_mesh(name)

print "Mesh name ", mesh_from_file.name

mesh_inside_wrapper = wrapper.add_mesh(mesh_from_file)

print "Case directory ", mesh_inside_wrapper.path

wrapper.run()
