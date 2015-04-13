"""Example to convert mesh from OpenFoam's format to H5CUDS

"""

from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from simphony.io.h5_cuds import H5CUDS

foam_controlwrapper = FoamControlWrapper()
name = 'poiseuille'
path = '.'
mesh_inside_wrapper = foam_controlwrapper.read_foammesh(name, path)
mesh_file = H5CUDS.open("poiseuille.cuds")
print 'Adding mesh ', mesh_inside_wrapper.name
mesh_file.add_mesh(mesh_inside_wrapper)
