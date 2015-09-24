"""Example to convert mesh from OpenFoam's format to H5CUDS

"""

from foam_controlwrapper.foam_controlwrapper import FoamControlWrapper
from foam_controlwrapper.foam_controlwrapper import read_foammesh
from simphony.io.h5_cuds import H5CUDS

foam_controlwrapper = FoamControlWrapper()
name = 'poiseuille_vof'
#name = 'poiseuille'
#name = 'simplemesh'
path = '.'
mesh_inside_wrapper = read_foammesh(name, path)
mesh_file = H5CUDS.open("poiseuille_vof.cuds")
#mesh_file = H5CUDS.open("poiseuille.cuds")
#mesh_file = H5CUDS.open("simplemesh.cuds")
mesh_file.add_dataset(mesh_inside_wrapper)
