"""Example to convert mesh from OpenFoam's format to SimPhony

"""

from simphony.engine import openfoam_file_io

wrapper = openfoam_file_io.Wrapper()

name = 'poiseuille'
path = '.'
mesh_inside_wrapper = openfoam_file_io.read_foammesh(name, path)
