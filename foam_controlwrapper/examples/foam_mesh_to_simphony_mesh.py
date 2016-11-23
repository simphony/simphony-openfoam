"""Example to convert mesh from OpenFoam's format to SimPhony

"""

import foam_controlwrapper

name = 'poiseuille'
path = '.'
mesh_inside_wrapper = foam_controlwrapper.read_foammesh(name, path)
