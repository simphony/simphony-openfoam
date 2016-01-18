""" Utility functions for blockMesh

"""

import os
import shutil

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

from .foam_files import write_default_files, remove_parser_files
from .foam_templates import blockMeshDict
from .foam_runner import FoamRunner
from .io_utils import read_foammesh


def create_quad_mesh(path, name, mesh_engine, corner_points,
                     nex, ney, nez):
    """ create and add mesh to engine

    Parameters
    ----------
    path : str
        path to mesh parent directory
    name : str
        name of mesh
    mesh_engine : ABCModelingEngine
        Mesh engine
    corner_points : list
        list of 8 [x,y,z] corner points
    nex : int
        number of elements in x -direction
    ney : int
        number of elements in y -direction
    nez : int
        number of elements in z -direction

    """
    file_name = 'blockMeshDict'
    case = os.path.join(path, name)
    templateName = 'simpleFoam'
    write_default_files(case, templateName, '0', False)
    full_name = os.path.join(os.path.join(
        os.path.join(case, 'constant'), 'polyMesh'), file_name)
    with open(full_name, 'w') as f:
        f.write(blockMeshDict)

    blockMesh = ParsedParameterFile(full_name)

    for i in range(8):
        corner_points[i] = str(corner_points[i]).replace(',', ' ')

    blockMesh["vertices"] = corner_points

    blockLines = [""]
    blockLines[0] = 'hex (0 1 2 3 4 5 6 7) (%i %i %i) simpleGrading (1 1 1)'\
        % (nex, ney, nez)
    blockMesh["blocks"] = blockLines

    blockMesh.writeFile()
    # remove PyFoam parser files
    remove_parser_files(os.getcwd())

    # this to overcome bug in blockMesh case attribute
    # blockMesh searches blockMeshDict -file from doubled case directory
    # copy file to that directory
    blockMesh_file_name =\
        os.path.join(case, 'constant', 'polyMesh', file_name)
    blockMesh_file_name = case + os.sep + blockMesh_file_name
    if not os.path.exists(os.path.dirname(blockMesh_file_name)):
        os.makedirs(os.path.dirname(blockMesh_file_name))
    shutil.copy(full_name, blockMesh_file_name)

    ncores = 1
    solver = 'blockMesh'
    runner = FoamRunner(solver, case, ncores)
    runner.run()

    foam_mesh = read_foammesh(name, path)

    # add mesh to engine
    mesh_engine.add_dataset(foam_mesh)


def create_block_mesh(path, name, mesh_engine, block_mesh_dict):
    """ create and add mesh to engine

    Parameters
    ----------
    path : str
        path to mesh parent directory
    name : str
        name of mesh
    mesh_engine : ABCModelingEngine
        Mesh engine
    block_mesh_dict : str
        blockMeshDict -file as a string

    """
    file_name = 'blockMeshDict'
    case = os.path.join(path, name)
    templateName = 'simpleFoam'
    write_default_files(case, templateName, '0', False)
    full_name = os.path.join(os.path.join(
        os.path.join(case, 'constant'), 'polyMesh'), file_name)
    with open(full_name, 'w') as f:
        f.write(block_mesh_dict)

    # this to overcome bug in blockMesh case attribute
    # blockMesh searches blockMeshDict -file from doubled case directory
    # copy file to that directory
    blockMesh_file_name =\
        os.path.join(case, 'constant', 'polyMesh', file_name)
    blockMesh_file_name = case + os.sep + blockMesh_file_name

    if not os.path.exists(os.path.dirname(blockMesh_file_name)):
        os.makedirs(os.path.dirname(blockMesh_file_name))
    shutil.copy(full_name, blockMesh_file_name)

    ncores = 1
    solver = 'blockMesh'
    runner = FoamRunner(solver, case, ncores)
    runner.run()

    foam_mesh = read_foammesh(name, path)

    # add mesh to engine
    mesh_engine.add_dataset(foam_mesh)
