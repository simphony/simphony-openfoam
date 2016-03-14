""" Utility functions for blockMesh

"""

import os
import shutil

from foam_internalwrapper.foam_dicts import (dictionaryMaps, parse_map,
                                             write_dictionary,
                                             create_directories)
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

    case = os.path.join(path, name)

    create_directories(case)

    solver = 'blockMesh'
    dictionary_name = 'controlDict'
    full_path = os.path.join(case, 'system')

    map_content = dictionaryMaps[solver]

    controlDict = parse_map(map_content[dictionary_name])
    write_dictionary(full_path, dictionary_name, controlDict)

    dictionary_name = 'fvSchemes'
    controlDict = parse_map(map_content[dictionary_name])
    write_dictionary(full_path, dictionary_name, controlDict)

    dictionary_name = 'fvSolution'
    controlDict = parse_map(map_content[dictionary_name])
    write_dictionary(full_path, dictionary_name, controlDict)

    dictionary_name = 'blockMeshDict'
    control = map_content[dictionary_name]

    dictionary_path = os.path.join('constant', 'polyMesh')
    full_name = os.path.join(os.path.join(case, dictionary_path),
                             dictionary_name)

    corners = "("
    for i in range(8):
        corners += str(corner_points[i]).replace(',', ' ')
    corners += ")"

    control['vertices'] = corners

    control['blocks'] =\
        '(hex (0 1 2 3 4 5 6 7) (%i %i %i) simpleGrading (1 1 1))'\
        % (nex, ney, nez)

    control['boundary'] =\
        """
    (
    walls
    {
        type patch;
        faces
        (
            (3 7 6 2)
            (1 5 4 0)
        );
    }
    inlet
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (2 6 5 1)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }
)
"""

    blockMeshDict = parse_map(control)
    write_dictionary(os.path.join(case, dictionary_path),
                     dictionary_name, blockMeshDict)

    # this to overcome bug in blockMesh case attribute
    # blockMesh searches blockMeshDict -file from doubled case directory
    # copy file to that directory
    blockMesh_file_name = case + os.sep + full_name

    if not os.path.exists(os.path.dirname(blockMesh_file_name)):
        os.makedirs(os.path.dirname(blockMesh_file_name))
    shutil.copy(full_name, blockMesh_file_name)

    ncores = 1
    solver = 'blockMesh'
    runner = FoamRunner(solver, name, case, ncores)
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

    case = os.path.join(path, name)

    create_directories(case)

    solver = 'blockMesh'
    dictionary_name = 'controlDict'
    full_path = os.path.join(case, 'system')

    map_content = dictionaryMaps[solver]

    controlDict = parse_map(map_content[dictionary_name])
    write_dictionary(full_path, dictionary_name, controlDict)

    dictionary_name = 'fvSchemes'
    controlDict = parse_map(map_content[dictionary_name])
    write_dictionary(full_path, dictionary_name, controlDict)

    dictionary_name = 'fvSolution'
    controlDict = parse_map(map_content[dictionary_name])
    write_dictionary(full_path, dictionary_name, controlDict)

    dictionary_name = 'blockMeshDict'
    dictionary_path = os.path.join('constant', 'polyMesh')
    full_name = os.path.join(os.path.join(case, dictionary_path),
                             dictionary_name)

    write_dictionary(os.path.join(case, dictionary_path),
                     dictionary_name, block_mesh_dict)

    # this to overcome bug in blockMesh case attribute
    # blockMesh searches blockMeshDict -file from doubled case directory
    # copy file to that directory
    blockMesh_file_name = case + os.sep + full_name

    if not os.path.exists(os.path.dirname(blockMesh_file_name)):
        os.makedirs(os.path.dirname(blockMesh_file_name))
    shutil.copy(full_name, blockMesh_file_name)

    ncores = 1
    solver = 'blockMesh'
    runner = FoamRunner(solver, name, case, ncores)
    runner.run()

    foam_mesh = read_foammesh(name, path)

    # add mesh to engine
    mesh_engine.add_dataset(foam_mesh)
