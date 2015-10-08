""" Utility functions for blockMesh

"""

import os
import tempfile

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

from .foam_files import write_default_files
from .foam_templates import blockMeshDict
from .foam_runner import FoamRunner
from .io_utils import read_foammesh


def create_quad_mesh(name, mesh_engine, corner_points,
                     nex, ney, nez):
    """ create and add mesh to engine to wrapper

    Parameters
    ----------
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
    path = os.path.join(tempfile.mkdtemp())
    case = os.path.join(path, name)
    templateName = 'simpleFoam'
    write_default_files(case, templateName, '0', True)
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

    ncores = 1
    solver = 'blockMesh'
    runner = FoamRunner(solver, case, ncores)
    runner.run()

    foam_mesh = read_foammesh(name, path)

    # add mesh to engine
    mesh_engine.add_dataset(foam_mesh)
