""" Utility functions for file IO

"""
import simphonyfoaminterface as foamface

from .foam_mesh import FoamMesh


def read_foammesh(name, path):
    """Read mesh from OpenFoam case files.

    Parameters
    ----------
    name : str
    name to give to mesh
    path : str
    case directory

    Raises
    ------
    Exception if some mesh from mesh names list not found

    """

    foamface.init_IO(name, path)
    foamface.readMesh(name)
    nPoints = foamface.getPointCount(name)
    nCells = foamface.getCellCount(name)
    nFaces = foamface.getFaceCount(name)
    nEdges = 0

    foamMesh = FoamMesh(name)
    foamMesh.generate_uuidmapping(nPoints, nEdges, nFaces, nCells)
    return foamMesh
