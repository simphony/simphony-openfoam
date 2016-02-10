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
    patchNames = foamface.getBoundaryPatchNames(name)
    patchFaces = foamface.getBoundaryPatchFaces(name)
    boundaries = {}
    i = 0
    k = 0
    while i < len(patchFaces):
        boundaries[patchNames[k]] = []
        start = i+1
        end = start+patchFaces[i]
        i += 1
        for j in range(start, end):
            boundaries[patchNames[k]].append(
                foamMesh._foamFaceLabelToUuid[patchFaces[j]])
            i += 1
        k += 1
    foamMesh._boundaries = boundaries

    return foamMesh
