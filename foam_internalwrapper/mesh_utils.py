""" Utility functions fo FoamMesh class

"""
import simphonyfoaminterface as foamface

from .foam_dicts import (dataTypeMap, dataKeyMap)


def create_dummy_celldata(name, data_name):
    """Created dummy cell data to OpenFoams objectRegistry

    Parameters
    ----------
    name : str
        Name of mesh
    data_name : str
        Name of data to be created
    dimensionset : list
        Data dimensionset

    """

    nCells = foamface.getCellCount(name)
    if dataTypeMap[dataKeyMap[data_name]] == 'vector':

        # this seems to break when setAllCellVectorData call is made (bug)
        values = [[0.0, 0.0, 0.0] for itemVector in range(nCells)]
        foamface.setAllCellVectorData(name,
                                      data_name,
                                      list(values))
    else:
        # this seems to break when setAllCellData call is made (bug)
        values = [0.0 for itemScalar in range(nCells)]
        foamface.setAllCellData(name,
                                data_name,
                                list(values))


def set_cells_data(name, cells, uuidToFoamLabel, dataNameKeyMap):
    """Set data to specific cells

   Parameters
    ----------
    name : str
        name of mesh
    cells : iterator Cell
        set of Cells
    """

    for dataName in dataNameKeyMap:
        dataKey = dataNameKeyMap[dataName]

        if dataTypeMap[dataKey] == "scalar":
            for cell in cells:
                foamface.setCellData(name,
                                     uuidToFoamLabel[cell.uid],
                                     dataName, cell.data[dataKey])
        else:
            for cell in cells:
                foamface.setCellVectorData(name,
                                           uuidToFoamLabel[cell.uid],
                                           dataName,
                                           list(cell.data[dataKey]))
