""" Utility functions fo FoamMesh class

"""
import simphonyfoaminterface as foamface

from .foam_dicts import (dataTypeMap, dataKeyMap)


def create_dummy_celldata(name, data_name, dimensionset):
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
        #                values = [(0.0, 0.0, 0.0) for item in range(nCells)]
        #                foamface.setAllCellVectorData(name,
        #                                              data_name,
        #                                              list(dimensionset),
        #                                              list(values))
        for item in range(nCells):
            foamface.setCellVectorData(name,
                                       item,
                                       data_name, list((0.0, 0.0, 0.0)))
    else:
        # this seems to break when setAllCellVectorData call is made (bug)
        #                values = [0.0 for item in range(nCells)]
        #                foamface.setAllCellData(name,
        #                                        data_name,
        #                                        list(dimensionset),
        #                                        list(values))
        for item in range(nCells):
            foamface.setCellData(name,
                                 item,
                                 data_name, 0.0)


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
