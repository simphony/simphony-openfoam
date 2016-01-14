""" Utility functions fo FoamMesh class

"""
import simphonyfoaminterface as foamface

from foam_controlwrapper.foam_variables import (dataTypeMap, dataKeyMap)
from foam_controlwrapper.foam_variables import cellDataTypes


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
    if dataTypeMap[dataKeyMap[data_name]] in cellDataTypes:
        if dataTypeMap[dataKeyMap[data_name]] == 'vector':
            values = [[0.0, 0.0, 0.0] for itemVector in range(nCells)]
            foamface.setAllCellVectorData(name,
                                          data_name,
                                          list(values))
        elif dataTypeMap[dataKeyMap[data_name]] == 'scalar':
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
        if dataTypeMap[dataKey] in cellDataTypes:
            if dataTypeMap[dataKey] == "scalar":
                for cell in cells:
                    foamface.setCellData(name,
                                         uuidToFoamLabel[cell.uid],
                                         dataName, cell.data[dataKey])
            elif dataTypeMap[dataKey] == "vector":
                for cell in cells:
                    foamface.setCellVectorData(name,
                                               uuidToFoamLabel[cell.uid],
                                               dataName,
                                               list(cell.data[dataKey]))
            elif dataTypeMap[dataKey] == "tensor":
                for cell in cells:
                    foamface.setCellTensorData(name,
                                               uuidToFoamLabel[cell.uid],
                                               dataName,
                                               list(cell.data[dataKey]))
