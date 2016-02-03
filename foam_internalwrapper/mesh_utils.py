""" Utility functions fo FoamMesh class

"""
import simphonyfoaminterface as foamface

from foam_controlwrapper.foam_variables import (dataTypeMap, dataKeyMap)
from foam_controlwrapper.foam_variables import dataDimensionMap
from foam_controlwrapper.foam_variables import cellDataTypes


def create_dummy_celldata(name, data_name, io=False):
    """Creates dummy cell data to OpenFoams objectRegistry and
    writes to case directory if path defined

    Parameters
    ----------
    name : str
        Name of mesh
    data_name : str
        Name of data to be created
    io : boolean
        if True write data to disk

    """

    nCells = foamface.getCellCount(name)
    if io:
        if dataTypeMap[dataKeyMap[data_name]] in cellDataTypes:
            dimension = dataDimensionMap[dataKeyMap[data_name]]
            if dataTypeMap[dataKeyMap[data_name]] == 'vector':
                values = [[0.0, 0.0, 0.0] for item in range(nCells)]
                foamface.setAllCellVectorData(name, data_name, 1, values,
                                              dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'scalar':
                values = [0.0 for item in range(nCells)]
                foamface.setAllCellData(name, data_name, 1, values,
                                        dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'tensor':
                values = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                          for item in range(nCells)]
                foamface.setAllCellTensorData(name, data_name, 1, values,
                                              dimension)

    else:
        if dataTypeMap[dataKeyMap[data_name]] in cellDataTypes:
            dimension = dataDimensionMap[dataKeyMap[data_name]]
            if dataTypeMap[dataKeyMap[data_name]] == 'vector':
                values = [[0.0, 0.0, 0.0] for item in range(nCells)]
                foamface.setAllCellVectorData(name, data_name, 0, values,
                                              dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'scalar':
                values = [0.0 for item in range(nCells)]
                foamface.setAllCellData(name, data_name, 0, values,
                                        dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'tensor':
                values = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                          for item in range(nCells)]
                foamface.setAllCellTensorData(name, data_name, 0, values,
                                              dimension)


def set_cells_data(name, cells, uuidToFoamLabel, dataNameKeyMap, io=False):
    """Set data to specific cells

   Parameters
    ----------
    name : str
        name of mesh
    cells : iterator Cell
        set of Cells
    io : boolean
        if True write data to disk

    """

    for dataName in dataNameKeyMap:

        dataKey = dataNameKeyMap[dataName]
        dimension = dataDimensionMap[dataKey]
        if io:
            if dataTypeMap[dataKey] in cellDataTypes:
                data = []
                for cell in cells:
                    data.append(cell.data[dataKey])
                if dataTypeMap[dataKey] == "scalar":
                    foamface.setAllCellData(name, dataName, 1,
                                            data, dimension)
                elif dataTypeMap[dataKey] == "vector":
                    foamface.setAllCellVectorData(name, dataName, 1,
                                                  data, dimension)
                elif dataTypeMap[dataKey] == "tensor":
                    foamface.setAllCellTensorData(name, dataName, 1,
                                                  data, dimension)
        else:
            if dataTypeMap[dataKey] in cellDataTypes:
                data = []
                for cell in cells:
                    data.append(cell.data[dataKey])
                if dataTypeMap[dataKey] == "scalar":
                    foamface.setAllCellData(name, dataName, 0,
                                            data, dimension)
                elif dataTypeMap[dataKey] == "vector":
                    foamface.setAllCellVectorData(name, dataName, 0,
                                                  data, dimension)
                elif dataTypeMap[dataKey] == "tensor":
                    foamface.setAllCellTensorData(name, dataName, 0,
                                                  data, dimension)
