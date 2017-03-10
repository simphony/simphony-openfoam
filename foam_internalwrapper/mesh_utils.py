""" Utility functions fo FoamMesh class

"""
import simphonyfoaminterface as foamface

from simphony.core.cuba import CUBA
from simphony.cuds.mesh import Cell

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
                values = [0.0] * (3 * nCells)
                foamface.setAllCellVectorData(name, data_name, 1, values,
                                              dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'scalar':
                values = [0.0] * nCells
                foamface.setAllCellData(name, data_name, 1, values,
                                        dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'tensor':
                values = [0.0] * (nCells * 9)
                foamface.setAllCellTensorData(name, data_name, 1, values,
                                              dimension)

    else:
        if dataTypeMap[dataKeyMap[data_name]] in cellDataTypes:
            dimension = dataDimensionMap[dataKeyMap[data_name]]
            if dataTypeMap[dataKeyMap[data_name]] == 'vector':
                values = [0.0] * (3 * nCells)
                foamface.setAllCellVectorData(name, data_name, 0, values,
                                              dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'scalar':
                values = [0.0] * nCells
                foamface.setAllCellData(name, data_name, 0, values,
                                        dimension)
            elif dataTypeMap[dataKeyMap[data_name]] == 'tensor':
                values = [0.0] * (nCells * 9)
                foamface.setAllCellTensorData(name, data_name, 0, values,
                                              dimension)


def set_cells_data(name, cells, dataNameKeyMap, io=False):
    """Set data to specific cells

   Parameters
    ----------
    name : str
        name of mesh
    cells : list Cell
        list of Cells
    io : boolean
        if True write data to disk

    """

    for dataName in dataNameKeyMap:

        dataKey = dataNameKeyMap[dataName]
        dimension = dataDimensionMap[dataKey]
        if io:
            if dataTypeMap[dataKey] in cellDataTypes:
                if dataTypeMap[dataKey] == "scalar":
                    data = []
                    for cell in cells:
                        data.append(cell.data[dataKey])
                    foamface.setAllCellData(name, dataName, 1,
                                            data, dimension)
                elif dataTypeMap[dataKey] == "vector":
                    data = []
                    for cell in cells:
                        for val in cell.data[dataKey]:
                            data.append(val)
                    foamface.setAllCellVectorData(name, dataName, 1,
                                                  data, dimension)
                elif dataTypeMap[dataKey] == "tensor":
                    data = []
                    for cell in cells:
                        for val in cell.data[dataKey]:
                            data.append(val)
                    foamface.setAllCellTensorData(name, dataName, 1,
                                                  data, dimension)
        else:
            if dataTypeMap[dataKey] in cellDataTypes:
                if dataTypeMap[dataKey] == "scalar":
                    data = []
                    for cell in cells:
                        data.append(cell.data[dataKey])
                    foamface.setAllCellData(name, dataName, 0,
                                            data, dimension)
                elif dataTypeMap[dataKey] == "vector":
                    data = []
                    for cell in cells:
                        for val in cell.data[dataKey]:
                            data.append(val)
                    foamface.setAllCellVectorData(name, dataName, 0,
                                                  data, dimension)
                elif dataTypeMap[dataKey] == "tensor":
                    data = []
                    for cell in cells:
                        for val in cell.data[dataKey]:
                            data.append(val)
                    foamface.setAllCellTensorData(name, dataName, 0,
                                                  data, dimension)


def get_cells_in_range(args):
    """ get list of cells on given label range

    Parameters
    ----------
    args: list
       list of parameters
       args[0] - cell start label
       args[1] - cell end label
       args[2] - packed list of all cells point indices
       args[3] - mesh
    """
    cell_start = args[0]
    cell_end = args[1]
    cells_puids = args[2]
    data_map = args[3]
    mesh = args[4]
    cells = []
    for cell_label in range(cell_start, cell_end + 1, 1):
        cell = Cell(cells_puids[cell_label],
                    mesh._foamCellLabelToUuid[cell_label])
        for dataKey, data in data_map.iteritems():
            if dataTypeMap[dataKey] == "scalar":
                if dataKey == CUBA.MATERIAL:
                    cell.data[dataKey] = \
                        mesh._foamMaterialLabelToUuid
                else:
                    cell.data[dataKey] = data[cell_label]
            elif dataTypeMap[dataKey] == "vector":
                cell.data[dataKey] = \
                    [data[cell_label * 3 + k] for k in range(3)]
            elif dataTypeMap[dataKey] == "tensor":
                cell.data[dataKey] = \
                    [data[cell_label * 9 + k] for k in range(9)]
        cells.append(cell)
    return cells
