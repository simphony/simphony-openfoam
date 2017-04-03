""" Utility functions fo FoamMesh class

"""
import simphonyfoaminterface as foamface

from simphony.core.cuba import CUBA
from simphony.cuds.mesh import Cell
from simphony.cuds.meta.api import PhaseVolumeFraction

from foam_controlwrapper.foam_variables import (dataTypeMap, dataKeyMap,
                                                phaseNames)
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
    iomode = 1 if io else 0
    if dataTypeMap[dataKeyMap[data_name]] in cellDataTypes:
        dimension = dataDimensionMap[dataKeyMap[data_name]]
        if dataTypeMap[dataKeyMap[data_name]] == 'vector':
            values = [0.0] * (3 * nCells)
            foamface.setAllCellVectorData(name, data_name, iomode, values,
                                          dimension)
        elif dataTypeMap[dataKeyMap[data_name]] == 'scalar':
            values = [0.0] * nCells
            if dataKeyMap[data_name] == CUBA.VOLUME_FRACTION:
                foamface.setAllCellData(name,
                                        data_name + '.' + phaseNames[0],
                                        iomode, values, dimension)
            else:
                foamface.setAllCellData(name, data_name, iomode, values,
                                        dimension)
        elif dataTypeMap[dataKeyMap[data_name]] == 'tensor':
            values = [0.0] * (nCells * 9)
            foamface.setAllCellTensorData(name, data_name, iomode, values,
                                          dimension)


def set_cells_data(name, cells, dataNameKeyMap, materials, io=False):
    """Set data to specific cells

   Parameters
    ----------
    name : str
        name of mesh
    cells : list Cell
        list of Cells
    dataNameKeyMap : dictionary
        variables name map to CUBA keys (only variables to be saved)
    materials: dictionary
        map from phase name to material
    io : boolean
        if True write data to disk

    """

    iomode = 1 if io else 0
    for dName in dataNameKeyMap:
        dataName, _, _ = dName.partition('.')
        dataKey = dataNameKeyMap[dataName]
        dimension = dataDimensionMap[dataKey]
        if dataTypeMap[dataKey] in cellDataTypes:
            if dataTypeMap[dataKey] == "scalar":
                data = []
                if dataKey == CUBA.VOLUME_FRACTION:
                    material1 = materials[phaseNames[0]]
                    for cell in cells:
                        if cell.data[dataKey][0].material == material1:
                            data.append(cell.data[dataKey][0]
                                        .volume_fraction)
                        else:
                            data.append(cell.data[dataKey][1]
                                        .volume_fraction)
                    dName = dataName + '.' + phaseNames[0]
                    foamface.setAllCellData(name, dName, iomode, data,
                                            dimension)
                else:
                    for cell in cells:
                        data.append(cell.data[dataKey])
                    foamface.setAllCellData(name, dataName, iomode,
                                            data, dimension)
            elif dataTypeMap[dataKey] == "vector":
                data = []
                for cell in cells:
                    for val in cell.data[dataKey]:
                        data.append(val)
                foamface.setAllCellVectorData(name, dataName, iomode,
                                              data, dimension)
            elif dataTypeMap[dataKey] == "tensor":
                data = []
                for cell in cells:
                    for val in cell.data[dataKey]:
                        data.append(val)
                foamface.setAllCellTensorData(name, dataName, iomode,
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
       args[3] - datamap (key, data) pair
       args[4] - foamCellLabelToUuid
       args[5] - foamPhaseNameToMaterial
    """
    cell_start = args[0]
    cell_end = args[1]
    cells_puids = args[2]
    data_map = args[3]
    foamCellLabelToUuid = args[4]
    foamPhaseNameToMaterial = args[5]
    cells = []
    for cell_label in range(cell_start, cell_end + 1, 1):
        cell = Cell(cells_puids[cell_label],
                    foamCellLabelToUuid[cell_label])
        for dataKey, data in data_map.iteritems():
            if dataTypeMap[dataKey] == "scalar":
                if dataKey == CUBA.VOLUME_FRACTION:
                    if foamPhaseNameToMaterial:
                        material1 = foamPhaseNameToMaterial[
                            phaseNames[0]]
                        material2 = foamPhaseNameToMaterial[
                            phaseNames[1]]
                        vol_frac1 = data[cell_label]
                        phase1_vol_frac = PhaseVolumeFraction(
                            material1, vol_frac1)
                        phase2_vol_frac = PhaseVolumeFraction(
                            material2, 1 - vol_frac1)
                        cell.data[dataKey] = [phase1_vol_frac,
                                              phase2_vol_frac]
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
