""" Utility functions for FoamMesh class

"""
import os

import simphonyfoaminterface as foamface

from .foam_files import set_all_cell_data
from .foam_templates import head
from .foam_templates import dataTemplates
from .foam_templates import (dataKeyMap, dataTypeMap, foamTypeMap)


def create_dummy_celldata(path, name, time, data_name, dimensionset):
    """Created dummy cell vector data

    Parameters
    ----------
    data_name : str
    Name of data to be created
    dimensionset : tuple
    Data dimensionset

    """

    version = '2.2'
    nCells = foamface.getCellCount(name)
    if dataTypeMap[dataKeyMap[data_name]] == 'scalar':
        values = [0.0 for item in range(nCells)]
    else:
        values = [(0.0, 0.0, 0.0) for item in range(nCells)]
    solver = 'interFoam'
    foamFile = os.path.join(time, data_name)
    foamClass = foamTypeMap[dataTypeMap[dataKeyMap[data_name]]]
    location = '\"' + os.path.dirname(foamFile) + '\"'
    foamObject = os.path.basename(foamFile)
    heading = head.format(version=version, foamclass=foamClass,
                          location=location, foamobject=foamObject)
    fileContent = heading + dataTemplates[solver][data_name]
    f = open(os.path.join(path, time, data_name), 'w')
    f.write(fileContent)
    f.close()
    set_all_cell_data(path, time, data_name, values,
                      dataTypeMap[dataKeyMap[data_name]])
