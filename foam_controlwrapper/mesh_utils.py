""" Utility functions fo FoamMesh class

"""
from foam_files import FoamFiles
import simphonyfoaminterface as foamface
from foam_templates import head
from foam_templates import scalarTemplates, vectorTemplates
import os


def create_dummy_cellvectordata(path, name, time, data_name, dimensionset):
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
    values = [(0.0, 0.0, 0.0) for item in range(nCells)]
    solver = 'interFoam'
    foamFile = os.path.join(time, data_name)
    foamClass = 'volVectorField'
    location = '\"' + os.path.dirname(foamFile) + '\"'
    foamObject = os.path.basename(foamFile)
    heading = head.format(version=version, foamclass=foamClass,
                          location=location, foamobject=foamObject)
    fileContent = heading + vectorTemplates[solver][data_name]
    try:
        f = open(os.path.join(path, time, data_name), 'w')
        f.write(fileContent)
        f.close()
    except IOError:
        error_str = "Can't write file: {}"
        raise IOError(error_str.format(os.path.join(path,
                                                    time,
                                                    data_name)))
    FoamFiles().set_all_cell_data(path, time,
                                  data_name,
                                  values, 'vector')


def create_dummy_celldata(path, name, time, data_name, dimensionset):
    """Created dummy cell data

    Parameters
    ----------
    data_name : str
    Name of data to be created
    dimensionset : tuple
    Data dimensionset

    """

    version = '2.2'
    nCells = foamface.getCellCount(name)
    values = [0.0 for item in range(nCells)]
    solver = 'interFoam'
    foamFile = os.path.join(time, data_name)
    foamClass = 'volScalarField'
    location = '\"' + os.path.dirname(foamFile) + '\"'
    foamObject = os.path.basename(foamFile)
    heading = head.format(version=version, foamclass=foamClass,
                          location=location, foamobject=foamObject)
    fileContent = heading + scalarTemplates[solver][data_name]
    try:
        f = open(os.path.join(path, time, data_name), 'w')
        f.write(fileContent)
        f.close()
    except IOError:
        error_str = "Can't write file: {}"
        raise IOError(error_str.format(os.path.join(path,
                                                    time,
                                                    data_name)))

    FoamFiles().set_all_cell_data(path, time,
                                  data_name,
                                  values, 'scalar')
