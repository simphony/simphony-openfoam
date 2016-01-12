""" foam_files module

Module for writing and modifying OpenFOAM files

"""
import os
from .foam_templates import (head, dictionaryTemplates)
from .foam_templates import dataTemplates
from .foam_templates import multiphase_solvers
from .foam_variables import (dataKeyMap, dataTypeMap)
from .foam_variables import (foamTypeMap, cellDataTypes)
from simphony.core.cuba import CUBA
from cuba_extension import CUBAExt
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import Vector


def create_file_content(path, solver, time, writeFields):
    """Create content mapping to files

    Parameters
    ----------
    path : str
        path to case directory
    solver : str
        name of the solver
    time : str
        time name
    writeFields : bool
        flag to tell whether field variables are written or not

    """
    version = '2.4'
    foamClass = 'dictionary'
    fileContent = {}
    for foamFile in dictionaryTemplates[solver]:
        foamClass = 'dictionary'
        location = '\"' + os.path.dirname(foamFile) + '\"'
        foamObject = os.path.basename(foamFile)
        heading = head.format(version=version, foamclass=foamClass,
                              location=location, foamobject=foamObject)
        fileContent[foamFile] = heading +\
            dictionaryTemplates[solver][foamFile]

    if writeFields:
        for foamFile in dataTemplates[solver]:
            if not os.path.exists(os.path.join(path, time, foamFile)):
                foamClass = foamTypeMap[dataTypeMap[dataKeyMap[foamFile]]]
                location = '\"' + time + '\"'
                foamObject = os.path.basename(foamFile)
                heading = head.format(version=version,
                                      foamclass=foamClass,
                                      location=location,
                                      foamobject=foamObject)
                file_name = os.path.join(time, foamFile)
                fileContent[file_name] = heading +\
                    dataTemplates[solver][foamFile]

    return fileContent


def create_directories(caseDirectory):
    """Create default directories

    Parameters
    ----------
    caseDirectory : str
        directory to create directories on.

    """

    directories = ('constant', 'system', '0',
                   os.path.join('constant', 'polyMesh'))
    for dir in directories:
        directory = os.path.join(caseDirectory, dir)
        if not os.path.exists(directory):
            os.makedirs(directory)


def write_default_files(caseDirectory, solver, time, writeFields):
    """Write default OpenFOAM -files base on solver
    attribute to given directory

    Parameters
    ----------
    caseDirectory : caseDirectory
        directory to write files on.
    solver : solver
        name of the solver

    """

    create_directories(caseDirectory)
    fileContent = create_file_content(caseDirectory, solver, time, writeFields)
    for file in fileContent:
        full_name = os.path.join(caseDirectory, file)
        with open(full_name, 'w') as f:
            f.write(fileContent[file])


def modify_files(case, startTime, SP, BC, solver, SPExt, CMExt):
    """Modify OpenFoam case files according to user settings

    Parameters
    ----------
    case : str
        case directory
    startTime : str
        start time
    SP : DataContainer
        System Parameters
    BC : DataContainer
        Boundary Conditions
    solver : str
        solver name
    SPExt : dictionary
        extension to SP

    CMExt : dictionary
        extension to CM

    Raises
    ------
    Exception if needed parameters not defined

    """

    nOfTimeSteps = SP[CUBA.NUMBER_OF_TIME_STEPS]
    deltaT = SP[CUBA.TIME_STEP]
    endTime = float(startTime) + nOfTimeSteps*deltaT
    writeInterval = nOfTimeSteps*deltaT
    # if empty type boundary condition is used this must be
    # changed in foam's polyMesh/boundary file to the same
    # type (default type is patch)

    if solver in multiphase_solvers:
        pressureBCs = BC[CUBA.CONCENTRATION]
    else:
        pressureBCs = BC[CUBA.PRESSURE]
    emptyBoundaries = []
    for boundary in pressureBCs:
        if pressureBCs[boundary] == "empty":
            emptyBoundaries.append(boundary)

    boundaryFile = os.path.join(case, 'constant', 'polyMesh', 'boundary')
    control = ParsedParameterFile(boundaryFile, boundaryDict=True)
    boundaries = control.content
    for boundary in emptyBoundaries:
        for bi in range(len(boundaries)):
            filb = boundaries[bi]
            if filb == boundary:
                boundaries[bi+1]['type'] = "empty"
                break

    control.writeFile()

    # parse system/controlDict -file in case directory
    parFile = os.path.join(case, 'system', 'controlDict')
    control = ParsedParameterFile(parFile)
    control["startTime"] = startTime
    control["endTime"] = endTime
    control["deltaT"] = deltaT
    control["writeInterval"] = writeInterval
    control.writeFile()

    if CUBAExt.NUMBER_OF_CORES in CMExt:
        # parse system/decomposeParDict -file in case directory
        parFile = os.path.join(case, 'system', 'decomposeParDict')
        control = ParsedParameterFile(parFile)
        control["numberOfSubdomains"] =\
            CMExt[CUBAExt.NUMBER_OF_CORES]
        control.writeFile()

    if solver == "interFoam":
        density = SP[CUBA.DENSITY]
        viscosity = SP[CUBA.DYNAMIC_VISCOSITY]

        # parse constant/transportProperties -file in case directory
        parFile = os.path.join(case, 'constant', 'transportProperties')
        control = ParsedParameterFile(parFile)
        control["phase1"]["nu"][2] = \
            viscosity[SPExt[CUBAExt.PHASE_LIST][0]] / \
            density[SPExt[CUBAExt.PHASE_LIST][0]]
        control["phase1"]["rho"][2] = \
            density[SPExt[CUBAExt.PHASE_LIST][0]]
        control["phase2"]["nu"][2] = \
            viscosity[SPExt[CUBAExt.PHASE_LIST][1]] / \
            density[SPExt[CUBAExt.PHASE_LIST][1]]
        control["phase2"]["rho"][2] = \
            density[SPExt[CUBAExt.PHASE_LIST][1]]
        phases = SPExt[CUBAExt.PHASE_LIST]
        if CUBAExt.SURFACE_TENSION in SPExt:
            if phases in SPExt[CUBAExt.SURFACE_TENSION]:
                control["sigma"][2] = SPExt[CUBAExt.SURFACE_TENSION][phases]
            else:
                phases = SPExt[CUBAExt.PHASE_LIST]
                if phases in SPExt[CUBAExt.SURFACE_TENSION]:
                    control["sigma"][2] = \
                        SPExt[CUBAExt.SURFACE_TENSION][phases]
                else:
                    error_str = "Surface tension not specified"
                    raise ValueError(error_str)
        else:
            error_str = "Surface tension not specified"
            raise ValueError(error_str)

        control.writeFile()
    elif solver == "driftFluxSimphonyFoam":
        density = SP[CUBA.DENSITY]
        viscosity = SP[CUBA.DYNAMIC_VISCOSITY]

        # parse constant/transportProperties -file in case directory
        parFile = os.path.join(case, 'constant', 'transportProperties')
        control = ParsedParameterFile(parFile)
        # phase1 uses mixtureViscosity models
        viscosity_model = SPExt[CUBAExt.VISCOSITY_MODEL][
            SPExt[CUBAExt.PHASE_LIST][0]]

        control["phase1"]["transportModel"] = viscosity_model
        vmc = viscosity_model + 'Coeffs'
        viscosity_model_coeffs =\
            SPExt[CUBAExt.VISCOSITY_MODEL_COEFFS][viscosity_model]
        for coeff in viscosity_model_coeffs:
            control["phase1"][vmc][coeff] = viscosity_model_coeffs[coeff]

        control["phase1"]["rho"] = \
            density[SPExt[CUBAExt.PHASE_LIST][0]]

        viscosity_model = SPExt[CUBAExt.VISCOSITY_MODEL][
            SPExt[CUBAExt.PHASE_LIST][1]]
        if viscosity_model == 'Newtonian':
            control["phase2"]["transportModel"] = viscosity_model
            control["phase2"]["nu"] = \
                viscosity[SPExt[CUBAExt.PHASE_LIST][1]] / \
                density[SPExt[CUBAExt.PHASE_LIST][1]]
            control["phase2"]["rho"] = \
                density[SPExt[CUBAExt.PHASE_LIST][1]]
        else:
            control["phase2"]["transportModel"] = viscosity_model
            vmc = viscosity_model + 'Coeffs'
            viscosity_model_coeffs =\
                SPExt[CUBAExt.VISCOSITY_MODEL_COEFFS][viscosity_model]
            for coeff in viscosity_model_coeffs:
                control["phase2"][vmc][coeff] = viscosity_model_coeffs[coeff]

        if CUBAExt.STRESS_MODEL in SPExt:
                    control["stressModel"] =\
                        SPExt[CUBAExt.STRESS_MODEL]
        else:
            error_str = "Stress model not specified"
            raise ValueError(error_str)

        if CUBAExt.RELATIVE_VELOCITY_MODEL in SPExt:
                    control["relativeVelocityModel"] =\
                        SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL]
        else:
            error_str = "Relative velocity model not specified"
            raise ValueError(error_str)
        if SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL] == "simple"\
                or SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL] == "general":
            relVelModelCoeffs = SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL]+'Coeffs'
            if CUBAExt.RELATIVE_VELOCITY_MODEL_COEFFS in SPExt:
                coeffs = SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL_COEFFS]
                V0 = coeffs["V0"]
                control[relVelModelCoeffs]["V0"] = "( " \
                    + str(V0[0]) + " " + str(V0[1]) + " " + str(V0[2]) + " )"
                control[relVelModelCoeffs]["a"] = coeffs["a"]
                control[relVelModelCoeffs]["a1"] = coeffs["a1"]
                control[relVelModelCoeffs]["residualAlpha"] =\
                    coeffs["residualAlpha"]

        control.writeFile()
    else:
        density = SP[CUBA.DENSITY]
        viscosity = SP[CUBA.DYNAMIC_VISCOSITY]
        kinematicViscosity = viscosity/density

        # parse constant/transportProperties -file in case directory
        parFile = os.path.join(case, 'constant', 'transportProperties')
        control = ParsedParameterFile(parFile)
        control["nu"][2] = kinematicViscosity
        control.writeFile()

    if CUBAExt.EXTERNAL_BODY_FORCE_MODEL in SPExt:
        if SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL] == 'gravitation':
            # parse constant/g -file in case directory
            parFile = os.path.join(case, 'constant', 'g')
            control = ParsedParameterFile(parFile)
            g = SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL_COEFFS]["g"]
            control["value"] =\
                "( " \
                + str(g[0]) + " " + str(g[1]) + " " + str(g[2]) + " )"
            control.writeFile()

    velocityBCs = BC[CUBA.VELOCITY]
    # parse startTime/U -file in case directory
    parFile = os.path.join(case, str(startTime), 'U')
    control = ParsedParameterFile(parFile)
    for boundary in velocityBCs:
        control["boundaryField"][boundary] = {}
        if velocityBCs[boundary] == "zeroGradient":
            control["boundaryField"][boundary]["type"] = \
                "zeroGradient"
            control["boundaryField"][boundary]["value"] = \
                "uniform (0 0 0)"
        elif velocityBCs[boundary] == "empty":
            control["boundaryField"][boundary]["type"] = \
                "empty"
        elif velocityBCs[boundary] == "slip":
            control["boundaryField"][boundary]["type"] = \
                "slip"
        elif isinstance(velocityBCs[boundary], tuple):
            if velocityBCs[boundary][0] == "pressureIOVelocity":
                control["boundaryField"][boundary]["type"] = \
                    "pressureInletOutletVelocity"
                valueString = "uniform ( " \
                    + str(velocityBCs[boundary][1][0]) + " " \
                    + str(velocityBCs[boundary][1][1]) + " " \
                    + str(velocityBCs[boundary][1][2]) + " )"
                control["boundaryField"][boundary]["value"] = \
                    valueString
            elif velocityBCs[boundary][0] == "fixedValue":
                control["boundaryField"][boundary]["type"] = \
                    "fixedValue"
                valueString = "uniform ( " \
                    + str(velocityBCs[boundary][1][0]) + " " \
                    + str(velocityBCs[boundary][1][1]) + " " \
                    + str(velocityBCs[boundary][1][2]) + " )"
                control["boundaryField"][boundary]["value"] = \
                    valueString
    control.writeFile()

    # parse startTime/p -file in case directory
    pname = 'p'
    if solver in multiphase_solvers:
        pname = 'p_rgh'
    parFile = os.path.join(case, str(startTime), pname)
    control = ParsedParameterFile(parFile)
    for boundary in pressureBCs:
        control["boundaryField"][boundary] = {}
        if pressureBCs[boundary] == "zeroGradient":
            control["boundaryField"][boundary]["type"] = \
                "zeroGradient"
            control["boundaryField"][boundary]["value"] = \
                "uniform 0"
        elif pressureBCs[boundary] == "fixedFluxPressure":
            control["boundaryField"][boundary]["type"] = \
                "fixedFluxPressure"
            control["boundaryField"][boundary]["value"] = \
                "uniform 0"
        elif pressureBCs[boundary] == "empty":
            control["boundaryField"][boundary]["type"] = \
                "empty"
        elif isinstance(pressureBCs[boundary], tuple):
            if pressureBCs[boundary][0] == "fixedValue":
                control["boundaryField"][boundary]["type"] = \
                    "fixedValue"
                valueString = "uniform "+str(pressureBCs[boundary][1])
                control["boundaryField"][boundary]["value"] = \
                    valueString
    control.writeFile()

    if solver in multiphase_solvers:
        volumeFractionBCs = BC[CUBA.VOLUME_FRACTION]
        # parse startTime/alpha.phase1 -file in case directory
        parFile = os.path.join(case, str(startTime), 'alpha.phase1')
        control = ParsedParameterFile(parFile)
        for boundary in volumeFractionBCs:
            control["boundaryField"][boundary] = {}
            if volumeFractionBCs[boundary] == "zeroGradient":
                control["boundaryField"][boundary]["type"] = \
                    "zeroGradient"
                control["boundaryField"][boundary]["value"] = \
                    "uniform 0"
            elif volumeFractionBCs[boundary] == "empty":
                control["boundaryField"][boundary]["type"] = \
                    "empty"
            elif isinstance(volumeFractionBCs[boundary], tuple):
                if volumeFractionBCs[boundary][0] == "fixedValue":
                    control["boundaryField"][boundary]["type"] = \
                        "fixedValue"
                    valueString = "uniform " + \
                        str(volumeFractionBCs[boundary][1])
                    control["boundaryField"][boundary]["value"] = \
                        valueString
                elif volumeFractionBCs[boundary][0] == "inletOutlet":
                    valueString = "uniform " + \
                        str(volumeFractionBCs[boundary][1])
                    control["boundaryField"][boundary]["inletValue"] = \
                        valueString
                    control["boundaryField"][boundary]["type"] = \
                        "inletOutlet"

        control.writeFile()

    if solver == "driftFluxSimphonyFoam":
        # take boundaries from volumefractions BC's
        volumeFractionBCs = BC[CUBA.VOLUME_FRACTION]
        # parse startTime/Sigma -file in case directory
        parFile = os.path.join(case, str(startTime), 'Sigma')
        control = ParsedParameterFile(parFile)
        for boundary in volumeFractionBCs:
            control["boundaryField"][boundary] = {}
            if volumeFractionBCs[boundary] == "empty":
                control["boundaryField"][boundary]["type"] = \
                    "empty"
            else:
                control["boundaryField"][boundary]["type"] = \
                    "zeroGradient"
        control.writeFile()

        # parse startTime/Vdj -file in case directory
        parFile = os.path.join(case, str(startTime), 'Vdj')
        control = ParsedParameterFile(parFile)
        for boundary in volumeFractionBCs:
            control["boundaryField"][boundary] = {}
            if volumeFractionBCs[boundary] == "empty":
                control["boundaryField"][boundary]["type"] = \
                    "empty"
            else:
                control["boundaryField"][boundary]["type"] = \
                    "zeroGradient"
        control.writeFile()


def set_all_cell_data(path, time_name, data_name,
                      values, value_type):
    """Set cell variable values


    Parameters
    ----------
    path : str
        path to case directory
    time_name : str
        name of time directory
    data_name : str
        name of data variable
    values : list
        list of new values
    value_type : str
        OpenFOAM variable type (volScalarField, volVectorField, ...)

    """

    dir_name = os.path.join(path, time_name)
    control = ParsedParameterFile(os.path.join(dir_name, data_name))

    field_str = 'nonuniform List<' + value_type + '>\n'
    field_str += str(len(values)) + '\n' + '(' + '\n'
    for value in values:
        field_str += str(value).replace(',', '') + '\n'
    field_str += ')'

    control['internalField'] = field_str

    control.writeFile()


def set_cell_data(path, time_name,
                  data_name, label, value, value_type):
    """Set data to specific cell

   Parameters
    ----------
    path : str
        path to case directory
    time_name : str
        name of time directory
    data_name : str
        name of data variable
    label : int
        cell label
    value : list
        new value
    value_type : str
        OpenFOAM variable type (volScalarField, volVectorField, ...)


    """

    dir_name = os.path.join(path, time_name)
    control = ParsedParameterFile(os.path.join(dir_name, data_name))

    values = control['internalField']
    values[label] = value

    field_str = 'nonuniform List<' + value_type + '>\n'
    field_str += str(len(values.val)) + '\n' + '(' + '\n'
    for value in values:
        field_str += str(value).replace(',', '') + '\n'
    field_str += ')'

    control['internalField'] = field_str

    control.writeFile()


def set_cells_data(path, time_name, cells, uuidToFoamLabel, dataNameKeyMap):
    """Set data to specific cells

   Parameters
    ----------
    path : str
        path to case directory
    time_name : str
        name of time directory
    cells : iterator Cell
        set of Cells
    """

    dir_name = os.path.join(path, time_name)

    for dataName in dataNameKeyMap:
        dataKey = dataNameKeyMap[dataName]
        if dataTypeMap[dataKey] in cellDataTypes:
            control = ParsedParameterFile(os.path.join(dir_name, dataName))
            if dataTypeMap[dataKey] == "vector":
                values = [(0, 0, 0) for cell in cells]
                for cell in cells:
                    if dataKey in cell.data:
                        values[uuidToFoamLabel[cell.uid]] = \
                            tuple(cell.data[dataKey])
            elif dataTypeMap[dataKey] == "tensor":
                values = [(0, 0, 0, 0, 0, 0, 0, 0, 0) for cell in cells]
                for cell in cells:
                    if dataKey in cell.data:
                        values[uuidToFoamLabel[cell.uid]] = \
                            tuple(cell.data[dataKey])
            elif dataTypeMap[dataKey] == "scalar":
                values = [0 for cell in cells]
                for cell in cells:
                    if dataKey in cell.data:
                        values[uuidToFoamLabel[cell.uid]] = \
                            cell.data[dataKey]
            field_str = 'nonuniform List<' + dataTypeMap[dataKey] + '>\n'
            field_str += str(len(values)) + '\n' + '(' + '\n'
            for value in values:
                field_str += str(value).replace(',', '') + '\n'
            field_str += ')'
            control['internalField'] = field_str
            control.writeFile()


def get_all_cell_data(path, time_name, data_name):
    """Get specific data from every cell

    Parameters
    ----------
    path : str
        path to case directory
    time_name : str
        name of time directory
    data_name : str
        name of data variable

    """

    dir_name = os.path.join(path, time_name)
    control = ParsedParameterFile(os.path.join(dir_name, data_name))

    values = control['internalField']
    return values.val


def get_cell_data(path, time_name, data_name, label):
    """Get specific cell data for specific cell


    Parameters
    ----------
    path : str
        path to case directory
    time_name : str
        name of time directory
    data_name : str
        name of data variable
    label : int
        cell label

    """

    dir_name = os.path.join(path, time_name)
    control = ParsedParameterFile(os.path.join(dir_name, data_name))

    values = control['internalField']
    if isinstance(values.val, list):
        if isinstance(values.val[label], Vector):
            return tuple([v for v in values.val[label].vals])
        else:
            return values.val[label]
    else:
        return values.val


def get_cell_data_names(path, time_name):
    """Get cell data names

    Parameters
    ----------
    path : str
        path to case directory
    time_name : str
        name of time directory

    """

    dir_name = os.path.join(path, time_name)
    if os.path.exists(dir_name):
        file_names = []
        for file_name in os.listdir(dir_name):
            if os.path.isfile(os.path.join(dir_name, file_name)):
                file_names.append(file_name)
        return file_names
    else:
        return []


def remove_parser_files(path):
    """Remove PyFoam parser files

    Parameters
    ----------
    path : str
        path to directory

    """

    fileName = os.path.join(path, 'PlyParser_FoamFileParser_parsetab.py')
    if os.path.exists(fileName):
        os.remove(fileName)
    fileName = os.path.join(path, 'PlyParser_FoamFileParser_parsetab.pyc')
    if os.path.exists(fileName):
        os.remove(fileName)
