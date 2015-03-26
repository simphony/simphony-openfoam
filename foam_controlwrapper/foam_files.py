""" foam_files module

Module for writing and modifying OpenFOAM files

"""
import os
from foam_templates import head, dictionaryTemplates
from foam_templates import scalarTemplates, vectorTemplates
from simphony.core.cuba import CUBA
from cuba_extension import CUBAExt
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import Vector


class FoamFiles():
    """ Module for OpenFOAM file management

    """

    def __init__(self):
        pass

    def create_file_content(self, path, solver, time, writeFields):
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
        version = '2.2'
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
            for foamFile in scalarTemplates[solver]:
                if not os.path.exists(os.path.join(path, time, foamFile)):
                    foamClass = 'volScalarField'
                    location = '\"' + time + '\"'
                    foamObject = os.path.basename(foamFile)
                    heading = head.format(version=version,
                                          foamclass=foamClass,
                                          location=location,
                                          foamobject=foamObject)
                    file_name = os.path.join(time, foamFile)
                    fileContent[file_name] = heading +\
                        scalarTemplates[solver][foamFile]

            for foamFile in vectorTemplates[solver]:
                if not os.path.exists(os.path.join(path, time, foamFile)):
                    foamClass = 'volVectorField'
                    location = '\"' + time + '\"'
                    foamObject = os.path.basename(foamFile)
                    heading = head.format(version=version,
                                          foamclass=foamClass,
                                          location=location,
                                          foamobject=foamObject)
                    file_name = os.path.join(time, foamFile)
                    fileContent[file_name] = heading +\
                        vectorTemplates[solver][foamFile]

        return fileContent

    def create_directories(self, caseDirectory):
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

    def write_default_files(self, caseDirectory, solver,
                            time, writeFields):
        """Write default OpenFOAM -files base on solver
        attribute to given directory

        Parameters
        ----------
        caseDirectory : caseDirectory
            directory to write files on.
        solver : solver
            name of the solver


        Raises
        ------
        Exception when IOError occurs.

        """

        self.create_directories(caseDirectory)
        fileContent = self.create_file_content(caseDirectory, solver,
                                               time, writeFields)
        for file in fileContent:
            try:
                full_name = os.path.join(caseDirectory, file)
                with open(full_name, 'w') as f:
                    f.write(fileContent[file])
            except IOError:
                error_str = "Can't write file: {}"
                raise IOError(error_str.format(os.path.join(caseDirectory,
                                                            file)))

    def modify_files(self, case, startTime, SP, BC, solver, SPExt):
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

        """

        dire = SolutionDirectory(case, archive="SimPhoNy")

        nOfTimeSteps = SP[CUBA.NUMBER_OF_TIME_STEPS]
        deltaT = SP[CUBA.TIME_STEP]
        endTime = int(startTime) + nOfTimeSteps*deltaT
        writeInterval = endTime-int(startTime)
        # if empty type boundary condition is used this must be
        # changed in foam's polyMesh/boundary file to the same
        # type (default type is patch)
        pressureBCs = BC[CUBA.PRESSURE]
        emptyBoundaries = []
        for boundary in pressureBCs:
            if pressureBCs[boundary] == "empty":
                emptyBoundaries.append(boundary)

        boundaryFile = os.path.join(case, 'constant', 'polyMesh', 'boundary')
        try:
            control = ParsedParameterFile(boundaryFile, boundaryDict=True)
            boundaries = control.content
            for boundary in emptyBoundaries:
                for bi in range(len(boundaries)):
                    filb = boundaries[bi]
                    if filb == boundary:
                        boundaries[bi+1]['type'] = "empty"
                        break

        except IOError:
            error_str = "File {} does not exist"
            raise IOError(error_str.format(boundaryFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise IOError(error_str.format(control))

# parse system/controlDict -file in case directory
        parFile = os.path.join(case, 'system', 'controlDict')
        try:
            control = ParsedParameterFile(parFile)
            control["startTime"] = startTime
            control["endTime"] = endTime
            control["deltaT"] = deltaT
            control["writeInterval"] = writeInterval
        except IOError:
            error_str = "File {} does not exist"
            raise IOError(error_str.format(parFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise IOError(error_str.format(control))

        if solver == "interFoam":
            density = SP[CUBA.DENSITY]
            viscosity = SP[CUBA.DYNAMIC_VISCOSITY]

            # parse constant/transportProperties -file in case directory
            parFile = os.path.join(case, 'constant', 'transportProperties')
            try:
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
                control["sigma"][2] = SPExt[CUBAExt.SURFACE_TENSION]
            except IOError:
                error_str = "File {} does not exist"
                raise IOError(error_str.format(parFile))
                try:
                    control.writeFile()
                except IOError:
                    error_str = "Can't write file with content: {}"
                    raise IOError(error_str.format(control))
        else:
            density = SP[CUBA.DENSITY]
            viscosity = SP[CUBA.DYNAMIC_VISCOSITY]
            kinematicViscosity = viscosity/density

            # parse constant/transportProperties -file in case directory
            parFile = os.path.join(case, 'constant', 'transportProperties')
            try:
                control = ParsedParameterFile(parFile)
                control["nu"][2] = kinematicViscosity
            except IOError:
                error_str = "File {} does not exist"
                raise IOError(error_str.format(parFile))
            try:
                control.writeFile()
            except IOError:
                error_str = "Can't write file with content: {}"
                raise IOError(error_str.format(control))

        velocityBCs = BC[CUBA.VELOCITY]
        # parse startTime/U -file in case directory
        parFile = os.path.join(case, str(startTime), 'U')
        try:
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
                else:
                    control["boundaryField"][boundary]["type"] = \
                        "fixedValue"
                    valueString = "uniform ( " \
                                  + str(velocityBCs[boundary][0]) + " " \
                                  + str(velocityBCs[boundary][1]) + " " \
                                  + str(velocityBCs[boundary][2]) + " )"
                    control["boundaryField"][boundary]["value"] = \
                        valueString
        except IOError:
            error_str = "File {} does not exist"
            raise IOError(error_str.format(parFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise IOError(error_str.format(control))

        pressureBCs = BC[CUBA.PRESSURE]

# parse startTime/p -file in case directory
        pname = 'p'
        if solver == "interFoam":
            pname = 'p_rgh'
        parFile = os.path.join(case, str(startTime), pname)
        try:
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
                else:
                    control["boundaryField"][boundary]["type"] = \
                        "fixedValue"
                    valueString = "uniform "+str(pressureBCs[boundary])
                    control["boundaryField"][boundary]["value"] = \
                        valueString
        except IOError:
            error_str = "File {} does not exist"
            raise IOError(error_str.format(parFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise IOError(error_str.format(control))

        if solver == "interFoam":
            volumeFractionBCs = BC[CUBA.VOLUME_FRACTION]
            # parse startTime/alpha1 -file in case directory
            parFile = os.path.join(case, str(startTime), 'alpha1')
            try:
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
                    else:
                        control["boundaryField"][boundary]["type"] = \
                            "fixedValue"
                        valueString = "uniform " + \
                                      str(volumeFractionBCs[boundary])
                        control["boundaryField"][boundary]["value"] = \
                            valueString
            except IOError:
                error_str = "File {} does not exist"
                raise IOError(error_str.format(parFile))
            try:
                control.writeFile()
            except IOError:
                error_str = "Can't write file with content: {}"
                raise IOError(error_str.format(control))

        return dire

    def set_all_cell_data(self, path, time_name, data_name,
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

        try:
            dir_name = os.path.join(path, time_name)
            control = ParsedParameterFile(os.path.join(dir_name, data_name))

            field_str = 'nonuniform List<' + value_type + '>\n'
            field_str += str(len(values)) + '\n' + '(' + '\n'
            for value in values:
                field_str += str(value).replace(',', '') + '\n'
            field_str += ')'

            control['internalField'] = field_str

        except IOError:
            error_str = "Can't read file: {}"
            raise IOError(error_str.format(os.path.join(dir_name,
                                                        data_name)))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file: {}"
            raise IOError(error_str.format(os.path.join(dir_name,
                                                        data_name)))

    def set_cell_data(self, path, time_name,
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

        try:
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

        except IOError:
            error_str = "Can't read file: {}"
            raise IOError(error_str.format(os.path.join(dir_name,
                                                        data_name)))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file: {}"
            raise IOError(error_str.format(os.path.join(dir_name,
                                                        data_name)))

    def get_all_cell_data(self, path, time_name, data_name):
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

        try:
            dir_name = os.path.join(path, time_name)
            control = ParsedParameterFile(os.path.join(dir_name, data_name))

            values = control['internalField']
            return values.val

        except IOError:
            error_str = "Can't read file: {}"
            raise IOError(error_str.format(os.path.join(dir_name,
                                                        data_name)))

    def get_cell_data(self, path, time_name, data_name, label):
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

        try:
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
        except IOError:
            error_str = "Can't read file: {}"
            raise IOError(error_str.format(os.path.join(dir_name,
                                                        data_name)))

    def get_cell_data_names(self, path, time_name):
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

    def remove_parser_files(self, path):
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
