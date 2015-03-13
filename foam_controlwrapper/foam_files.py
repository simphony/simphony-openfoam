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


class FoamFiles():

    def __init__(self):
        pass

    def create_file_content(self, solver):
        """ create content mapping to files

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

        for foamFile in scalarTemplates[solver]:
            foamClass = 'volScalarField'
            location = '\"' + os.path.dirname(foamFile) + '\"'
            foamObject = os.path.basename(foamFile)
            heading = head.format(version=version, foamclass=foamClass,
                                  location=location, foamobject=foamObject)
            fileContent[foamFile] = heading +\
                scalarTemplates[solver][foamFile]

        for foamFile in vectorTemplates[solver]:
            foamClass = 'volVectorField'
            location = '\"' + os.path.dirname(foamFile) + '\"'
            foamObject = os.path.basename(foamFile)
            heading = head.format(version=version, foamclass=foamClass,
                                  location=location, foamobject=foamObject)
            fileContent[foamFile] = heading +\
                vectorTemplates[solver][foamFile]

        return fileContent

    def create_directories(self, caseDirectory):
        """ create default directories

        """
        directories = ('constant', 'system', '0',
                       os.path.join('constant', 'polyMesh'))
        for dir in directories:
            directory = os.path.join(caseDirectory, dir)
            if not os.path.exists(directory):
                os.makedirs(directory)

    def write_default_files(self, caseDirectory, solver):
        """ write default OpenFOAm -files base on solver attribute to given directory

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
        fileContent = self.create_file_content(solver)
        for file in fileContent:
            try:
                f = open(os.path.join(caseDirectory, file), 'w')
                f.write(fileContent[file])
                f.close()
            except IOError:
                error_str = "Can't write file: {}"
                raise ValueError(error_str.format(os.path.join(caseDirectory,
                                                               file)))

    def modify_files(self, case, SP, BC, solver,SPExt):

        dire = SolutionDirectory(case, archive="SimPhoNy")
        dire.clearResults()

#        startTime = SP["startTime"]
        startTime = 0
        nOfTimeSteps = SP[CUBA.NUMBER_OF_TIME_STEPS]
        deltaT = SP[CUBA.TIME_STEP]
        endTime = nOfTimeSteps*deltaT
        writeInterval = endTime-startTime
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
            raise ValueError(error_str.format(boundaryFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise ValueError(error_str.format(control))

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
            raise ValueError(error_str.format(parFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise ValueError(error_str.format(control))

        if solver == "interFoam":
            density = SP[CUBA.DENSITY]
            viscosity = SP[CUBA.DYNAMIC_VISCOSITY]

            # parse constant/transportProperties -file in case directory
            parFile = os.path.join(case, 'constant', 'transportProperties')
            try:
                control = ParsedParameterFile(parFile)
#                control["twoPhase"]["phase1"] = SPExt[CUBAExt.PHASE_LIST][0]
#                control["twoPhase"]["phase2"] = SPExt[CUBAExt.PHASE_LIST][1]
                
                control["phase1"]["nu"][2] = viscosity[SPExt[CUBAExt.PHASE_LIST][0]]/density[SPExt[CUBAExt.PHASE_LIST][0]]
                control["phase2"]["nu"][2] = viscosity[SPExt[CUBAExt.PHASE_LIST][1]]/density[SPExt[CUBAExt.PHASE_LIST][1]]
            except IOError:
                error_str = "File {} does not exist"
                raise ValueError(error_str.format(parFile))
                try:
                    control.writeFile()
                except IOError:
                    error_str = "Can't write file with content: {}"
                    raise ValueError(error_str.format(control))
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
                raise ValueError(error_str.format(parFile))
            try:
                control.writeFile()
            except IOError:
                error_str = "Can't write file with content: {}"
                raise ValueError(error_str.format(control))
        
            
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
            raise ValueError(error_str.format(parFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise ValueError(error_str.format(control))


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
            raise ValueError(error_str.format(parFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise ValueError(error_str.format(control))

 
        if solver == "interFoam":
            volumeFractionBCs = BC[CUBA.VOLUME_FRACTION]
            # parse startTime/alpha1 -file in case directory
            parFile = os.path.join(case, str(startTime), 'alpha1')
            try:
                control = ParsedParameterFile(parFile)
                for boundary in pressureBCs:
                    control["boundaryField"][boundary] = {}
                    if pressureBCs[boundary] == "zeroGradient":
                        control["boundaryField"][boundary]["type"] = \
                                                                     "zeroGradient"
                        control["boundaryField"][boundary]["value"] = \
                                                                      "uniform 0"
                    elif pressureBCs[boundary] == "empty":
                        control["boundaryField"][boundary]["type"] = \
                                                                     "empty"
                    else:
                        control["boundaryField"][boundary]["type"] = \
                                                                     "fixedValue"
                        valueString = "uniform "+str(volumeFractionBCs[boundary])
                        control["boundaryField"][boundary]["value"] = \
                                                                      valueString
            except IOError:
                error_str = "File {} does not exist"
                raise ValueError(error_str.format(parFile))
            try:
                control.writeFile()
            except IOError:
                error_str = "Can't write file with content: {}"
                raise ValueError(error_str.format(control))


        return dire
