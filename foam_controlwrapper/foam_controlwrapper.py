""" foam_controlwrapper module

Wrapper module for OpenFOAM control using pyFoam -package

"""
from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA
from PyFoam.Execution.ConvergenceRunner import ConvergenceRunner
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.LogAnalysis.BoundingLogAnalyzer import BoundingLogAnalyzer
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import os


class FoamControlWrapper(object):
    def __init__(self):

        self.CM = DataContainer()
        self.BC = DataContainer()
        self.SP = DataContainer()

    def run(self):

        case = self.CM[CUBA.NAME]
        dire = SolutionDirectory(case, archive="SimPhoNy")
        dire.clearResults()
#        solver = self.CM[CUBA.SOLVER]
        solver = "simpleFoam"
#        startTime = self.SP["startTime"]
        startTime = 0
        nOfTimeSteps = self.SP[CUBA.NUMBER_OF_TIME_STEPS]
        deltaT = self.SP[CUBA.TIME_STEP]
        endTime = nOfTimeSteps*deltaT
        writeInterval = endTime-startTime
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

        density = self.SP[CUBA.DENSITY]
        viscosity = self.SP[CUBA.DYNAMIC_VISCOSITY]
        kinematicViscosity = viscosity/density

# parse constant/transportProperties -file in case directory
        parFile = os.path.join(case, 'constant', 'transportProperties')
        try:
            control = ParsedParameterFile(parFile)
            nu = control["nu"]
            nu[2] = kinematicViscosity
        except IOError:
            error_str = "File {} does not exist"
            raise ValueError(error_str.format(parFile))
        try:
            control.writeFile()
        except IOError:
            error_str = "Can't write file with content: {}"
            raise ValueError(error_str.format(control))

        velocityBCs = self.BC[CUBA.VELOCITY]

# parse startTime/U -file in case directory
        parFile = os.path.join(case, str(startTime), 'U')
        try:
            control = ParsedParameterFile(parFile)
            for boundary in control["boundaryField"]:
                if boundary in velocityBCs:
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

        velocityBCs = self.BC[CUBA.PRESSURE]

# parse startTime/p -file in case directory
        parFile = os.path.join(case, str(startTime), 'p')
        try:
            control = ParsedParameterFile(parFile)
            for boundary in control["boundaryField"]:
                if boundary in velocityBCs:
                    if velocityBCs[boundary] == "zeroGradient":
                        control["boundaryField"][boundary]["type"] = \
                            "zeroGradient"
                        control["boundaryField"][boundary]["value"] = \
                            "uniform 0"
                    elif velocityBCs[boundary] == "empty":
                        control["boundaryField"][boundary]["type"] = \
                            "empty"
                    else:
                        control["boundaryField"][boundary]["type"] = \
                            "fixedValue"
                        valueString = "uniform "+str(velocityBCs[boundary])
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

        run = ConvergenceRunner(BoundingLogAnalyzer(),
                                argv=[solver, "-case", case],
                                logname="SimPhoNy",
                                silent=True,
                                noLog=True)
        run.start()

        os.remove('PlyParser_FoamFileParser_parsetab.py')

        return dire.getLast()
