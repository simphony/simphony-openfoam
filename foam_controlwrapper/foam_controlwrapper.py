""" foam_controlwrapper module

Wrapper module for OpenFOAM control using pyFoam -package

"""
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from PyFoam.Execution.ConvergenceRunner import ConvergenceRunner
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.LogAnalysis.BoundingLogAnalyzer import BoundingLogAnalyzer
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import os
from foam_files import FoamFiles
from cuba_extension import CUBAExt


class FoamControlWrapper(ABCModelingEngine):
    """ Wrapper to OpenFOAM control using pyFoam -package

    """

    def __init__(self):

        self.CM = DataContainer()
        self.BC = DataContainer()
        self.SP = DataContainer()
#: to be able to use CUBAExt keywords, which are not in accepted
#  CUBA keywords this extension to CM is used
        self.CM_extensions = {}

    def run(self):
        """ run OpenFoam based on CM, BC and SP data

        Returns
        -------
        lastTime : lastTime
            Last time step taken.

        Raises
        ------
        Exception when solver not supported.
        Exception when file IO occurs.

        """

        case = self.CM[CUBA.NAME]

        GE = self.CM_extensions[CUBAExt.GE]
        solver = "simpleFoam"
        if CUBAExt.LAMINAR_MODEL in GE and not(CUBAExt.VOF in GE):
            solver = "simpleFoam"
        else:
            error_str = "GE does not define supported solver: GE = {}"
            raise NotImplementedError(error_str.format(GE))

        turbulent = 'Turbulent' if not (CUBAExt.LAMINAR_MODEL in GE) else ''
        FoamFiles().writeFoamFiles(case, (solver+turbulent))

        dire = SolutionDirectory(case, archive="SimPhoNy")
        dire.clearResults()

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

        pressureBCs = self.BC[CUBA.PRESSURE]

# parse startTime/p -file in case directory
        parFile = os.path.join(case, str(startTime), 'p')
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

        run = ConvergenceRunner(BoundingLogAnalyzer(),
                                argv=[solver, "-case", case],
                                logname="SimPhoNy",
                                silent=True,
                                noLog=True)
        run.start()

        os.remove('PlyParser_FoamFileParser_parsetab.py')

        return dire.getLast()

    def get_particle_container(self, name):
        raise NotImplementedError()

    def add_particle_container(self, particle_container):
        raise NotImplementedError()

    def delete_particle_container(self, name):
        raise NotImplementedError()

    def iter_particle_containers(self, names=None):
        raise NotImplementedError()

    def add_lattice(self, lattice):
        raise NotImplementedError()

    def add_mesh(self, mesh):
        raise NotImplementedError()

    def delete_lattice(self, name):
        raise NotImplementedError()

    def delete_mesh(self, name):
        raise NotImplementedError()

    def get_lattice(self, name):
        raise NotImplementedError()

    def get_mesh(self, name):
        raise NotImplementedError()

    def iter_lattices(self, names=None):
        raise NotImplementedError()

    def iter_meshes(self, names=None):
        raise NotImplementedError()
