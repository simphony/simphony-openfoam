""" foam_runner module

Module for running OpenFoam case

"""
import subprocess
import os
import re

import simphonyfoaminterface as foamface


class FoamRunner():
    """ Module for running OpenFOAM case

    """

    def __init__(self, solver, name, case, ncores):
        self.solver = solver
        self.name = name
        self.case = case
        self.ncores = ncores

    def run(self):
        """Run specified OpenFoam solver in specified case directory
        as external command

        """
        if self.ncores < 1:
            raise ValueError('Number of cores must be greater than zero')
        elif self.ncores == 1:
            cmd = self.solver + " -case "+self.case + " > " +\
                os.path.join(self.case, 'log')
            subprocess.check_call(cmd, shell=True)
        else:
            # write and modify decomposeParDict file
            dict =\
                """
            numberOfSubdomains {ncores};
            method          scotch;
            """
            foamface.modifyDictionary(self.name, "decomposeParDict",
                                      dict.format(ncores=self.ncores))
            foamface.writeDictionary(self.name, "decomposeParDict", False)
            # decompose
            cmd = "decomposePar -force -case " + self.case + " > " +\
                  os.path.join(self.case, 'log.decompose')

            subprocess.check_call(cmd, shell=True)
            cmd = "mpirun -np " + str(self.ncores) + " " + self.solver +\
                  " -parallel -case " + self.case + " > " +\
                  os.path.join(self.case, 'log')
            subprocess.check_call(cmd, shell=True)

            # reconstruct
            cmd = "reconstructPar  -case " + self.case + " > " +\
                  os.path.join(self.case, 'log.reconstruct')
            subprocess.check_call(cmd, shell=True)

    def get_last_time(self):
        """Return last time step in OpenFoam case directory

        """

        time_dirs = [f for f in os.listdir(self.case) if
                     re.match(r'[0-9]+.*', f)]
        time_dirs.sort(key=lambda time_dirs: float(time_dirs.split()[0]))
        return float(time_dirs.pop())
