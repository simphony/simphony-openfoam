""" foam_runner module

Module for running OpenFoam case

"""
import subprocess
import os
import re


class FoamRunner():
    """ Module for running OpenFOAM case

    """

    def __init__(self, solver, case):
        self.solver = solver
        self.case = case

    def run(self):
        """Run specified OpenFoam solver in specified case directory
        as external command

        """

        cmd = self.solver + " -case "+self.case + " > " +\
            os.path.join(self.case, 'log') + '2>&1'
        subprocess.call(cmd, shell=True)

    def get_last_time(self):
        """Return last time step in OpenFoam case directory

        """

        time_dirs = [f for f in os.listdir(self.case) if
                     re.match(r'[0-9]+.*', f)]
        time_dirs.sort()
        return time_dirs.pop()
