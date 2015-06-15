import os

from setuptools import setup, find_packages, Command
from setuptools.command.install import install

with open('README.rst', 'r') as readme:
    README_TEXT = readme.read()

VERSION = '0.1.2.dev0'


class CleanCommand(Command):
    description = "custom clean command \
    that forcefully removes dist/build directories"
    user_options = []

    def initialize_options(self):
        self.cwd = None

    def finalize_options(self):
        self.cwd = os.getcwd()

    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in root: %s' % self.cwd
        os.system('./Allwclean')


class InstallCommand(install):
    """Customized setuptools install command - prints a friendly greeting."""
    def run(self):
        os.system("./install_foam_interface.sh")
        install.run(self)


def write_version_py(filename=None):

    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__), 'foam_controlwrapper', 'version.py')
    ver = """\
version = '%s'
"""
    fh = open(filename, 'wb')
    try:
        fh.write(ver % VERSION)
    finally:
        fh.close()

write_version_py(os.path.join(os.path.dirname(__file__),
                              'foam_controlwrapper', 'version.py'))
write_version_py(os.path.join(os.path.dirname(__file__),
                              'foam_internalwrapper', 'version.py'))

setup(
    name='foam_wrappers',
    version=VERSION,
    author='SimPhoNy FP7 European Project',
    description='Implementation of OpenFoam wrappers',
    long_description=README_TEXT,
    packages=find_packages(),
    install_requires=['simphony == 0.1.3'],
    entry_points={'simphony.engine':
                  ['openfoam_file_io = foam_controlwrapper',
                   'openfoam_internal = foam_internalwrapper']},
    cmdclass={
        'clean': CleanCommand,
        'install': InstallCommand}
)
