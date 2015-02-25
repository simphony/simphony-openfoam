import os

from setuptools import setup, find_packages

with open('README.rst', 'r') as readme:
    README_TEXT = readme.read()

VERSION = '0.001'


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

write_version_py()


# Get all example files
data_files = []
examples_folder = os.path.join('foam_controlwrapper', 'examples')
for root, dirnames, filenames in os.walk(examples_folder):
    base = os.path.relpath(root, 'foam_controlwrapper')
    for filename in filenames:
        data_files.append(os.path.join(base, filename))


setup(
    name='foam_controlwrapper',
    version=VERSION,
    author='SimPhoNy FP7 European Project',
    description='Implementation of OpenFoam wrappers',
    long_description=README_TEXT,
    packages=find_packages(),
    package_data={'foam_controlwrapper': data_files},
    install_requires=['simphony'],
    entry_points={
        'simphony.engine': ['openfoam = foam_controlwrapper']}
)
