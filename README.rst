simphony-openfoam
===============

The implementation of the SimPhoNy OpenFOAM -wrappers.

.. image:: https://travis-ci.org/simphony/simphony-openfoam.svg?branch=master
    :target: https://travis-ci.org/simphony/simphony-openfoam

Repository
----------

Simphony-openfoam is hosted on github: https://github.com/simphony/simphony-openfoam

Installation
------------

Package foam_controlwrapper requires python 2.7.x, OpenFOAM 2.2.2 and pyFoam 0.6.4
 
Installing the OpenFoam interface wrappers
$ sh ./install_foam_interface.sh 

Installation is based on setuptools::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop


Testing
-------

To run the full test-suite run::

    python -m unittest discover


Directory structure
-------------------

Subpackages:


- foam_controlwrapper --  wrapper class and tests for OpenFOAM -wrapping using pyFoam 
- foam_internalwrapper --  wrapper class and tests for OpenFOAM -wrapping using internal interfaces
