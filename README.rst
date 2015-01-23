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

Package foam_wrapper requires python 2.7.x, OpenFOAM 2.1.0 and pythonFLU.
Package foam_tokenwrapper requires python 2.7.x, OpenFOAM 2.2.x
 
Installation is based on setuptools::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop


or::

    # build for in-place development
    python setupTokewrapper.py develop

Testing
-------

To run the full test-suite run::

    python -m unittest discover


Directory structure
-------------------

Subpackages:

- foam_wrapper --  wrapper class and tests for OpenFOAM -wrapper based on pythonFLU 
- foam_tokenwrapper --  wrapper class and tests for native OpenFOAM -wrapping of data transfer(mesh and solution variables) 
