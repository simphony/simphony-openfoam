simphony-openfoam
===============

The implementation of the SimPhoNy OpenFOAM -wrappers.

.. image:: https://travis-ci.org/simphony/simphony-openfoam.svg?branch=master
    :target: https://travis-ci.org/simphony/simphony-openfoam

Repository
----------

Simphony-openfoam is hosted on github: https://github.com/simphony/simphony-openfoam

Requirements
------------

- `simphony-common`_ == 0.1.3

.. _simphony-common: https://github.com/simphony/simphony-common

Installation
------------

Package foam_controlwrapper requires python 2.7.x, OpenFOAM 2.2.2 and pyFoam 0.6.4
 
Installation is based on setuptools::

    # build and install
    python setup.py install

Testing
-------

To run the full test-suite run::

    python -m unittest discover


Directory structure
-------------------

Subpackages:


- foam_controlwrapper --  wrapper class and tests for OpenFOAM using IO wrapping 
- foam_internalwrapper --  wrapper class and tests for OpenFOAM wrapping using internal interfaces

Cleaning
-------

To run the cleaner run::

    python setup.py clean

