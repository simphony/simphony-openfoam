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

- `simphony-common`_ ~= 0.5.0

.. _simphony-common: https://github.com/simphony/simphony-common

Installation
------------

Package foam_controlwrapper requires python 2.7.x, OpenFOAM 2.3.0 or OpenFOAM 2.4.0
 

Before installing or using simphony-openfoam , make sure the OpenFOAM environment variables are set by the following command for OpenFOAM 2.3.0
    source /opt/openfoam230/etc/bashrc 
or for OpenFOAM 2.4.0
    source /opt/openfoam240/etc/bashrc 



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


- foam_controlwrapper --  wrapper class, tests and examples for OpenFOAM using IO wrapping 
- foam_internalwrapper --  wrapper class, tests and examples for OpenFOAM wrapping using internal interfaces
- openfoam-interface -- OpenFOAM interface and modified solver codes


Cleaning
-------

To run the cleaner run::

    python setup.py clean

