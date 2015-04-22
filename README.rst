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

Package foam_controlwrapper requires python 2.7.x, OpenFOAM 2.2.x and pyFoam 0.6.4
 
First, install the OpenFoam internal interface wrapper
    ./install_foam_interface.sh 

Installation is based on setuptools::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop

Testing
-------

An example case can be run
    cd foam_controlwrapper/examples/
    python poiseuille.py 


Directory structure
-------------------

Subpackages:


- foam_controlwrapper --  wrapper class and tests for OpenFOAM -wrapping using pyFoam 
