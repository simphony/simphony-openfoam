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

- `simphony-common`_ == 0.2.0

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



EDM egg build and uploading
---------------------------

The repository supports building of EDM packages. The command::
    
    python edmsetup.py egg

Creates the egg, and the command::

    python edmsetup.py upload_egg

uploads it to EDM repository. An important difference with respect to other EDM-enhanced repositories
is that due to the current layout of the repository, we are building **two** eggs. The consequence is that
there are **two** packageinfo.py with different package names, release tags and build numbers. 
Check under ``openfoam-interface/wrapper/`` for the secondary packageinfo.py.

The automatic Jenkins based builder creates new eggs for branches named ``release-<version>-<build>``.
To do a new build (or release), you should modify the relevant packageinfo.py, merge the PR to master, then 
create the above branch to trigger the build and upload of the new eggs.

For additional information about the EDM build system, check the ``simphonyproject/buildrecipes-common`` 
repository.
