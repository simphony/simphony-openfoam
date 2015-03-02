#!/bin/bash

# install pyFoam as a regular user
# fetch pyFoam
cd ~
svn co https://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder/other/scripting/PyFoam/
# move to installation directory
cd PyFoam
python setup.py install


