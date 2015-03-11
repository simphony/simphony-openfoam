#!/bin/bash
set -e
# install pyFoam 
# fetch pyFoam
svn co https://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder/other/scripting/PyFoam/
# move to installation directory
cd PyFoam
python setup.py install
cd ..
source /opt/openfoam222/etc/bashrc
python check_PyFoam.py
# make and setup interface library
cd openfoam-interface
wmake libso
python setup.py install
