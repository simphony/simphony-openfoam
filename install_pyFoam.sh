#!/bin/bash
# fetch pyFoam
svn co https://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder/other/scripting/PyFoam/
# move to installation directory
cd PyFoam
# install as a root
sudo python setup.py install
