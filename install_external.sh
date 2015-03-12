#!/bin/bash
set -e

# Fetch PyFoam
svn co https://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder/other/scripting/PyFoam/

# Build PyFoam
pushd PyFoam
python setup.py install
popd

# Check PyFoam
source /opt/openfoam222/etc/bashrc
python check_PyFoam.py

# Install openfoam-interface
pushd openfoam-interface
wmake libso
python setup.py install
popd
