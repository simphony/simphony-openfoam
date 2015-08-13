#!/bin/bash
set -e
# install pyFoam
mkdir -p src/simphony-openfoam/pyfoam
wget https://openfoamwiki.net/images/3/3b/PyFoam-0.6.4.tar.gz -O src/simphony-openfoam/pyfoam/pyfoam.tgz --no-check-certificate
tar -xzf src/simphony-openfoam/pyfoam/pyfoam.tgz -C src/simphony-openfoam/pyfoam
pip install --upgrade src/simphony-openfoam/pyfoam/PyFoam-0.6.4
rm -Rf src/simphony-openfoam/pyfoam
source /opt/openfoam222/etc/bashrc
python check_PyFoam.py
