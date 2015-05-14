#!/bin/bash
set -e

pushd openfoam-interface
wmake libso
python setup.py install
popd

pushd .
cd openfoam-interface/internal-interface/libs
./Allwmake
cd ../wrapper
wclean
wmake libso
python setup.py install
popd
