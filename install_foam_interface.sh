#!/bin/bash
set -e

pushd .
cd openfoam-interface/libs
./Allwmake
cd ../wrapper
wclean
wmake libso
python setup.py install
popd
