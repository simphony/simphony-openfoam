#!/bin/bash
set -e

pushd openfoam-interface
wmake libso
python setup.py install
popd

export FOAM_MPI_INCLUDE=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$FOAM_MPI/include

pushd .
cd openfoam-interface/internal-interface/libs
./Allwmake
cd ../wrapper
wclean
wmake libso
python setup.py install
popd

export PATH=$PATH:$PWD/openfoam-interface/internal-interface/bin

pushd .
cd openfoam-interface/internal-interface
mpicc worker.c -o bin/worker
popd

