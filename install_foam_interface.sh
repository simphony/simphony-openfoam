#!/bin/bash
set -e

export FOAM_MPI_INCLUDE=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$FOAM_MPI/include

pushd .
cd openfoam-interface/simphony-solvers/driftFluxSimphonyFoam
./Allwmake
popd
pushd .
cd openfoam-interface/internal-interface/libs
./Allwmake
popd
pushd .
cd openfoam-interface/internal-interface/wrapper/relativeVelocityModels
wmake libso
popd
pushd .
cd openfoam-interface/internal-interface/wrapper
wmake libso
popd
pushd .
cd openfoam-interface/simphony-solvers/shearStressPowerLawSlipVelocity/io
wmake libso
popd
pushd .
cd openfoam-interface/simphony-solvers/shearStressPowerLawSlipVelocity/internal
wmake libso
popd
python setup.py install $1
popd
export PATH=$PATH:$PWD/openfoam-interface/internal-interface/bin
