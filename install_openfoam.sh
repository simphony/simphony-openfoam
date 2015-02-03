#!/bin/sh

# Download openfoam package
wget "http://sourceforge.net/projects/foam/files/foam/2.0.1/OpenFOAM-2.0.1.gtgz/download" -O OpenFoam-2.0.1.gtgz

popd ~

# Extract sources
mkdir OpenFoam
tar -xvf OpenFoam-2.0.1.gtgz

# build OpenFoam
cd OpenFoam-2.0.1
source /etc/bashrc
./Allwmake

pushd
