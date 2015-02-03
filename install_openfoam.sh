#!/bin/sh

mkdir ~/OpenFoam
cd ~/OpenFoam

# Download and extract sources
wget "http://sourceforge.net/projects/foam/files/foam/2.0.1/OpenFOAM-2.0.1.gtgz/download" -O OpenFoam-2.0.1.gtgz
tar -xf OpenFoam-2.0.1.gtgz

# Build OpenFoam
cd OpenFoam-2.0.1
./Allwmake
