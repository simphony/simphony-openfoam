#!/bin/sh
set -ex

mkdir ~/OpenFoam
cd ~/OpenFoam

# Download and extract sources
wget "http://sourceforge.net/projects/foam/files/foam/2.0.1/OpenFOAM-2.0.1.gtgz/download" -O OpenFoam-2.0.1.gtgz
tar -xf OpenFoam-2.0.1.gtgz

# Build OpenFoam
ls -a
cd OpenFoam-2.0.1
source ./etc/bashrc
export WM_MPLIB = SYSTEMOPENMPI
export WM_NCOMPPROCS=3
./Allwmake
