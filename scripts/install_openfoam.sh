#!/bin/bash
set -ex

mkdir ~/OpenFOAM
cd ~/OpenFOAM

# Download and extract sources
wget "http://sourceforge.net/projects/foam/files/foam/2.0.1/OpenFOAM-2.0.1.gtgz/download" -O OpenFOAM-2.0.1.gtgz
tar -xf OpenFOAM-2.0.1.gtgz

# Build OpenFOAM
cd OpenFOAM-2.0.1
source ./etc/bashrc
export WM_MPLIB=SYSTEMOPENMPI
export WM_NCOMPPROC=2
./Allwmake | grep -v 'g++'
