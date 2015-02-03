#!/bin/sh

# sudo apt-get install build-essential flex cmake zlib1g-dev qt4-dev-tools libqt4-dev gnuplot libreadline-dev libxt-dev libscotch-dev libopenmpi-dev


# Download openfoam package
wget "http://sourceforge.net/projects/foam/files/foam/2.0.1/OpenFOAM-2.0.1.gtgz/download" -O OpenFoam-2.0.1.gtgz
tar -xvf OpenFoam-2.0.1.gtgz


# Download pythonflu packages
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/confFlu_r2.1-swig_1.tar.gz/download" -O confFlu_r2.1-swig_1.tar.gz
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/managedFlu_r1.1_1.tar.gz/download" -O managedFlu_r1.1_1.tar.gz
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/pythonFlu_r9.1-SWIG_1.tar.gz/download" -O pythonFlu_r9.1-SWIG_1.tar.gz
