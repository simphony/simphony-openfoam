#!/bin/sh

# Download pythonflu packages
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/confFlu_r2.1-swig_1.tar.gz/download" -O confFlu_r2.1-swig_1.tar.gz
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/managedFlu_r1.1_1.tar.gz/download" -O managedFlu_r1.1_1.tar.gz
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/pythonFlu_r9.1-SWIG_1.tar.gz/download" -O pythonFlu_r9.1-SWIG_1.tar.gz

mkdir Flu
cd Flu
cp ../scripts/env.sh .


# Build confFlu
mkdir confFlu
cd confFlu
tar -xf ../../confFlu_r2.1-swig_1.tar.gz
./build_configure
./configure
make install

# Build managedFlu
mkdir managedFlu
cd managedFlu
tar -xf ../../managedFlu_r1.1_1.tar.gz
cp env.sh.example env.sh
source ./env.sh
./build_configure
./configure
make
bash -c "source ./bashrc  && make install"

# Build pythonFlu
mkdir pythonFlu
cd pythonFlu
tar -xf ../../managedFlu_r1.1_1.tar.gz
source ../Flu env.sh
./build_configure
./configure
make
bash -c "source ./bashrc  && make install"
