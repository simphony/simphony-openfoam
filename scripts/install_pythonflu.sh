#!/bin/bash

# Download pythonflu packages
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/confFlu_r2.1-swig_1.tar.gz/download" -O confFlu_r2.1-swig_1.tar.gz
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/managedFlu_r1.1_1.tar.gz/download" -O managedFlu_r1.1_1.tar.gz
wget "http://sourceforge.net/projects/pythonflu/files/9.1-SWIG/pythonFlu_r9.1-SWIG_1.tar.gz/download" -O pythonFlu_r9.1-SWIG_1.tar.gz

mkdir Flu
cd Flu


# Build confFlu
mkdir confFlu
cd confFlu
tar -xf ../../confFlu_r2.1-swig_1.tar.gz
./build_configure
./configure
make install

# Build managedFlu
mkdir ../managedFlu
cd ../managedFlu
tar -xf ../../managedFlu_r1.1_1.tar.gz
cp ../../scripts/env_managedFlu.sh  env.sh
bash -c "source env.sh && ./build_configure && ./configure && make"
bash -c "source ./bashrc  && make install"

# Build pythonFlu
mkdir ../pythonFlu
cd ../pythonFlu
tar -xf ../../pythonFlu_r9.1-SWIG_1.tar.gz
cp ../../scripts/env_pythonFlu.sh  env.sh
bash -c "source env.sh && ./build_configure && ./configure && make"
bash -c "source ./bashrc  && make install"
