#!/bin/bash
set -e

pushd openfoam-interface
wmake libso
python setup.py install
popd