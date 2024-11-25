#!/bin/sh

# source ./AmberClassic.sh
export MKLROOT=/opt/intel/oneapi/mkl/latest

./configure --openmp --mkl

cd src
make serial
cd ..

rsync -av bin lib $PREFIX
