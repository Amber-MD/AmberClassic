#!/bin/sh

export MKLROOT=/opt/intel/oneapi/mkl/latest

./configure --openmp --mkl

cd src
make serial
cd ..

rsync -av README.md LICENSE AmberClassic.sh bin lib $PREFIX
