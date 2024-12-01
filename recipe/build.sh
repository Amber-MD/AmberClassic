#!/bin/sh

export MKLROOT=/opt/intel/oneapi/mkl/latest

./configure --conda --openmp --openblas

cd src
make serial
cd ..

rsync -av README.md LICENSE config_testing.h AmberClassic.sh dat bin lib test $PREFIX
mkdir -p $PREFIX/doc
rsync -av doc/AmberClassic.pdf $PREFIX/doc
