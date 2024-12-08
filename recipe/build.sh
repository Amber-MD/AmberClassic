#!/bin/sh

export CUDA_HOME=/usr/local/cuda

./configure --conda --openmp
(cd src && make conda)

./configure --conda --mpi --openmp
(cd src && make clean && make parallel)

if [ "`uname`" == "Linux" ]; then
   ./configure --conda --cuda
   (cd src && make clean && make cuda)

   ./configure --conda --cuda --mpi
   (cd src && make clean && make mpicuda)
fi

rsync -a README.md LICENSE AmberClassic.sh config_testing.h include dat bin lib $PREFIX
mkdir -p $PREFIX/doc
rsync -a doc/AmberClassic.pdf $PREFIX/doc
