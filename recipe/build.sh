#!/bin/sh

export CUDA_HOME=/usr/local/cuda

# ./configure --conda --openmp
# (cd src && make conda)

# ./configure --conda --mpi
# (cd src && make clean && make parallel)

./configure --conda --cuda
(cd src && make clean && make cuda)

rsync -a README.md LICENSE AmberClassic.sh config_testing.h include dat bin lib $PREFIX
mkdir -p $PREFIX/doc
rsync -a doc/AmberClassic.pdf $PREFIX/doc
