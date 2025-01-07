#!/bin/sh

export CUDA_HOME=/usr/local/cuda

./configure --conda
(cd src && make conda)

./configure --conda --mpi
(cd src && make clean2 && make parallel)

if [ "`uname`" == "Linux" ]; then
   ./configure --conda --cuda
   (cd src && make clean2 && make cuda)

   ./configure --conda --cuda --mpi
   (cd src && make clean2 && make mpicuda)
fi

rsync -a README.md LICENSE AmberClassic.sh config_testing.h include dat bin lib $PREFIX
mkdir -p $PREFIX/doc
rsync -a doc/AmberClassic.pdf $PREFIX/doc

# Export AMBERCLASSICHOME automatically
mkdir -p ${PREFIX}/etc/conda/{activate,deactivate}.d
cp ${RECIPE_DIR}/activate.sh ${PREFIX}/etc/conda/activate.d/amberclassic.sh
cp ${RECIPE_DIR}/deactivate.sh ${PREFIX}/etc/conda/deactivate.d/amberclassic.sh
