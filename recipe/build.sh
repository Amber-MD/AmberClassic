#!/bin/sh

./configure --conda --openblas

cd src
make serial
cd ..

rsync -a README.md LICENSE config_testing.h AmberClassic.sh dat bin lib test $PREFIX
mkdir -p $PREFIX/doc
rsync -a doc/AmberClassic.pdf $PREFIX/doc
