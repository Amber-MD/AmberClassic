#!/bin/sh

./configure --conda --openmp

(cd src && make conda)

rsync -a README.md LICENSE AmberClassic.sh config_testing.h include dat bin lib $PREFIX
mkdir -p $PREFIX/doc
rsync -a doc/AmberClassic.pdf $PREFIX/doc
