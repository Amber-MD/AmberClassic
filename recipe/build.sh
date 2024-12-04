#!/bin/sh

./configure --conda 

cd src && make serial

rsync -a README.md LICENSE config.h AmberClassic.sh include dat bin lib test $PREFIX
mkdir -p $PREFIX/doc
rsync -a doc/AmberClassic.pdf $PREFIX/doc
