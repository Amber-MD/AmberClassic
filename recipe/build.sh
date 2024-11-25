#!/bin/sh

source ./AmberClassic.sh
./configure --openmp

cd src
make serial
cd ..

rsync -av bin lib $PREFIX
