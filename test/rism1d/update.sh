#!/bin/bash

#change rism1d deprecated variables to current names

for f in $@; do
    echo $f
    sed -e 's/outlst/outlist/i'   -e 's/closur\([^e]\)/closure\1/i' -e 's/routup/nrout/i'    \
        -e 's/toutup/nkout/i'     -e 's/nis/mdiis_nvec/i'           -e 's/delvv/mdiis_del/i' \
        -e 's/tolvv/tolerance/i'  -e 's/kshow/progress/i'            -e 's/maxste\([^p]\)/maxstep\1/i'  \
        -e 's/temper\([^ature]\)/temperature\1/i'         -e 's/ksave=-1/ksave=0/i' \
        -e 's/psen_order/closure_order/i' \
        $f > temp.out
    mv temp.out $f
    chmod u+x $f
done