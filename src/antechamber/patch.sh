#!/bin/sh

for program in am1bcc \
               antechamber \
               atomtype \
               bondtype \
               espgen \
               parmchk2 \
               prepgen \
               residuegen \
               respgen; do
    cat bin_template | sed "s/REPLACE_ME/$program/" > $BINDIR/$program
    chmod +x $BINDIR/$program
done
