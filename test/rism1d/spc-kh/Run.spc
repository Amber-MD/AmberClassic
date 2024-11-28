#!/bin/bash


cat > spc.inp <<EOF
&PARAMETERS
outlist='uxghcbtensq', THEORY='DRISM', closure='KH',
exchem_SC=1, exchem_SM=1
extra_precision=0
!grid
NR=16384, DR=0.025, rout=100, kout=30.,
!MDIIS
mdiis_nvec=20, mdiis_del=0.3, tolerance=1.e-12,
!iter
ksave=0, progress=1, maxstep=10000,
!ElStat
SMEAR=1, ADBCOR=0.5,
!bulk solvent properties
temperature=298, DIEps=78.497,
NSP=1
/
&SPECIES
!corresponds very closely to 0.0333 1/A3
DENSITY=55.296d0,
MODEL="../cSPCE.mdl"
/
EOF

../../../bin/rism1d spc > spc.out || { echo "  $0:  Program error"; exit 1; }

../check1d spc
/bin/rm -f spc.inp spc.out spc.sav
