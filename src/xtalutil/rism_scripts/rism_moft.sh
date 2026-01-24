#!/bin/sh

# take a single conformer from fd_helix, run RISM to get a ccp4 map file of
#     the solvent, then feed that into moft to get a water distribution

if false; then
# ================== set up in Amber
cat <<EOF > leap.in
set default nocenter on
source leaprc.RNA.OL3
x = loadpdb $1.pdb
saveamberparm x $1.parm7 $1.rst7
quit
EOF

tleap -f leap.in
ambpdb -p $1.parm7 -c $1.rst7 > $1.amber.pdb
fi

if false; then
# ============= Run rism3d.snglpnt
solvent=KCL_0.140_pse3.xvv

rism3d.snglpnt --prmtop $1.parm7 --rst $1.rst7 \
    --xvv $solvent \
    --closure kh,pse2,pse3 --tolerance 0.1,0.001,0.0000001 --solvcut 9. \
    --grdspc 0.5,0.5,0.5 --buffer 40.0 --verbose 2 --volfmt dx \
    --mdiis_del 0.3  --mdiis_nvec 10 --mdiis_restart 10. \
    --maxstep 2000 --guv KCL_pse3  > $1.KCL_pse3.r3d
fi

if false; then
# ============= Do the Laplacian convolution
  metatwist --dx KCL_pse3.O.1.dx \
      --species O  --convolve 4 --sigma 1.0 --odx KCL_pse3.O.1.dx \
     > pse3.lp
fi

if false; then
# ============= Estimate water oxygen positions
  metatwist --dx KCL_pse3.O.1.dx \
      --ldx convolution-KCL_pse3.O.1.dx --map blobsper \
      --species O WAT --bulk 55.55 --threshold 0.3  \
     > pse3.blobs
  grep -v TER KCL_pse3.O.1-convolution-KCL_pse3.O.1-blobs-centroid.pdb > wats_0.3.pdb
fi

if true; then
# ============= Add hydrogens using gwh

  sed 's/END/TER/' $1.amber.pdb > $1.wat_0.5.pdb
  gwh -p $1.parm7 -w wats_0.5.pdb < $1.amber.pdb >> $1.wat_0.5.pdb
  
fi

