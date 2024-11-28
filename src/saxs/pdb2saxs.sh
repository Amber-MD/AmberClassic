#!/bin/sh

#=============================================================================
#  Convert $1.pdb (via rism and rism-saxs) to $1.saxs.dat and $1.saxst.dat.
#
#  Usage:  pdb2saxs.sh <file-prefix>
#
#  N.B.: This script cannot be used to process more than one pdb file
#        in a given directory at a time!  (The guv* files would conflict,
#        and fixing this requires fixing the underlying saxs_rism code.)
#
#  Method is described here:
#
#     H. Nguyen, S.A. Pabitt, S. Meisburger, L. Pollack and D.A. Case.
#     Accurate small and wide angle X-ray scattering profiles from 
#     atomic models of proteins and nucleic acids.
#     J. Chem. Phys. 114, 22D508 (2014).
#
#  N.B.: very little error handling (yet); comment out the clean instructions
#        below if you are getting errors
#
#        This is set up for AmberTools16 and later: change the leaprc file name
#        from "leaprc.RNA.OL3" to "leaprc.ff14SB" for AmberTools15.
#
#        This example assumes a solvent of 0.1 M NaCl in water. To change
#        this, check the location and content of the "solvent" variable 
#        in step 2, below.
#
#=============================================================================

#=============================================================================
#  Step 1: use tleap to create the parm7 and rst7 files:
#=============================================================================

cat <<EOF >leap.in
set default PBRadii mbondi3
source leaprc.RNA.OL3
x = loadpdb $1.pdb
saveAmberParm x $1.parm7 $1.rst7
savePdb x $1.amber.pdb
quit
EOF

tleap -f leap.in > tleap.out
if ! [ -s $1.parm7 ]; then
   echo "tleap step failed: check tleap.out"
   exit 1
fi

#=============================================================================
#  Step 2: run 3D-rism
#=============================================================================

# Right now (June, 2016), the following file is only in the git master branch:
solvent=$AMBERHOME/dat/rism1d/cSPCE/cSPCE_KH_NaCl_0.1M.xvv

(time rism3d.snglpnt --pdb $1.amber.pdb --prmtop $1.parm7 \
               --xvv $solvent --tolerance 1e-6 --verbose 2 \
               --guv  guv --mdiis_nvec 10 --buffer 20.0 )  > $1.r3d 2>&1

if ! [ -s guv.O.1.dx ]; then
   echo "rism3d.snglpnt step failed: check $1.r3d"
   exit 1
fi

#=============================================================================
#  Step 3: run rims_saxs conversion script
#          (assumes that the only guv*.dx files in this directory are for
#           the current run; you should remove or compress these files
#           after running this, as in Step 4, below)
#=============================================================================

(time saxs_rism --grid_dir . --solute $1.amber.pdb --conc_salt 0.1 --qcut  1.0 \
     --dq 0.01 --tight --expli --exclV --output  $1.saxs.dat ) \
     > $1.saxs_rism.log 2>&1

#   file with just to q,total data:
awk 'NR>20{printf "%6.3f\t%10.5f\n", $1,log($5)/2.303}' $1.saxs.dat > $1.saxst.dat

#=============================================================================
#  Clean up:
#=============================================================================

gzip -c guv.O.1.dx > $1_kh.O.1.dx.gz
gzip -c guv.H1.1.dx > $1_kh.H1.1.dx.gz
gzip -c guv.Na+.1.dx > $1_kh.Na+.1.dx.gz
gzip -c guv.Cl-.1.dx > $1_kh.Cl-.1.dx.gz
/bin/rm -f leap.in leap.log tleap.out \
           $1.parm7 $1.rst7 $1.amber.pdb guv*.dx

