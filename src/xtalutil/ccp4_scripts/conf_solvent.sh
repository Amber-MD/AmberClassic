#!/bin/sh

#######################################################################
# take a solute structure ($2), add a solvent map ($3);               #
#     test these with SFALL/REFMAC against the observed SF ($4)       #
#                                                                     #
# Usage: conf_solvent.sh <id> <solute-pdb> <solvent-map> <exp-mtz>    #
#######################################################################

id=$1

grid=`echo | mapdump mapin $3 | awk '/Grid sampling on x, y, z ../{print $8,$9,$10}'`


# =================== run SFALL
sfall mapin $3 hklout sfalled.mtz << EOF
mode sfcalc mapin
reso 0.94
EOF

# =================== rename column labels for refmac
cad hklin1 sfalled.mtz hklout solvent.mtz << EOF
labin file 1 E1=FC E2=PHIC
labou file 1 E1=Fpart1 E2=PHIpart1
EOF

# =================== convert back into a map for sanity check
fft hklin solvent.mtz mapout fftback.map << EOF
labin F1=Fpart1 PHI=PHIpart1
grid $grid
EOF

# =========== compute Pearson CC between original and re-calculated map
correlate.com $3 fftback.map | tee correlation.$id.log


# ========  combine observed data with solvent density as Fpart for refmac
cad hklin1 $4 hklin2 solvent.mtz hklout refme.mtz << EOF
labin file 1 all
labin file 2 all
resolution over_all 0.965
EOF

# ========== run refmac with the md-derived solvent and 
#            disable the default flat bulk solvent
refmac5 hklin refme.mtz xyzin $2 xyzout refmacout_$id.pdb \
   hklout refmacout_$id.mtz  << EOF | tee refmac_$id.log 
LABIN FP=FP SIGFP=SIGFP FPART1=Fpart1 PHIP1=PHIpart1  FREE=FreeR_flag
SOLVENT NO
SCPART 1
NCYC 40
EOF

# ========== as a control, run refmac with the defaults
refmac5 hklin refme.mtz xyzin $2 xyzout refmacout_flat.pdb \
    hklout refmacout_flat.mtz  << EOF | tee refmac_flat.log
LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag
NCYC 40
EOF

if false; then
# =========== try using both flat and expl bulk solvent at the same time, 
#             their relative weights will be their partial-strucutre scales
refmac5 hklin refme.mtz xyzin $2 xyzout refmacout_both.pdb \
   hklout refmacout_both.mtz  << EOF | tee refmac_both.log
LABIN FP=FP SIGFP=SIGFP FPART1=Fpart1 PHIP1=PHIpart1  FREE=FreeR_flag
SCPART 1
NCYC 40
EOF
fi

# --- clean up:
/bin/rm -f fftback.map refme.mtz sfalled.mtz solvent.mtz
