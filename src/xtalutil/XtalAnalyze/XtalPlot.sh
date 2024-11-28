#! /bin/bash

########################################################################
#
# Script to produce some plots of the data output by XtalAnalyze.sh. A
# new directory names 'plots' is created and all images are place there.
#
# Instructions:
# 1. Run XtalAnalyze.sh.
# 2. Run GetBfactors.py.  (optional: set BFAC_CRYST_CALPHA to "" to skip)
# 3. Copy XtalPlot.py to working directory (where you ran XtalAnalyze.sh)
#    and run.
#
########################################################################
 
#======================================================================#
#                                                                      #
#  USER: SET VARIABLES                                                 #
#======================================================================#
PLOT_PFX="lys_ortho"                     # filename prefix (for plots, etc.)
BFAC_CRYST_CALPHA=../calpha.bfactors     # experimental Bfacs
BFAC_CRYST_SDCH=../meanresidue.bfactors
PDB_FILE=asu.pdb       # same as XtalAnalyze.sh
SC_PROP="2 2 2"        # same as XtalAnalyze.sh
ROWS=8                 # number of panel rows to use on individual asu plots
COLS=4                 # number of panel cols to use on individual asu plots
                       # Note: ROWS x COLUMNS must be equal to 
                       #       (# of unit cells) x (# of symmetry elements)


#======================================================================#
# USER: NOTHING TO SET BELOW THIS LINE                                 #
#                                                                      #
#======================================================================#

XTAL_ANALYSIS_PATH=$AMBERHOME/AmberTools/src/xtalutil/Analysis

#======================================================================#
# get absolute paths and check existence of files
echo "Using files:"
# for i in BFAC_CRYST_CALPHA BFAC_CRYST_SDCH PDB_FILE
for i in PDB_FILE
do
  x=${!i}
  x=$(cd $(dirname ${x}); pwd)/$(basename ${x})
  eval "$i=\$x"
  if [ ! -e ${!i} ]; then 
    echo "Error: File ${!i} does not exist."
    exit
  else
    echo "$i = ${!i}"  
  fi
done
echo

#======================================================================#
arr=($SC_PROP)
let "UNITCELLS = ${arr[0]}* ${arr[2]}*${arr[1]}"
ASUS=`$XTAL_ANALYSIS_PATH/GetSym.py $PDB_FILE`

# volume
echo "plotting volume"
${XTAL_ANALYSIS_PATH}/plotVolume.py -Title $PLOT_PFX

# rmsd
echo "plotting rmsd"
${XTAL_ANALYSIS_PATH}/plotRmsd.py $PLOT_PFX
${XTAL_ANALYSIS_PATH}/plotRmsd_v4.py $PLOT_PFX

# per asu rmsd
echo "plotting individual asu rmsd"
${XTAL_ANALYSIS_PATH}/plotRmsdPerASU.py $PLOT_PFX $ASUS $UNITCELLS $ROWS $COLS

# bfactors
if [ ! -z "$BFAC_CRYST_CALPHA" ]; then
   echo "plotting bfactors"
   ${XTAL_ANALYSIS_PATH}/plotBfac.py \
     $PLOT_PFX \
     $BFAC_CRYST_CALPHA \
     $BFAC_CRYST_SDCH

   # individual asu's
   echo "plotting individual asu data"
   cd splittrajectories
   ${XTAL_ANALYSIS_PATH}/plotIndivASU.py \
     -Title $PLOT_PFX -u $UNITCELLS -a $ASUS \
     -b $BFAC_CRYST_CALPHA -row $ROWS -col $COLS
   cd ..
fi

# move plots to plots directory
mkdir -p plots; mv *.pdf plots
