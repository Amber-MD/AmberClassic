#! /bin/bash

# Super-brief instructions: Copy this file to where your crystal
#      trajectory is.  In the copy, edit the variables below (per PAJ's 
#      instructions), then run this script. Outputs will contain 
#      a basic analysis of structures and B-factors.
#
# More complete instructions (PAJ, 1/14):
#
#  1. Prepare supercell trajectory for analysis. Such preparation could 
#     involve, for example, stripping unnecessary atoms or filtering for
#     a sparser selection of frames to speed up analysis. 
#     Also make sure to create a corresponding topology file. This is 
#     easily done in cpptraj when stripping by adding the "outprefix" 
#     flag.
#
#  2. Prepare supercell reference coordinate (rst7) file. This should have 
#     the same set of atoms as the trajectory in step 1 above and therefore 
#     should work with the topology file created in step 1. The coordinates 
#     of the first asymmetric unit in this file must be the same as the 
#     experimentallly obtained (pdb) coordinates. Easiest way to do this:
#     when first using tleap to set up the supercell for MD, use the 
#     'set default nocenter on' option. Then take the supercell rst7 file
#     before any minimization/dynamics and use cpptraj to strip desired
#     atoms to match the trajectory file from step 1.
#
#  3. Prepare a PDB file of the asymmetric unit that contains experimental pdb 
#      SMTRY and CRYST1 cards.  This should identified as ${PDB_FILE} below,
#      and the corresponding prmtop file is identified as ${PARM_FILE}.
#
#  4. Copy this script to some working directory (where output will be
#     written, edit the variables below, then run this script.Afterwards 
#     XtalPlot.py can be run to produce plots, if you wish.

#======================================================================#
#======================================================================#
#                                                                      #
#  USER: SET VARIABLES                                                 #
#======================================================================#

#  (paths may be absolute or relative to working directory;
#  files here are for a lysozyme example: change to match your filenames)

SC_PRMTOP=../alt1_uc.parm7         # supercell topology (possibly stripped)
SC_TRAJECTORY=../alt1_009.min1.x   # supercell trajectory
SC_REFERENCE=../alt1_uc.rst7       # supercell with pdb (experimental) coords 
                                   #   and correct supercell box size
PDB_FILE=asu.pdb                   # must contain SMTRY and CRYST1 cards
PARM_FILE=asu.parm7                #  (also need a matching "asu.parm7" file)
SC_PROP="2 4 2"                    # (x,y,z) used to create supercell 
ASU_NRESIDUES=11                   # number of residues in an asymmetric unit
Timestep=1.0   	                   # Trajectory frame timestep

# Amber masks for b-factor and rmsd calculations. These can be modified 
# according to the user's needs.

BM1=":1-6@CA"                 # B-factor "Calpha" mask
BM2=":1-6&!@H="               # B-factor "meanresidue" mask
RM1=":1-6&!(@H=)"             # RMSD heavy atom mask, for asu 
RM2=":1-6@CA,C,N,O"           # RMSD backbone atom mask, for asu
SCMASK="@CA,C,N"              # supercell mask for drift calculation

# location of the analysis scripts:
XTAL_ANALYSIS_PATH=$MSANDERHOME/src/xtalutil/XtalAnalyze

# List of scripts used by XtalAnalyze.sh
#       RevSym.py
#       RevSym_com.py
#       SplitTrajectory.py
#       AnalyzeRevSym.py
#       AnalyzeIndivASU.py
#       GetVolume.py
#       GetSym.py

#======================================================================#
# USER: NOTHING TO SET BELOW THIS LINE                                 #
#======================================================================#

#======================================================================#
# SETUP                                                                # 
#======================================================================#
echo
echo '#################################################################'
echo '#                          setup                                #'
echo '#################################################################'

# working directory variable, for use below:
WD=`pwd`

# get absolute paths and check existence of files
#~ python -c "import os; print os.path.abspath('$SC_PRMTOP')"
echo "Using files:"
for i in SC_PRMTOP SC_TRAJECTORY SC_REFERENCE PDB_FILE PARM_FILE
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

#calculate number of unit cells in supercell
arr=($SC_PROP)
let "UNITCELLS = ${arr[0]}* ${arr[2]}*${arr[1]}"

#calculate number of asymmetric units per unit cell
ASUS=`$XTAL_ANALYSIS_PATH/GetSym.py $PDB_FILE`

#check existence of variables
for i in "$UNITCELLS" "$ASUS"; do
  if [ -z "$i" ]; then
    echo "ERROR: variable $i is null"
    exit
  fi
done

echo "UNITCELLS=$UNITCELLS  ASUS=$ASUS"

#======================================================================#
# FIT: translate each trajectory frame to align its backbone           #
#      with SC_REFERENCE backbone.                                     #
#======================================================================#
echo
echo '#################################################################'
echo '#              Fitting supercell backbone                       #'
echo '#################################################################'
cat <<EOF > ctraj.translate.in
parm ${SC_PRMTOP}
trajin ${SC_TRAJECTORY}
trajout fit.nc netcdf
reference ${SC_REFERENCE}
rmsd reference raw "$SCMASK" nofit out drift.dat
rmsd reference no-drift "$SCMASK"  norotate out drift.dat
go
EOF
echo "TRANSLATING TRAJECTORY" >> XtalAnalyze.log
cpptraj < ctraj.translate.in >> XtalAnalyze.log
if [ $? -ne 0 ]; then
	exit
fi
echo "fit.nc netcdf supercell trajectory created"


#======================================================================#
# SPLIT: split supercell trajectories into multiple trajectories of    #
#      individual asymmetric units.                                    #
#======================================================================#
echo
echo '#################################################################'
echo '#          Splitting into asymmetric unit trajectories          #'
echo '#################################################################'
echo "SPLITTING TRAJECTORY" >> XtalAnalyze.log
rm -rf splittrajectories; mkdir splittrajectories; cd splittrajectories
${XTAL_ANALYSIS_PATH}/SplitTrajectory.py \
	-p ${SC_PRMTOP} -t ../fit.nc \
	-u ${UNITCELLS} -a ${ASUS} -r ${ASU_NRESIDUES}
cd ..

#======================================================================#
# REVSYM: reverse crystal symmetr/translation to bring each asymmetric #
#      unit into the space of the original asu.                        #
#======================================================================#
echo
echo '#################################################################'
echo '#          Reversing symmetry on split trajectories             #'
echo '#################################################################'
echo "MAKING REVSYM FILES" >> XtalAnalyze.log
rm -rf revsym; mkdir revsym; cd revsym
${XTAL_ANALYSIS_PATH}/RevSym.py \
	-p ${PDB_FILE} -r ${SC_REFERENCE} \
	-prop "$SC_PROP"
cd ..

#======================================================================#
# REVSYM_COM: translate each revsym'ed asu trajectory so that its      #
#      of mass aligns with ideal crystal lattice center of mass.       #
#======================================================================#
echo
echo '#################################################################'
echo '#   Aligning rev sym traj center of mass with crystal lattice   #'
echo '#################################################################'
echo "MAKING REVSYM FILES" >> XtalAnalyze.log
cd revsym
${XTAL_ANALYSIS_PATH}/RevSym_com.py \
	-p ${PARM_FILE} -pdb ${PDB_FILE} \
	-pdbs ${PDB_FILE} -prop "$SC_PROP"
cd ..

#======================================================================#
# ANALYZE REVSYM: analyze the revsym trajectories for rmsd, bfactors   #
#      average coordinates.                                            #
#======================================================================#
cd revsym
echo
echo '##################################################################'
echo '#             Analyzing revsym trajectories                      #'
echo '##################################################################'
echo "ANALYZING REVSYM FILES" >> XtalAnalyze.log
${XTAL_ANALYSIS_PATH}/AnalyzeRevSym.py  \
	-p ${PARM_FILE} -pdb ${PDB_FILE}  \
	-u ${UNITCELLS} -a ${ASUS}  \
	-Bm1 ${BM1}  \
	-Bm2 ${BM2}  \
	-Rm1 ${RM1}  \
	-Rm2 ${RM2}  \
	-tm ${Timestep}
cd ..

#======================================================================#
# VOLUME: get volume of each frame and compare to experimental.        #
#                                                                      #
#======================================================================#
echo 
echo '##################################################################'
echo '#                   Calculating volume                           #'
echo '##################################################################'
echo "CALCULATING VOLUME" >> XtalAnalyze.log
${XTAL_ANALYSIS_PATH}/GetVolume.py \
	-t fit.nc \
	-r ${SC_REFERENCE}


#======================================================================#
# INDIV_ASU: calculate average structure, rmsd of average structures,  #
#     and b-factors for each individual asu trajectory                 #
#======================================================================#
echo 
echo '##################################################################'
echo '#                  Analyzing individual asu trajectories         #'
echo '##################################################################'
echo "ANALYZING INDIVIDUAL ASU TRAJECTORIES" >> XtalAnalyze.log
cd splittrajectories
${XTAL_ANALYSIS_PATH}/AnalyzeIndivASU.py \
	-p ${PARM_FILE} -pdb ${PDB_FILE} \
	-u ${UNITCELLS} -a ${ASUS}  \
	-Bm1 ${BM1}  \
	-Rm1 ${RM1}
cd ..


#======================================================================#
# CLEAN                                                                #
#======================================================================#
rm -f ctraj.translate.in
rm -f splittrajectories/ctraj.*.in
rm -f revsym/ctraj.*.in

