#!/bin/bash

###############################################################################
#
# Calculate average electron density map and structure factors from an 
# MD trajectory run under constant unit cell volume conditions.
#
# Returns:
#   md_avg.map - average ASU electron density map in ccp4 map format
#   md_avg.mtz - structure factors of md_avg.map
#
# Authors v1:James Holton
#         v2: Pawel Janowski
#         v3: Dave Case (for constant volume simulations)
#
###############################################################################
#
#   Usage: make a copy of this file in your working directory (generally
#          where you store the trajectory files), then follow the instructions
#          below, editing the copied file.
#
###############################################################################
#
#   Step 1:  Use cpptraj to create pdb snapshots, e.g.:
#
#   cpptraj <<EOF
#   parm prmtop
#   reference <exp. struct in original crystal frame>
#   trajin  md1-56.nc  .......
#   rms reference <mask> norotate out fit.dat
#   trajout <$sim_pdb_dir/name.pdb> pdb multi
#   EOF
#
#   The only files in the $sim_pdb_dir (set below) should be these pdb files
#
###############################################################################
#
#   Step 2: set the following variables:
#

# directory with simulation pdb files:
sim_pdb_dir=PDBData/

# supercell multiplicity (x y z):
MD_mult=( 1 1 1 )

# ccp4 space group symbol:
SG="P1"

# grid spacing for SFALL map
GRID="GRID 120 120 120"

# flat B-factor added to frames in SFALL map calculation:
B=15

# structure factor calculation resolution limit:
reso=1.0

#  CRYST1 record from experimental PDB:
CRYST="CRYST1   27.240   31.870   34.230  88.52 108.53 111.89 P 1           1 "

###############################################################################
#
#   Step 3: run this script:  "./mdv2map.sh"
#
#  You should not need to edit anything below this line!
#
###############################################################################

function MakeMap {
####################################################################
#                      CALCULATE ED MAP                            #
####################################################################

# INPUT VARIABLES:
local pdb2map=$1        # input pdb file name
local mapname=$2        # output map name

# REFORMAT PDB FILE (unit cell cryst and space group, 
#   set b-factors, set occupancy=1, reformat atom names:

cat <<EOF > awk.in
BEGIN{B=$B; printf("$CRYST\n");}
  /^ATOM/||/^HETATM/{
     RES= substr(\$0, 18, 9);
     X = substr(\$0, 31, 8);
     Y = substr(\$0, 39, 8);
     Z = substr(\$0, 47, 8);
     Ee = \$NF;
     printf("ATOM %6d %2s   %9s    %8.3f%8.3f%8.3f  1.00%6.2f%12s\n", \
            ++n,Ee,RES,X,Y,Z,B,Ee);}
  END{print "END"}
EOF
awk -f awk.in $pdb2map > sfallme.pdb

# GET UNIT CELL PARAMETERS
pcell=( `head -1 sfallme.pdb | awk '{print $2,$3,$4,$5,$6,$7}'` )

# ADD PDB SCALE MATRIX INFORMATION
pdbset xyzin sfallme.pdb xyzout new.pdb << EOF > /dev/null
SPACE $SG
CELL ${pcell[*]}
EOF
mv new.pdb sfallme.pdb

# CALCULATE MAP
sfall xyzin sfallme.pdb mapout ${mapname} << EOF > /dev/null
mode atmmap
CELL ${pcell[*]}
SYMM $SG
FORMFAC NGAUSS 5
$GRID
EOF

#CLEAN
rm -f new.xyz awk.in
rm -f sfallme.xyz 
rm -f sfallme.pdb 
}

########################################################################
#                             MAIN                                     #
########################################################################

# SET-UP:
frames=`ls -1 $sim_pdb_dir | wc -l`  #no. of frames
symops=`awk -v SG=$SG '$4==SG{print $2;exit}' ${CLIBD}/symop.lib` #no. of sym operations
rm -rf maps
mkdir -p maps

# GET MAP FOR EACH TRAJECTORY FRAME:
for ((frame=1; frame<=$frames; frame++)); do
   sim_pdb=`ls -1 $sim_pdb_dir | sed $frame'q;d'`
   echo "Calculating map: ${sim_pdb_dir}${sim_pdb}"
   MakeMap ${sim_pdb_dir}${sim_pdb} maps/${sim_pdb}.map
done 

####################################################################
#     CALCULATE AVERAGE DENSITY MAP AND ITS STRUCTURE FACTORS      #
####################################################################

# SUM ALL MAPS:
echo "Summing maps"
rm -f sum.map
for ((frame=1; frame<=$frames; frame++)); do
   frame_map=`ls -1 maps | sed $frame'q;d'`
   if [ ! -e sum.map ]; then
      cp maps/${frame_map} sum.map
   else
      mapmask mapin1 sum.map mapin2 maps/${frame_map} mapout new.map <<EOF | grep -e "Mean density" |head -1
maps add
EOF
      mv new.map sum.map
   fi
done

# NORMALIZE SUM MAP:
scale=`echo $frames ${MD_mult[*]} $symops | awk '{print $1*$2*$3*$4*$5}'`
mapmask mapin sum.map mapout md_avg.map <<EOF | grep -e "Mean density"
scale factor `echo "1/$scale" | bc -l` 0
EOF
echo -e "=======================\n"
echo "Scaling map by 1/$scale"
echo | mapdump mapin md_avg.map | egrep dens

# CALCULATE SF FROM AVERAGE DENISTY MAP:
sfall mapin md_avg.map hklout md_avg_1.mtz << EOF > /dev/null
MODE SFCALC MAPIN
RESOLUTION $reso
EOF

# EDIT MTZ FILE TO INCLUDE SIGFP AND FOM COLUMNS:
sftools << EOF > /dev/null
read md_avg_1.mtz
set labels
FP
PHI
calc Q COL SIGFP = 0.1
calc W COL FOM = 0.9
write md_avg_2.mtz col FP SIGFP PHI FOM
y
stop
EOF

# REMOVE THE BFACTOR THAT WAS ADDED FOR SF CALCULATION TO AVOID SINGULARITY:
cad hklin1 md_avg_2.mtz hklout md_avg.mtz << EOF >/dev/null
scale file 1 1 -$B
labin file 1 all
EOF

# CLEAN:
rm -f md_avg_1.mtz 
rm -f md_avg_2.mtz
rm -f sfall.mtz
rm -f sum.map
