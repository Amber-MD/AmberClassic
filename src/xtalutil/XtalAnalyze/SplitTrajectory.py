#! /usr/bin/env python
import os
import sys
import glob
import argparse

#========================================================================#
#                                                                        #
# Splits supercell trajectory into multiple trajectories, one for each   #
# asymmetric unit. Split trajectories are placed in "splittrajectories"  #
# dir which is created in the current workig dir.                        #
#                                                                        #
#========================================================================#

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--Topology", help="Amber topology file")
parser.add_argument("-t", "--Trajectory", help="topology file")
parser.add_argument("-u", "--UnitCells", help="number of unit cells in supercell")
parser.add_argument("-a", "--ASUs", help="number of asymmetric units per unit cell")
parser.add_argument("-r", "--Residues", help="number of residues per assymetric unit")
args = parser.parse_args()

# Setup variables
UC = int(args.UnitCells)
ASU = int(args.ASUs)
RES = int(args.Residues)

# Set-up cpptraj scripts
for i in range(UC):
    for j in range(ASU):
        f = open('ctraj.split_%02d_%02d.in' % (i + 1, j + 1), 'w')
        f.write('parm %s\n' % args.Topology)
        f.write('trajin %s\n' % args.Trajectory)
        start1 = 1
        end1 = RES * j + i * ASU * RES
        start2 = (RES * j + RES + 1) + i * ASU * RES
        end2 = 99999
        if end1 < start1:
            f.write('strip \':%d-%d\' \n' % (start2, end2))
        elif start2 > RES * UC * ASU:
            f.write('strip \':%d-%d\'\n' % (start1, end1))
        else:
            f.write('strip \':%d-%d | :%d-%d\'\n' % (start1, end1, start2, end2))
        f.write('trajout %02d_%02d.nc netcdf\n' % (i + 1, j + 1))
        f.close()

# Execute ptraj to create multiple trajectories
for i in range(UC):
    for j in range(ASU):
        print("UC: %2d  ASU: %2d" % (i + 1, j + 1))
        if os.system('cpptraj <ctraj.split_%02d_%02d.in >>../XtalAnalyze.log' %
                     (i + 1, j + 1)) != 0:
            sys.exit()
