#! /usr/bin/env python
import sys
import os
from numpy import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--AsuTopology", help="Amber single ASU topology file")
parser.add_argument("-pdb", "--AsuPDB", help="Single ASU PDB file")
parser.add_argument("-u", "--UnitCells", help="number of unit cells in crystal")
parser.add_argument("-a", "--ASUs", help="number of asymmetric units per unit cell")
parser.add_argument("-Bm1", help="Mask for Calpha B-factor analysis")
parser.add_argument("-Rm1", help="Suffix for heavy atom rmsd")
args = parser.parse_args()


###########################
# SETUP                   #
###########################
asymunits = int(args.ASUs)
unitcells = int(args.UnitCells)
topo = args.AsuTopology
pdb = args.AsuPDB


########################################################################
#                                                                      #
#   AVERAGE STRUCTURE calculate for each asym unit                     #
#                                                                      #
########################################################################
print("Calculating average coordinates of each asu trajectory.")
for i in range(unitcells):
    for j in range(asymunits):
        f = open('ctraj.AverageAsu.in', 'w')
        f.write('parm %s\n' % topo)
        f.write('reference %s\n' % pdb)
        f.write('trajin %02d_%02d.nc\n' % (i + 1, j + 1))
        f.write('rms reference mass %s \n' % args.Rm1)
        f.write('average average_%02d_%02d.mdcrd\n' % (i + 1, j + 1))
        f.close()
        os.system('cpptraj <ctraj.AverageAsu.in >../XtalAnalyze.log')


########################################################################
#                                                                      #
#   RMSD calculate of each average structure to pdb and to the other   #
#		average structures                                             #
#                                                                      #
########################################################################
print("Calculating rmsd between all average structures.")
# calculate the rmsd's using ptraj
f = open('ctraj.RmsdPerAsu.in', 'w')
f.write('parm %s\n' % topo)
f.write('reference %s\n' % pdb)
f.write('trajin %s\n' % pdb)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin average_%02d_%02d.mdcrd\n' % (i + 1, j + 1))
f.write('rms reference mass out RMSD_UC_0_0.dat %s time 0.01 \n' % args.Rm1)
f.close()
os.system('cpptraj <ctraj.RmsdPerAsu.in >XtalAnalyze.log')

for i in range(unitcells):
    for j in range(asymunits):
        f = open('ctraj.RmsdPerAsu.in', 'w')
        f.write('parm %s\n' % topo)
        f.write('reference average_%02d_%02d.mdcrd\n' % (i + 1, j + 1))
        f.write('trajin %s\n' % pdb)
        for k in range(unitcells):
            for l in range(asymunits):
                f.write('trajin average_%02d_%02d.mdcrd\n' % (k + 1, l + 1))
        f.write(
            'rms reference mass out RMSD_UC_%02d_%02d.dat %s time 0.01 \n' %
            (i + 1, j + 1, args.Rm1))
        f.close()
        os.system('cpptraj <ctraj.RmsdPerAsu.in >XtalAnalyze.log')

### create matrix combining the rmsd's of each average str to the others ###
mat = zeros((unitcells * asymunits + 1, unitcells * asymunits + 1))
f = genfromtxt('RMSD_UC_0_0.dat')
mat[:, 0] = f[:, 1]

for i in range(unitcells):
    for j in range(asymunits):
        f = genfromtxt('RMSD_UC_%02d_%02d.dat' % (i + 1, j + 1))
        mat[:, (asymunits * i + (j + 1))] = f[:, 1]
savetxt('RMSDUC.dat', mat)

########################################################################
#                                                                      #
#   B-FACTORS  calculate                                               #
#                                                                      #
########################################################################
print("Calculating B-factors for each asu trajectory")
# IRL B-factors
for i in range(unitcells):
    for j in range(asymunits):
        f = open('ctraj.BfactorPerAsu.in', 'w')
        f.write('parm %s\n' % topo)
        f.write('trajin %02d_%02d.nc\n' % (i + 1, j + 1))
        f.write(
            'atomicfluct out bfac_lattice_%02d_%02d.dat %s byatom bfactor\n' %
            (i + 1, j + 1, args.Bm1))
        f.close()
        os.system('cpptraj <ctraj.BfactorPerAsu.in >XtalAnalyze.log')

# I B-factors
for i in range(unitcells):
    for j in range(asymunits):
        f = open('ctraj.BfactorPerAsu.in', 'w')
        f.write('parm %s\n' % topo)
        f.write('trajin %02d_%02d.nc\n' % (i + 1, j + 1))
        f.write('rms first mass %s \n' % args.Rm1)
        f.write(
            'atomicfluct out bfac_monomer_%02d_%02d.dat %s byatom bfactor\n' %
            (i + 1, j + 1, args.Bm1))
        f.close()
        os.system('cpptraj <ctraj.BfactorPerAsu.in >XtalAnalyze.log')

# average of IRL B-factors from each ASU (quasi IR B-factors)
for i in range(unitcells):
    for j in range(asymunits):
        x = genfromtxt('bfac_lattice_%02d_%02d.dat' % (i + 1, j + 1), skip_header=1)
        if i == 0 and j == 0:
            avg_bfacs = x[:]
        else:
            avg_bfacs[:, 1] += x[:, 1]

avg_bfacs[:, 1] /= (unitcells * asymunits)
with open('bfacs_calpha_lattice_INDIVavg.dat', 'w') as f:
    f.write('#Atom       B-factors\n')
    savetxt(f, avg_bfacs, fmt='%7d    %9.4f')
