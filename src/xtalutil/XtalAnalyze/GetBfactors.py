#! /usr/bin/env python
import sys
import os
from numpy import *

#======================================================================#
#                                                                      #
# Obtain bfactors from a pdb file (usually the original pdb-deposited  #
# file) and print them in the order that they appear in another file   #
# (usually the pdb obtained from an Amber topology. When a topology/crd#
# file is created by tleap, the atom order is rearranged making it     #
# difficult to compare properties such as B-factors. This script       #
# extracts B-factors from the original pdb and prints them in the amber#
# happy order.                                                         #
#                                                                      #
# Arguments:                                                           #
#     [original.pdb] - original pdb file name. Usually the pdb fed to  #
#                      tleap to make amber topology. Atom names and    #
#                      residue order should correspond to Amber. No alt#
#                      confs.
#     [amber.pdb]    - amber-derived pdb file (produced by             #
#                      XtalAnalysis.sh as asu.pdb                      #
#     [nres]         - number of residues to use (ignore all residues  #
#                      such as solvent after this residue seqid)       #
# Output:                                                              #
#   all.bfactors - list of all bfactors in amber order                 #
#   calpha.bfactors - list of only the c-alpha b-factors               #
#   sdch.bfactors - for each residue, mean bfactor of all side-chain   #
#                   atoms. (Glycines skipped)                          #
#   meanresidue.bfactors - for each residue, mean bfactor of all atoms #
#                          (side chain and backbone)                   #
# Usage:                                                               #
#       GetBfactors.py [original.pdb] [amber.pdb] [nres]               #
#                                                                      #
# v1 Pawel Janowski 16.12.2013                                         #
#======================================================================#


gfile = sys.argv[1]
ffile = sys.argv[2]
noofresidues = int(sys.argv[3])


glines = []
with open(gfile) as file_:
    for line in file_:
        if not line[0:4] == 'ATOM':
            continue
        glines.append(line)

flines = []
with open(ffile) as file_:
    for line in file_:
        if not line[0:4] == 'ATOM':
            continue
        flines.append(line)

out = open('all.bfactors', 'w')

lnum = 0

for linef in flines:
    atom = linef[12:16].strip()
    resid = linef[22:26].strip()
    for lineg in glines:
        if lineg[12:16].strip() == atom and lineg[22:26].strip() == resid:
            lnum += 1
            out.write('%4d  %6.2f  %4s  %4s \n' %
                      (lnum, float(lineg[60:66]), atom, resid))
out.close()

########################################################################
# Now produce a file with just C-alpha and C1' b-factors               #
########################################################################

out = open('calpha.bfactors', 'w')
lnum = 0
with open('all.bfactors') as file_:
    for line in file_:
        line = line.strip().split()
        if line[2] == 'CA' or line[2] == 'C1\'':
            lnum += 1
            out.write('%4d  %5.2f  %4s  %4s \n' %
                      (lnum, float(line[1]), line[2], line[3]))
out.close()


########################################################################
# Now produce a file with an avg. side chain b-factor for each residue #
########################################################################

glines = []
with open('all.bfactors') as file_:
    for line in file_:
        line = line.strip().split()
        glines.append(line)
out = open('sdch.bfactors', 'w')
lnum = 0

for i in range(1, noofresidues + 1):
    residue = [x for x in glines if int(x[3]) == i]
    bfacs = array([])
    for atom in residue:
        if atom[2] in ['C', 'O', 'N', 'CA']:
            continue
        elif atom[2][0] == 'H':
            continue
        else:
            bfacs = append(bfacs, float(atom[1]))
    with errstate(invalid='ignore'):
        avgbfac = mean(bfacs)
    if isnan(avgbfac):
        continue
    lnum += 1
    out.write('%4d  %6.2f  %4d\n' % (lnum, avgbfac, i))

out.close()

########################################################################
# Now produce a file with an average b-factor for each entire residue  #
########################################################################


glines = []
with open('all.bfactors') as file_:
    for line in file_:
        line = line.strip().split()
        glines.append(line)
out = open('meanresidue.bfactors', 'w')
lnum = 0

for i in range(1, noofresidues + 1):
    residue = [x for x in glines if int(x[3]) == i]
    bfacs = array([])
    for atom in residue:
        if atom[2][0] == 'H':
            continue
        else:
            bfacs = append(bfacs, float(atom[1]))
    with errstate(invalid='ignore'):
        avgbfac = mean(bfacs)
    if isnan(avgbfac):
        continue
    lnum += 1
    out.write('%4d  %6.2f  %4d \n' % (lnum, avgbfac, i))

out.close()
