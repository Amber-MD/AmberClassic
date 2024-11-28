#! /usr/bin/env python
import os
import sys
from numpy import *
from ReadAmberFiles import *
import argparse

#======================================================================#
#                                                                      #
# Bfactors, RMSD and average structure from Reverse Symmetried ASU     #
# trajectories.                                                        #
#                                                                      #
#======================================================================#


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--AsuTopology", help="Amber single ASU topology file")
parser.add_argument("-pdb", "--AsuPDB", help="Single ASU PDB file")
parser.add_argument("-u", "--UnitCells", help="number of unit cells in crystal")
parser.add_argument("-a", "--ASUs", help="number of asymmetric units per unit cell")
parser.add_argument("-Bm1", help="Mask for Calpha B-factor analysis")
parser.add_argument("-Bm2", help="Mask for side chain B-factor analysis")
parser.add_argument("-Rm1", help="Suffix for heavy atom rmsd")
parser.add_argument("-Rm2", help="Suffix for backbone atom rmsd")
parser.add_argument("-tm", help="Timestep for rmsd output")
args = parser.parse_args()


#############################
# SETUP                     #
#############################
unitcells = int(args.UnitCells)
asymunits = int(args.ASUs)
topo = args.AsuTopology
pdb = args.AsuPDB
traj = nc('RevSym_01_01.nc')
frames = traj.Get_Frames()


######################################
#   Calculate Lattice rmsd vs. time  #
######################################
print("Calculating lattice rmsd.")
f = open('ctraj.LatRmsd.in', 'w')
f.write('parm %s\n' % topo)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin RevSym_%02d_%02d.nc\n' % (i + 1, j + 1))
f.write('reference %s\n' % pdb)
f.write(
    'rms reference mass out rmsd_lat_heavy.dat %s nofit time %s \n' %
    (args.Rm1, args.tm))
f.write(
    'rms reference mass out rmsd_lat_bkbn.dat %s nofit time %s \n' %
    (args.Rm2, args.tm))
f.close()
os.system('cpptraj <ctraj.LatRmsd.in >>../XtalAnalyze.log')


##################################
#   Calculate ASU rmsd  vs. time #
##################################
print("Calculating asu rmsd.")
f = open('ctraj.ASURmsd.in', 'w')
f.write('parm %s\n' % topo)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin RevSym_%02d_%02d.nc\n' % (i + 1, j + 1))
f.write('reference %s\n' % pdb)
f.write('rms reference mass out rmsd_asu_heavy.dat %s time %s \n' % (args.Rm1, args.tm))
f.write('rms reference mass out rmsd_asu_bkbn.dat %s time %s \n' % (args.Rm2, args.tm))
f.close()
os.system('cpptraj <ctraj.ASURmsd.in >>../XtalAnalyze.log')


##################################
#   Process rmsd data            #
##################################
print("Processing rmsd data.")
#  This function averages the msd from each ASU into a single value


def ASUrmsd(ifile, ofile, frames, unitcells):
    f = open(ofile, 'w')
    data = genfromtxt(ifile)
    x = array([])
    for i in range(frames):
        for j in range(unitcells):
            x = append(x, data[i + j * frames, 1])
        rmsd = sqrt(mean(x**2))
        f.write('%7.2f  %7.5f \n' % (data[i, 0], rmsd))
        x = array([])
    f.close()
ASUrmsd('rmsd_lat_heavy.dat', 'rmsd_lat_heavy_ASU.dat', frames, unitcells * asymunits)
ASUrmsd('rmsd_asu_heavy.dat', 'rmsd_asu_heavy_ASU.dat', frames, unitcells * asymunits)
ASUrmsd('rmsd_lat_bkbn.dat', 'rmsd_lat_bkbn_ASU.dat', frames, unitcells * asymunits)
ASUrmsd('rmsd_asu_bkbn.dat', 'rmsd_asu_bkbn_ASU.dat', frames, unitcells * asymunits)

#  This function processes the rmsd's into a single table.
#    Ncolumns: 3 + the number of asus
#    Nrows:    1 for each frame of the trajectory
#    Col1: time
#    Col2: mean rmsd over all asu's at given time point
#    Col3: std of rmsd over all asu's at given time point
#    Col4+: rmsd of each asu at given time point


def ProcessRmsd(ifile, ofile, frames, total_monomers):
    f = open(ofile, 'w')
    data = genfromtxt(ifile)
    assert len(data) % total_monomers == 0, "number of frames not a multiple of total_monomers"
    x = array([])
    f.write('    time     mean      std     rmsds_of_each_monomer\n')
    for i in range(frames):
        for j in range(total_monomers):
            x = append(x, data[i + j * frames, 1])
        f.write('%7.2f   ' % data[i, 0])
        f.write('%7.3f   %7.3f   ' % (mean(x), std(x, ddof=1)))
        # To take mean of msd instad of mean of rmsd:
        # f.write('%7.3f   %7.3f   ' %(sqrt(mean(x**2)), std(x,ddof=1)))
        for value in x:
            f.write('%7.3f   ' % value)
        f.write('\n')
        x = array([])
    f.close()
ProcessRmsd(
    'rmsd_lat_heavy.dat',
    'rmsd_lat_heavy_table.dat',
    frames,
    unitcells *
    asymunits)
ProcessRmsd(
    'rmsd_asu_heavy.dat',
    'rmsd_asu_heavy_table.dat',
    frames,
    unitcells *
    asymunits)
ProcessRmsd('rmsd_lat_bkbn.dat', 'rmsd_lat_bkbn_table.dat', frames, unitcells * asymunits)
ProcessRmsd('rmsd_asu_bkbn.dat', 'rmsd_asu_bkbn_table.dat', frames, unitcells * asymunits)


#######################################################
# Calculate asu average structure and its rmsd        #
#######################################################
print("Calculating average asu-rmsd structure.")
f = open('ctraj.AvgCoord_asu.in', 'w')
f.write('parm %s\n' % topo)
f.write('reference %s\n' % pdb)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin RevSym_%02d_%02d.nc\n' % (i + 1, j + 1))
f.write('rms reference mass %s\n' % args.Rm2)
f.write('average AvgCoord_asu.rst7 restart \n')
f.close()
os.system('cpptraj <ctraj.AvgCoord_asu.in >>../XtalAnalyze.log')

f = open('ctraj.rmsd_asu.in', 'w')
f.write('parm %s\n' % topo)
f.write('reference %s\n' % pdb)
f.write('trajin AvgCoord_asu.rst7\n')
f.write(
    'rms reference mass out AvgCoord_asu_rmsd_heavy.dat %s time %s \n' %
    (args.Rm1, args.tm))
f.write(
    'rms reference mass out AvgCoord_asu_rmsd_bkbn.dat  %s time %s \n' %
    (args.Rm2, args.tm))
f.close()
os.system('cpptraj <ctraj.rmsd_asu.in >>../XtalAnalyze.log')


#######################################################
# Calculate lattice average structure and its rmsd    #
#######################################################
print("Calculating average lattice-rmsd structure.")
f = open('ctraj.AvgCoord_lat.in', 'w')
f.write('parm %s\n' % topo)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin RevSym_%02d_%02d.nc\n' % (i + 1, j + 1))
f.write('average AvgCoord_lat.rst7 restart \n')
f.close()
os.system('cpptraj <ctraj.AvgCoord_lat.in >>../XtalAnalyze.log')

f = open('ctraj.rmsd_lat.in', 'w')
f.write('parm %s\n' % topo)
f.write('reference %s\n' % pdb)
f.write('trajin AvgCoord_lat.rst7\n')
f.write(
    'rms reference mass nofit out AvgCoord_lat_rmsd_heavy.dat %s time %s \n' %
    (args.Rm1, args.tm))
f.write(
    'rms reference mass nofit out AvgCoord_lat_rmsd_bkbn.dat  %s time %s \n' %
    (args.Rm2, args.tm))
f.close()
os.system('cpptraj <ctraj.rmsd_lat.in >>../XtalAnalyze.log')


################################
# Calculate I+R+L B-factors    #
################################
print("Calculating \'Internal+RigidBody+Lattice\' B-factors")
f = open('ctraj.bfactor_lat.in', 'w')
f.write('parm %s\n' % topo)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin RevSym_%02d_%02d.nc\n' % (i + 1, j + 1))
f.write('atomicfluct out bfac_lat_calpha.dat %s byatom bfactor\n' % (args.Bm1))
f.write('atomicfluct out bfac_lat_sdch.dat %s byres bfactor\n' % (args.Bm2))
f.close()
os.system('cpptraj <ctraj.bfactor_lat.in >>../XtalAnalyze.log')
#~ os.unlink('ctraj_bfactor_lat')


#############################
#   Calculate I+R B-factors #
#############################
print("Calculating \'Internal+RigidBody\' B-factors")
f = open('ctraj.bfactor_com.in', 'w')
f.write('parm %s\n' % topo)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin RevSym_com_%02d_%02d.nc\n' % (i + 1, j + 1))
f.write('atomicfluct out bfac_com_calpha.dat %s byatom bfactor\n' % (args.Bm1))
f.write('atomicfluct out bfac_com_sdch.dat %s byres bfactor\n' % (args.Bm2))
f.close()
os.system('cpptraj <ctraj.bfactor_com.in >>../XtalAnalyze.log')
#~ os.unlink('ctraj_bfactor_lat')


###########################
#   Calculate I B-factors #
###########################
print("Calculating \'Internal motion\' B-factors")
f = open('ctraj.bfactor_asu.in', 'w')
f.write('parm %s\n' % topo)
for i in range(unitcells):
    for j in range(asymunits):
        f.write('trajin RevSym_%02d_%02d.nc\n' % (i + 1, j + 1))
f.write('reference AvgCoord_asu.rst7\n')
f.write('rms reference mass %s\n' % args.Rm1)
f.write('atomicfluct out bfac_asu_calpha.dat %s byatom bfactor\n' % (args.Bm1))
f.write('atomicfluct out bfac_asu_sdch.dat %s  byres bfactor\n' % (args.Bm2))
f.close()
os.system('cpptraj <ctraj.bfactor_asu.in >>../XtalAnalyze.log')
#~ os.unlink('ctraj_bfactor_asu')
