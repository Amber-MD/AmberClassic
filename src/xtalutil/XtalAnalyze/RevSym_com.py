#!/usr/bin/env python
import os
from numpy import *
from parmed.utils import netcdf as NetCDF
import ReadAmberFiles as raf
import argparse

#======================================================================#
#                                                                      #
# Translate each ASU split trajectory so that it's average center of   #
# mass aligns with the experimental center of mass.                    #
#                                                                      #
#======================================================================#


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--AsuTopology", help="Amber single ASU topology file")
parser.add_argument("-pdb", "--AsuPDB", help="Single ASU PDB file")
parser.add_argument(
    "-pdbs",
    "--PdbWithSymmetry",
    help="PDB file containing SMTRY and CRYST1 records")
parser.add_argument(
    "-prop",
    help="Propagation used to create the supercell. Format: \"x y z\"")
args = parser.parse_args()

############################################
#  SETUP                                   #
############################################
topo = raf.prmtop(args.AsuTopology)
pdb = raf.pdb(args.AsuPDB)
masses = topo.Get_Masses()
xtal_com = raf.COM(pdb.Get_Coords(), masses)
ix = int(args.prop.split()[0])
iy = int(args.prop.split()[1])
iz = int(args.prop.split()[2])
pdbs = raf.pdb(args.PdbWithSymmetry)
smtry, tr = pdbs.Get_SMTRY()
asymunits = len(smtry)

###########################################
#  MAIN                                   #
###########################################

i = 1  # counter for the names of the unitcell trajectories
# for each unit cell in the supercell
for x in range(ix):
    for y in range(iy):
        for z in range(iz):
            for j in range(asymunits):
                # frames stores the number of frames in the trajectory
                orig_filename = 'RevSym_%02d_%02d.nc' % (i, j + 1)
                filename = 'RevSym_com_%02d_%02d.nc' % (i, j + 1)
                print("Translating %s" % filename)
                os.system('cp %s %s' % (orig_filename, filename))
                ofile = NetCDF.NetCDFFile(filename, 'a', mmap=False)
                coords = ofile.variables['coordinates']
                avg = average(coords[:], axis=0)
                sim_com = raf.COM(avg, masses)
                shift = (xtal_com - sim_com).astype(float32)
                coords[:, :, :] = coords[:, :, :] + shift
                ofile.close()
            i += 1
