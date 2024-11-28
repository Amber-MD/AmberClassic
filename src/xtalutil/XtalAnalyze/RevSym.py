#!/usr/bin/env python
import os
from numpy import *
from parmed.utils import netcdf as NetCDF
from ReadAmberFiles import *
import argparse

#======================================================================#
#                                                                      #
# Reverse symmetry operations on each ASU split trajectory             #
#                                                                      #
#======================================================================#


parser = argparse.ArgumentParser()
parser.add_argument(
    "-p",
    "--PdbWithSymmetry",
    help="PDB file containing SMTRY and CRYST1 records")
parser.add_argument(
    "-prop",
    help="Propagation used to create the supercell. Format: \"x y z\"")
parser.add_argument(
    "-r",
    "--SCRestart",
    help="Supercell restart file with correct box info on last line")
args = parser.parse_args()


############################################
#  SETUP                                   #
############################################
# Get propagations
ix = int(args.prop.split()[0])
iy = int(args.prop.split()[1])
iz = int(args.prop.split()[2])
# Get box info
rst7file = rst7(args.SCRestart)
SCBox = rst7file.Get_Box()
expbox = hstack((SCBox[0:3] / [ix, iy, iz], SCBox[3:]))
exp_u, exp_invu = CompXfrm(expbox)
# Get symmetry rotation matrices (smtry) and translation vectors (tr)
pdbfile = pdb(args.PdbWithSymmetry)
smtry, tr = pdbfile.Get_SMTRY()
asymunits = len(smtry)
# Sanity check
print('Number of symmetry operations in pdb: %d' % (asymunits))
print('Unit cell propagation was %dx%dx%d' % (ix, iy, iz))
print('Total number of unit cells: %d' % (ix * iy * iz))
print('Supercell box dimensions are: ')
print(SCBox)

###########################################
#  MAIN                                   #
###########################################

i = 1  # counter for the names of the unitcell trajectories
# for each unit cell in the supercell
for x in range(ix):
    for y in range(iy):
        for z in range(iz):

            # Move vector in fractional coordinates
            FracVector = array([x, y, z], dtype=float32)
            exp_MoveVector = dot(exp_invu, FracVector).astype(float32)

            for j in range(asymunits):
                # get symmetry operations for this asym unit
                s = float32(transpose(smtry[j]))
                t = float32(tr[j])
                print('\n\nStarting file %02d_%02d.nc' % (i, j + 1))
                print('Symmetry rotation matrix is:')
                print(s)
                print('Symmetry translation vector is %s' % (str(t)))
                print("Unit cell translation vector is %s" % (str(exp_MoveVector)))
                # copy the unit cell trajectory to a new taget file
                # read the netcdf file of each unitcell trajectory
                # coords stores the entire coordinates variable which is
                #    [frames,atoms,xyzcoords]
                # frames stores the number of frames in the trajectory
                orig_filename = '../splittrajectories/%02d_%02d.nc' % (i, j + 1)
                filename = 'RevSym_%02d_%02d.nc' % (i, j + 1)
                os.system('cp %s %s' % (orig_filename, filename))
                ofile = NetCDF.NetCDFFile(filename, 'a', mmap=False)
                coords = ofile.variables['coordinates']
                frames = coords.shape[0]
                celllen = ofile.variables['cell_lengths']

                # iterate over each frame: calculate move scaled to frame
                # box size and apply reverse move
                for frame in range(frames):
                    # get frame box
                    box = celllen[frame, :]
                    # convert to unit cell box
                    box = box / [ix, iy, iz]
                    # add angle information as provided by user
                    UCbox = hstack((box, SCBox[3:]))
                    # calc orthogonalization matrix
                    u, invu = CompXfrm(UCbox)
                    # matrix product
                    MoveVector = dot(invu, FracVector).astype(float32)
                    coords[frame, :, :] = coords[frame, :, :] - MoveVector
                # reverse symmetry operate to original asym unit by
                # applying the symmetry translation and rotation
                coords[:, :, :] = dot((coords[:, :, :] - t), linalg.inv(s))

                # close file (write to netcdf file)
                ofile.close()
            i += 1
