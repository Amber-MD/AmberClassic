#! /usr/bin/env python
import os
import numpy as np
from parmed.utils import netcdf as NetCDF
import argparse
import copy

#######################################################################
# Calculate ADP's and B-factors from an MD trajectory.                #
#                                                                     #
# Arguments:                                                          #
#     names of netcdf format MD trajectories                          #
#     -s: skip this interval of frames when reading trajectories      #
#         (default=1)                                                 #
#     -ipdb: if set will replace/create B-factor and ANISOU entries   #
#            this pdb file with the newly calculated values           #
#            (No error checking: make sure atom order in provided pdb #
#             and MD trajectories match.)                             #
#     -opdb: name of output pdb file (default="out.pdb")     		  #
# Example usage:                                                      #
#     xyz2adp.py md1.nc md2.nc -s 10 -ipdb UC.pdb -opdb UC_B.pdb      #
#######################################################################


# get coords array
def get_coords(file, skip_frames):
    ofile = NetCDF.NetCDFFile(file, 'r')
    coords = ofile.variables['coordinates'][::skip_frames, :, :]
    return coords


# calculate adps
def get_adp(coords):

    # DAC: using /usr/bin/python from El Capitan, the following lines
    #   fail; for now, I'll assume we have a reasonable recent version
    #   of numpy....
    #v = np.version.version
    #v = tuple(map(int, (v.split("."))))

    #if v >= (1, 7, 0):
    means = np.mean(coords, axis=0, keepdims=True)
    #else:
    #    tmp_means = np.mean(coords, axis=0)
    #    means = np.zeros((1, tmp_means.shape[0], tmp_means.shape[1]))
    #    means[0, :, :] = tmp_means

    fluct = (coords - means)
    sq_fluct = fluct**2
    u11 = np.mean(sq_fluct[:, :, 0], axis=0) * 10000
    u22 = np.mean(sq_fluct[:, :, 1], axis=0) * 10000
    u33 = np.mean(sq_fluct[:, :, 2], axis=0) * 10000

    u12 = np.mean(fluct[:, :, 0] * fluct[:, :, 1], axis=0) * 10000
    u13 = np.mean(fluct[:, :, 0] * fluct[:, :, 2], axis=0) * 10000
    u23 = np.mean(fluct[:, :, 1] * fluct[:, :, 2], axis=0) * 10000

    B = np.mean(sq_fluct[:, :, 0] + sq_fluct[:, :, 1] + sq_fluct[:, :, 2], axis=0)
    B = B * 8 / 3 * np.pi**2
    return B, u11, u22, u33, u12, u13, u23

# print adps to pdb


def print_adp(B, u11, u22, u33, u12, u13, u23, ipdb, opdb):
    ofile = open(opdb, 'w')
    cnt = 0
    with open(ipdb) as f:
        for line in f:
            line = line.strip()
            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                ofile.write('%s%6.2f%-14s\n' % (line[0:60],B[cnt],line[66:80]))
                ofile.write('ANISOU%s%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f%-10s\n' 
                    % (line[6:28], u11[cnt], u22[cnt], u33[cnt], u12[cnt], 
                    u13[cnt], u23[cnt], line[70:80]))
                cnt += 1
            elif line[0:6] != 'ANISOU':
                ofile.write("%-80s\n" % line)
    f.close()

# main run


def run(files, skip_frames, ipdb, opdb):
    coords = 0
    for file in files:
        try:
            if coords == 0:
                coords = get_coords(file, skip_frames)
        except:
            coords = np.vstack((coords, get_coords(file, skip_frames)))
    B, u11, u22, u33, u12, u13, u23 = get_adp(coords)
    if (ipdb):
        print_adp(B, u11, u22, u33, u12, u13, u23, ipdb, opdb)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "files",
        nargs='*',
        help="names of netcdf trajectory files to get coords")
    parser.add_argument(
        "-s",
        "--skip_frames",
        help="trajectory frame skip step (for big datasets)",
        default=1)
    parser.add_argument("-ipdb", help="if set, will fill the B-factor and ANISOU \
        records of the provided pdb with the calculated values. The pdb \
        must have same number of atoms as the trajectory.")
    parser.add_argument("-opdb", help="output pdb file name", default="out.pdb")
    args = parser.parse_args()
    run(args.files, int(args.skip_frames), args.ipdb, args.opdb)
