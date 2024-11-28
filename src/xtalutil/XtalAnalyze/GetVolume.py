#! /usr/bin/env python
from numpy import *
from ReadAmberFiles import *
import argparse
from parmed.utils import netcdf as NetCDF

parser = argparse.ArgumentParser()
parser.add_argument(
    "-r",
    "--SCRestart",
    help="Supercell restart file with correct box info on last line")
parser.add_argument("-t", "--Trajectory", help="Supercell Trajectory")
args = parser.parse_args()

# GET CRYSTAL VOLUME
rst7file = rst7(args.SCRestart)
SCBox = rst7file.Get_Box()
crystvol = Get_volume(SCBox)
angles = SCBox[3:] * 180 / pi
# GET TRAJECTORY VOLUMES
ofile = NetCDF.NetCDFFile(args.Trajectory, 'r', mmap=False)
coords = ofile.variables['coordinates']
frames = coords.shape[0]
lengths = ofile.variables['cell_lengths']
data1 = zeros((frames, 3))
for frame in range(frames):
    box = UCbox = hstack((lengths[frame, :], angles))
    data1[frame, 0] = frame
    current_volume = Get_volume(box)
    data1[frame, 1] = current_volume
    data1[frame, 2] = (current_volume / crystvol) * 100
f = open('volume.dat', 'w')
f.write('Frame       Volume(Ang)    Volume(%_of_exp)\n')
for frame in data1:
    f.write('%8.2f   %12.5f   %12.5f\n' % (frame[0], frame[1], frame[2]))

# CONVERT TO PERCENT
frames = len(data1[:, 1])
minvol = min(data1[:, 1])
maxvol = max(data1[:, 1])
meanvol = mean(data1[:, 1])
minperc = min(data1[:, 2])
maxperc = max(data1[:, 2])
meanperc = mean(data1[:, 2])


# SAVE STATISTICS TO FILE
f = open('volume.txt', 'w')
f.write('%s\n' % (os.getcwd()))
f.write('crystal volume = %10f\n' % crystvol)
f.write('min volume = %10f   %6.2f\n' % (minvol, minperc))
f.write('max volume = %10f   %6.2f\n' % (maxvol, maxperc))
f.write('mean volume = %10f   %6.2f\n' % (meanvol, meanperc))
f.close()
