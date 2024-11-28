#! /usr/bin/env python
import sys
import os
from numpy import *
import matplotlib as mpl
# avoid using any X-windows backend
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as patches
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-u", "--UnitCells", help="number of unit cells in crystal")
parser.add_argument("-a", "--ASUs", help="number of asymmetric units per unit cell")
parser.add_argument("-Title", help="Title prefix for plots")
parser.add_argument(
    "-b",
    "--CrystBfac",
    help="File with experimental B-factors for plotting")
parser.add_argument("-row", help="Rows of panels")
parser.add_argument("-col", help="Columns of panels")
args = parser.parse_args()


###########################
# SETUP                   #
###########################
asymunits = int(args.ASUs)
unitcells = int(args.UnitCells)
title_txt = args.Title
crystbfac = args.CrystBfac
rows = args.row
columns = args.col

########################################################################
#                                                                      #
#   PLOT B-FACTORS on individual plots                                 #
#                                                                      #
########################################################################
fig = plt.figure(figsize=(12, 8))
fig.subplots_adjust(wspace=.15, right=.95, left=.05, top=.93, bottom=.05, hspace=.35)
plt.suptitle('%s per monomer B-factors' % title_txt, fontsize=28)
data1 = genfromtxt(crystbfac)
for i in range(unitcells):
    for j in range(asymunits):
        ax = plt.subplot(rows, columns, i * asymunits + j + 1)
        lines1 = ax.plot(data1[:, 0], data1[:, 1], 'k-')

        data2 = genfromtxt('bfac_monomer_%02d_%02d.dat' % (i + 1, j + 1))
        lines2 = ax.plot(data1[:, 0], data2[:, 1], 'r-')

        data3 = genfromtxt('bfac_lattice_%02d_%02d.dat' % (i + 1, j + 1))
        lines3 = ax.plot(data1[:, 0], data3[:, 1], 'b-')

        for label in ax.xaxis.get_ticklabels():
            label.set_fontsize(4)
        for label in ax.yaxis.get_ticklabels():
            label.set_fontsize(4)
        ax.xaxis.labelpad = 0
        ax.yaxis.labelpad = 0
        plt.title(str(i * j + 1), fontsize=6)
        plt.xlabel('Time(ns)', fontsize=6)
        plt.ylabel('B-factor', fontsize=6)
        minorLocator = MultipleLocator(5)
        ax.xaxis.set_minor_locator(minorLocator)
        plt.ylim(0, 100)
plt.legend(["Experimental", "RMSD method", "RevSym method"], loc='upper right', bbox_to_anchor=(
    0, 0, 1, 1), bbox_transform=plt.gcf().transFigure, prop=dict(size='x-small'))
plt.savefig('../bfac_INDIV.pdf')


########################################################################
#                                                                      #
#   PLOT INDIVIDUAL ASYMUNIT RMSD's                                    #
#                                                                      #
########################################################################
### This function rescales the colorbar to be centered at zero. ###
def cmap_center_at_zero(cmap, array):
    array_range = min(array), max(array)
    center = 0.
    if not ((array_range[0] < center) and (center < array_range[1])):
        return cmap
    center_ratio = abs(center - array_range[0]) / abs(array_range[1] - array_range[0])
    if not (0. < center_ratio) & (center_ratio < 1.):
        return cmap
    a = math.log(center_ratio) / math.log(0.5)
    if a < 0.:
        return cmap
    cdict = copy.copy(cmap._segmentdata)
    fn = lambda x: (x[0]**a, x[1], x[2])
    for key in ('red', 'green', 'blue'):
        cdict[key] = list(map(fn, cdict[key]))
        cdict[key].sort()
        assert (cdict[key][0] < 0 or cdict[key][-1] >
                1), "Resulting indices extend out of the [0, 1] segment."
    return matplotlib.colors.LinearSegmentedColormap('colormap', cdict, 1024)

### Just a 2D color plot matrix ###
fig = plt.figure()
plt.suptitle(title_txt, fontsize=18)
ax = plt.subplot(111)
dist = genfromtxt('RMSDUC.dat')
plt.imshow(dist, cmap=plt.cm.PuOr, origin='lower', interpolation='nearest')
plt.title('RMSD of all monomer average structures (crystal is 0)', fontsize=14)
plt.colorbar()
majorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_locator(majorLocator)
plt.savefig('../RmsdPerAsu_matrix.pdf')


### Matrix, bar graph and min/max/mean field ###
fig = plt.figure()
plt.suptitle(title_txt, fontsize=18)

# 2D color plot matrix
#~ majorLocator   = MultipleLocator(16)
#~ ax.yaxis.set_major_locator(majorLocator)
dist = genfromtxt('RMSDUC.dat')
ax1 = plt.subplot(211)
majorLocator = MultipleLocator(1)
ax1.xaxis.set_major_locator(majorLocator)
#~ ax1.title.set_y(.975)
plt.imshow(dist, cmap=plt.cm.PuOr, origin='lower', interpolation='nearest')
plt.title('rmsd matrix (crystal is 0)', fontsize=10)
plt.colorbar()

# Create bar graph
y = 'RMSD_UC_0_0.dat'
dist = genfromtxt(y)
ax2 = plt.subplot(212)
plt.bar(dist[:, 0] * 100 - 1, dist[:, 1], align='center')
#~ ax2.title.set_y(.975)
majorLocator = MultipleLocator(1)
ax2.xaxis.set_major_locator(majorLocator)
ax2.yaxis.set_major_locator(majorLocator)
#~ ax2.set_xticklabels([1,2,3,4,5,6,7,8,9,10,11,12])
plt.title('rmsd of each unit cell to crystal', fontsize=10)
plt.xlim(0.5, (int(asymunits) * int(unitcells) + .5))

# Create text giving mean, max, min of the matrix
mmin = dist[1:, 1].min()
mmax = dist[1:, 1].max()
mmean = dist[1:, 1].mean()
#~ txt=plt.figtext(0.15,0.85, 'max = 1.6236 \n\nmin = 0.7028 \n\nmean = 0.9342', backgroundcolor='white', verticalalignment='top')
ax3 = fig.add_axes([0, 0, 1, 1])
p = patches.Rectangle((0.1, .65), .3, .25, fill=True, color='white', clip_on=True)
ax3.add_patch(p)
ax3.set_axis_off()
ax3.text(.15, .7, 'max = ' + str(mmax) + ' \n\nmin = ' +
         str(mmin) + ' \n\nmean = ' + str(mmean))

#~ plt.show()
plt.savefig('../RmsdPerAsu_bar_matrix.pdf')
