#! /usr/bin/env python
import sys
import os
from numpy import *
import matplotlib as mpl
# avoid using any X-windows backend
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

title_txt = ('%s rmsd' % sys.argv[1])

data1 = genfromtxt('revsym/rmsd_lat_bkbn_ASU.dat')
data2 = genfromtxt('revsym/rmsd_asu_bkbn_ASU.dat')
data3 = genfromtxt('revsym/rmsd_asu_heavy_ASU.dat')

frames = data1.shape[0]

# padding for tick labels
plt.rc(('xtick.major', 'xtick.minor', 'ytick.major', 'ytick.minor'), pad=15)
# thicker axes frame
plt.rc('axes', linewidth=4)
plt.rc('legend', fontsize=20)

fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)


x = []
y = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data1[10 * i, 0]))
    y.append(average(data1[10 * i:10 * i + 10, 1]))
ax.plot(x, y, '#F86606', linewidth=4)

x = []
y = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data2[10 * i, 0]))
    y.append(average(data2[10 * i:10 * i + 10, 1]))
ax.plot(x, y, 'b', linewidth=4)

x = []
y = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data3[10 * i, 0]))
    y.append(average(data3[10 * i:10 * i + 10, 1]))
ax.plot(x, y, '#9B2DE7', linewidth=4)

for label in ax.xaxis.get_ticklabels():
    label.set_fontsize(24)
for label in ax.yaxis.get_ticklabels():
    label.set_fontsize(24)
plt.title(title_txt, fontsize=28)
plt.xlabel('Time (ns)', fontsize=28, labelpad=10)
plt.ylabel(r"RMSD ($\AA$)", fontsize=28, labelpad=10)

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(4)
    line.set_markersize(10)

from matplotlib.font_manager import fontManager, FontProperties
font = FontProperties(size=18)
ax.legend(["lattice backbone", "asu backbone", "asu all heavy"],
          bbox_to_anchor=(0, 0, .95, .2), prop=font)

plt.savefig('rmsd.pdf')
