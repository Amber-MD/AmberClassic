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

data2 = genfromtxt('revsym/rmsd_asu_bkbn_table.dat', skip_header=1)
data3 = genfromtxt('revsym/rmsd_asu_heavy_table.dat', skip_header=1)
data4 = genfromtxt('revsym/rmsd_lat_bkbn_table.dat', skip_header=1)
frames = data2.shape[0]

# padding for tick labels
plt.rc(('xtick.major', 'xtick.minor', 'ytick.major', 'ytick.minor'), pad=15)
# thicker axes frame
plt.rc('axes', linewidth=4)
plt.rc('legend', fontsize=20)

fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)


x = []
y = []
x4errs = []
y4errs = []
yerrs = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data2[10 * i, 0]))
    y.append(average(data2[10 * i:10 * i + 10, 1]))
for i in range(int(frames/150)):
    x4errs.append((data2[150 * i, 0]))
    y4errs.append(average(data2[150 * i:150 * i + 150, 1]))
    yerrs.append(average(data2[150 * i, 2]))
p1 = ax.plot(x, y, 'k', linewidth=2)
ax.errorbar(x4errs, y4errs, yerr=yerrs, fmt='r.', color='k')

x = []
y = []
x4errs = []
y4errs = []
yerrs = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data3[10 * i, 0]))
    y.append(average(data3[10 * i:10 * i + 10, 1]))
for i in range(int(frames/150)):
    x4errs.append((data3[150 * i, 0]))
    y4errs.append(average(data3[150 * i:150 * i + 150, 1]))
    yerrs.append(average(data3[150 * i, 2]))
p2 = ax.plot(x, y, '#9B2DE7', linewidth=2)
ax.errorbar(x4errs, y4errs, yerr=yerrs, fmt='r.', color='#9B2DE7')

x = []
y = []
x4errs = []
y4errs = []
yerrs = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data4[10 * i, 0]))
    y.append(average(data4[10 * i:10 * i + 10, 1]))
for i in range(int(frames/150)):
    x4errs.append((data4[150 * i, 0]))
    y4errs.append(average(data4[150 * i:150 * i + 150, 1]))
    yerrs.append(average(data4[150 * i, 2]))
p3 = ax.plot(x, y, '#37AECD', linewidth=2)
ax.errorbar(x4errs, y4errs, yerr=yerrs, fmt='r.', color='#37AECD')


for label in ax.xaxis.get_ticklabels():
    label.set_fontsize(24)
for label in ax.yaxis.get_ticklabels():
    label.set_fontsize(24)
plt.xlabel('Time (ns)', fontsize=28, labelpad=10)
plt.ylabel('RMSD ($\AA$)', fontsize=28, labelpad=10)
plt.ylim(ymin=0)

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(4)
    line.set_markersize(10)
ax.legend((p1[0], p2[0], p3[0]), ("ASU backbone", "ASU heavy",
                                  "lattice backbone"), bbox_to_anchor=(0.0, 0.0, 1, .24	))

plt.savefig('rmsd_simplemean.pdf')
