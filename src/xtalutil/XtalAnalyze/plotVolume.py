#! /usr/bin/env python
from numpy import *
import matplotlib as mpl
# avoid using any X-windows backend
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from ReadAmberFiles import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-Title", help="Title prefix for plots")
args = parser.parse_args()

# get data and process
data1 = genfromtxt('volume.dat', skip_header=1)
frames = data1.shape[0]
x = []
y = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data1[10 * i, 0]) / 1000)
    y.append(average(data1[10 * i:10 * i + 10, 2]))

# plot settings
plt.rc(('xtick.major', 'xtick.minor', 'ytick.major', 'ytick.minor'), pad=15)
plt.rc('axes', linewidth=4)
plt.rc('legend', fontsize=20)
plt.rc('mathtext', default='regular')

# make plot
fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)
ax.plot(x, y, 'b', linewidth=4)

# more plot settings
plt.ylim((99.5, 100.5))

for label in ax.xaxis.get_ticklabels():
    label.set_fontsize(24)
for label in ax.yaxis.get_ticklabels():
    label.set_fontsize(24)
plt.title(args.Title, fontsize=28)
plt.xlabel('Time (ns)', fontsize=28, labelpad=10)
# ax.set_xticklabels([0,1.0,2.0,3.0,4.0,5.0])
plt.ylabel(r'Volume ($\%$)', fontsize=28, labelpad=10)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(4)
    line.set_markersize(10)

# plt.show()
plt.savefig('volume.pdf')
