#! /usr/bin/env python
import sys
import os
from numpy import *
import matplotlib as mpl
# avoid using any X-windows backend
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


### SET VARIABLES  ###
rows = sys.argv[4]
columns = sys.argv[5]


####################################
# set up general variables         #
####################################
title_txt = ('%s rmsd' % sys.argv[1])
totasymunits = int(sys.argv[2]) * int(sys.argv[3])
data1 = genfromtxt('revsym/rmsd_lat_bkbn.dat')
data2 = genfromtxt('revsym/rmsd_asu_bkbn.dat')
frames = data1.shape[0] / totasymunits

x = []
for i in range(int(frames/10)):  # this is to get average over each 10 frames
    x.append((data1[10 * i, 0]))

#####################################
# GENERAL PLOT SETTINGS             #
#####################################
plt.rc(('xtick.major', 'xtick.minor', 'ytick.major', 'ytick.minor'), pad=15)
plt.rc('axes', linewidth=4)
plt.rc('legend', fontsize=20)
colors = ['#F86606',
          '#CC0000',
          '#F3EF02',
          '#5DF304',
          '#4E9A06',
          '#C4A000',
          '#729FCF',
          '#0618F4',
          '#06EFF4',
          '#A406F4',
          '#F4069D',
          '#936F70']


########################################################################
#                                                                      #
#       PLOT Lattice RMSD per monomer on oneplot                       #
#                                                                      #
########################################################################
fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)
for unit in range(totasymunits):
    dataunit = data1[unit * frames:unit * frames + frames]
    y = []
    for i in range(int(frames/10)):  # this is to get average over each 10 frames
        y.append(average(dataunit[10 * i:10 * i + 10, 1]))
    ax.plot(x, y, colors[unit % 12], linewidth=2, label=str(unit + 1))
for label in ax.xaxis.get_ticklabels():
    label.set_fontsize(24)
for label in ax.yaxis.get_ticklabels():
    label.set_fontsize(24)
plt.title('%s per monomer lattice rmsd' % title_txt, fontsize=28)
plt.xlabel('Time (ns)', fontsize=28, labelpad=10)
# ax.set_xticklabels([0,1.0,2.0,3.0,4.0,5.0])
plt.ylabel(r"RMSD ($\AA$)", fontsize=28, labelpad=10)
#~ plt.ylim((0,6))
# plt.xlim(xmax=51)	#modify to trajectory length
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(4)
    line.set_markersize(10)
from matplotlib.font_manager import fontManager, FontProperties
font = FontProperties(size=18)
ax.legend()
ax.legend(bbox_to_anchor=(0, 0, .95, .2), prop=font, ncol=4)
plt.savefig('rmsd_permonomer_lattice.pdf')

########################################################################
#                                                                      #
#       PLOT Lattice RMSD per monomer on separate plots                #
#                                                                      #
########################################################################
fig = plt.figure(figsize=(12, 8))
fig.subplots_adjust(wspace=.15, right=.95, left=.05, top=.93, bottom=.05, hspace=.35)
plt.suptitle('%s per monomer lattice rmsd' % title_txt, fontsize=28)
for unit in range(totasymunits):
    ax = plt.subplot(rows, columns, unit + 1)

    dataunit = data1[unit * frames:unit * frames + frames]
    y = []
    for i in range(int(frames/10)):  # this is to get average over 10 frames
        y.append(average(dataunit[10 * i:10 * i + 10, 1]))
    ax.plot(x, y, colors[unit % 12], linewidth=1, label=str(unit + 1))

    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(10)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(10)
    ax.xaxis.labelpad = 0
    ax.yaxis.labelpad = 0
    plt.title(str(unit + 1), fontsize=12)
    plt.xlabel('Time(ns)', fontsize=12)
    plt.ylabel(r"RMSD ($\AA$)", fontsize=12)
    minorLocator = MultipleLocator(5)
    # ax.xaxis.set_minor_locator(minorLocator)
    #~ plt.ylim(0,5)
    #~ plt.xlim(0,120)
#~ plt.legend(["Simulation", "Experiment"], loc='upper right', bbox_to_anchor=(0, 0, 1, 1),bbox_transform=plt.gcf().transFigure, prop=dict(size='x-small'))
# plt.show()
plt.savefig('rmsd_permonomer_lattice_INDIV.pdf')


########################################################################
#                                                                      #
#       PLOT monomer RMSD per monomer on one plot                      #
#                                                                      #
########################################################################
fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)
for unit in range(totasymunits):
    dataunit = data2[unit * frames:unit * frames + frames]
    y = []
    for i in range(int(frames/10)):  # this is to get average over 10 frames
        y.append(average(dataunit[10 * i:10 * i + 10, 1]))
    ax.plot(x, y, colors[unit % 12], linewidth=2, label=str(unit + 1))
for label in ax.xaxis.get_ticklabels():
    label.set_fontsize(24)
for label in ax.yaxis.get_ticklabels():
    label.set_fontsize(24)
plt.title('%s per monomer monomer rmsd' % title_txt, fontsize=28)
plt.xlabel('Time (ns)', fontsize=28, labelpad=10)
# ax.set_xticklabels([0,1.0,2.0,3.0,4.0,5.0])
plt.ylabel(r"RMSD ($\AA$)", fontsize=28, labelpad=10)
#~ plt.ylim((0,6))
# plt.xlim(xmax=51)	#modify to trajectory length
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(4)
    line.set_markersize(10)
from matplotlib.font_manager import fontManager, FontProperties
font = FontProperties(size=18)
ax.legend()
ax.legend(bbox_to_anchor=(0, 0, .95, .2), prop=font, ncol=4)
plt.savefig('rmsd_permonomer_monomer.pdf')

########################################################################
#                                                                      #
#       PLOT monomer RMSD per monomer on separate plots                #
#                                                                      #
########################################################################
fig = plt.figure(figsize=(12, 8))
fig.subplots_adjust(wspace=.15, right=.95, left=.05, top=.93, bottom=.05, hspace=.35)
plt.suptitle('%s per monomer monomer rmsd' % title_txt, fontsize=28)
for unit in range(totasymunits):
    ax = plt.subplot(rows, columns, unit + 1)

    dataunit = data2[unit * frames:unit * frames + frames]
    y = []
    for i in range(int(frames/10)):  # this is to get average over 10 frames
        y.append(average(dataunit[10 * i:10 * i + 10, 1]))
    ax.plot(x, y, colors[unit % 12], linewidth=1, label=str(unit + 1))

    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(10)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(10)
    ax.xaxis.labelpad = 0
    ax.yaxis.labelpad = 0
    plt.title(str(unit + 1), fontsize=12)
    plt.xlabel('Time(ns)', fontsize=12)
    plt.ylabel(r"RMSD ($\AA$)", fontsize=12)
    minorLocator = MultipleLocator(5)
    #~ ax.xaxis.set_minor_locator(minorLocator)
    #~ plt.ylim(0,2)
    #~ plt.xlim(0,120)
#~ plt.legend(["Simulation", "Experiment"], loc='upper right', bbox_to_anchor=(0, 0, 1, 1),bbox_transform=plt.gcf().transFigure, prop=dict(size='x-small'))
# plt.show()
plt.savefig('rmsd_permonomer_monomer_INDIV.pdf')
