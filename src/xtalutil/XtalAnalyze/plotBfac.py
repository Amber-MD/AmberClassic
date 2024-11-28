#! /usr/bin/env python
import sys
import os
from numpy import *
import matplotlib as mpl
# avoid using any X-windows backend
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#######################
# SET Input variables #
#######################

x = [
    sys.argv[2],
    'revsym/bfac_asu_calpha.dat',
    'revsym/bfac_com_calpha.dat',
    'revsym/bfac_lat_calpha.dat',
    sys.argv[3],
    'revsym/bfac_asu_sdch.dat',
    'revsym/bfac_com_sdch.dat',
    'revsym/bfac_lat_sdch.dat'
]
figname = 'bfac_calpha.pdf'

###########################
# Matlplotlib rc settings #
###########################
# padding for tick labels
plt.rc(('xtick.major', 'xtick.minor', 'ytick.major', 'ytick.minor'), pad=10)
# thicker axes frame
plt.rc('axes', linewidth=4)
plt.rc('legend', fontsize=20)
# This will set width of all tickmarks. I can then reset the major ticks
# so this is specifically for setting the minor tick width.
plt.rc('lines', markeredgewidth=2)
plt.rc('xtick.minor', size=9)
plt.rc('lines', linewidth=4)


###############
# Import data #
###############
data1 = genfromtxt(x[0])
data2 = genfromtxt(x[1])
data3 = genfromtxt(x[2])
data4 = genfromtxt(x[3])
# data5=genfromtxt(x[4])
# data6=genfromtxt(x[5])

########
# Plot #
########
fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)
lines1 = ax.plot(data1[:, 0], data1[:, 1], 'k-')
lines2 = ax.plot(data1[:, 0], data2[:, 1], 'b-')
lines3 = ax.plot(data1[:, 0], data3[:, 1], 'g-')
lines4 = ax.plot(data1[:, 0], data4[:, 1], 'r-')
# lines3=ax.plot(data1[:,0], data3[:,1],'y-')


##########
# Labels #
##########
#plt.xlabel('Atom',fontsize=28, labelpad=10)
plt.ylabel('B-factor', fontsize=28, labelpad=5)
st = fig.suptitle(r"%s Bfactors for C$_\alpha$" % sys.argv[1], fontsize=28)
# plt.title('',fontsize=28)

########
# Axis #
########
ax.autoscale(enable=True, axis='x', tight='True')
# plt.ylim(0,20)
# plt.xlim(0,200)
# ticks locations
minorLocator = MultipleLocator(1)
majorLocator = MultipleLocator(10)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_minor_locator(minorLocator)
# tick label
for label in ax.xaxis.get_ticklabels():
    label.set_fontsize(20)
for label in ax.yaxis.get_ticklabels():
    label.set_fontsize(20)
#~ ax.set_xticklabels([0,0.5,1.0,1.5,2.0,2.5])

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
# tick sizes, for minor labels set rc file above
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(4)
    line.set_markersize(12)

############################################
# LEGEND
############################################
from matplotlib.font_manager import fontManager, FontProperties
font = FontProperties(size=18)
ax.legend([r"Experimental", r"I Bfac", r"I+R Bfac", r"I+R+L Bfac"],
          bbox_to_anchor=(0, 0.0, .95, .95), prop=font)

############################################
# Annotations
#############################################
#~ plt.annotate('Phe4',xy=(26,18), xytext=(26,20),arrowprops=dict(facecolor='black',shrink=.1,width=1,headwidth=5))
#~ plt.annotate('Slope=$7.2 x 10^{-4}$',xy=(5000,200), xytext=(340000,200),fontsize=20)
#plt.grid(which='both',color='b', linestyle='--', linewidth=.5)


#~ plt.show()
plt.savefig(figname)

#=========================================================================

figname = 'bfac_residue.pdf'

#############################
# Matlplotlib rc settings
#############################
# padding for tick labels
plt.rc(('xtick.major', 'xtick.minor', 'ytick.major', 'ytick.minor'), pad=10)
# thicker axes frame
plt.rc('axes', linewidth=4)
plt.rc('legend', fontsize=20)
# This will set width of all tickmarks. I can then reset the major ticks
# so this is specifically for setting the minor tick width.
plt.rc('lines', markeredgewidth=2)
plt.rc('xtick.minor', size=9)
plt.rc('lines', linewidth=4)


##############################
# Import data
##############################
#~ data1=genfromtxt(x[0])
#~ data2=genfromtxt(x[1])
#~ data3=genfromtxt(x[2])
data5 = genfromtxt(x[4])
data6 = genfromtxt(x[5])
data7 = genfromtxt(x[6])

###############################
# Plot
##############################
fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)
#~ lines1=ax.plot(data1[:,0], data1[:,1],'k-')
#~ lines2=ax.plot(data1[:,0], data2[:,1],'b-')
#~ lines3=ax.plot(data1[:,0], data3[:,1],'#2FF101')
lines5 = ax.plot(data5[:, 0], data5[:, 1], 'm-')
lines6 = ax.plot(data5[:, 0], data6[:, 1], 'r-')
lines7 = ax.plot(data5[:, 0], data7[:, 1], 'c-')

# Set line width already set in rc above
# plt.setp(lines1,linewidth=4)
# plt.setp(lines2,linewidth=4)


###############################
# Labels
################################
#plt.xlabel('Atom',fontsize=28, labelpad=10)
plt.ylabel('B-factor', fontsize=28, labelpad=5)
fig.suptitle(r"%s Sidechain B-factors " % sys.argv[1], fontsize=28)
# plt.title('',fontsize=28)

############################
# Axis
############################
# scale
ax.autoscale(enable=True, axis='x', tight='True')
#~ plt.ylim(0,3)
#~ plt.xlim(0,200)
# ticks locations
minorLocator = MultipleLocator(1)
majorLocator = MultipleLocator(10)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_minor_locator(minorLocator)
# tick label
for label in ax.xaxis.get_ticklabels():
    label.set_fontsize(20)
for label in ax.yaxis.get_ticklabels():
    label.set_fontsize(20)
#~ ax.set_xticklabels([0,0.5,1.0,1.5,2.0,2.5])

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
# tick sizes, for minor labels set rc file above
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(4)
    line.set_markersize(12)

############################################
# LEGEND
############################################
from matplotlib.font_manager import fontManager, FontProperties
font = FontProperties(size=18)
ax.legend([r"crystal", r"rmsd", r"revsym"], bbox_to_anchor=(0, 0.0, .95, .95), prop=font)

############################################
# Annotations
#############################################
#~ plt.annotate('Phe4',xy=(26,18), xytext=(26,20),arrowprops=dict(facecolor='black',shrink=.1,width=1,headwidth=5))
#~ plt.annotate('Slope=$7.2 x 10^{-4}$',xy=(5000,200), xytext=(340000,200),fontsize=20)
#plt.grid(which='both',color='b', linestyle='--', linewidth=.5)

#~ plt.show()
plt.savefig(figname)
