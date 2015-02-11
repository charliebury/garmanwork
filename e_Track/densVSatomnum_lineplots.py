# -*- coding: utf-8 -*-

# this set of functions plots a line graph of density metric over the whole 
# range of atom numbers

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker
import matplotlib as mpl
import sys
import math

# this function is used in the plotting function 'damage_wrtatmnum' below
def annotate_group(name, xspan, ax=None):
    """Annotates a span of the x-axis"""
    def annotate(ax, name, left, right, y, pad):
        arrow = ax.annotate(name,
                xy=(left, y), xycoords='data',
                xytext=(right, y-pad), textcoords='data',
                annotation_clip=False, verticalalignment='top',
                horizontalalignment='center', linespacing=2.0,
                arrowprops=dict(arrowstyle='-', shrinkA=0, shrinkB=0,
                        connectionstyle='angle,angleB=90,angleA=0,rad=5')
                )
        return arrow
    if ax is None:
        ax = plt.gca()
    ymin = ax.get_ylim()[0]
    ypad = 0.01 * np.ptp(ax.get_ylim())
    xcenter = np.mean(xspan)
    left_arrow = annotate(ax, name, xspan[0], xcenter, ymin, ypad)
    right_arrow = annotate(ax, name, xspan[1], xcenter, ymin, ypad)
    return left_arrow, right_arrow

# this function is used in the plotting function 'damage_wrtatmnum' below
def make_second_bottom_spine(ax=None, label=None, offset=0, labeloffset=20):
    """Makes a second bottom spine"""
    if ax is None:
        ax = plt.gca()
    second_bottom = mpl.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
    second_bottom.set_position(('outward', offset))
    ax.spines['second_bottom'] = second_bottom

    if label is not None:
        # Make a new xlabel
        ax.annotate(label, 
                xy=(0.5, 0), xycoords='axes fraction', 
                xytext=(0, -labeloffset), textcoords='offset points', 
                verticalalignment='top', horizontalalignment='center')


# this function plots the damage sites wrt to atom number in the structure,
# with the atom numbers split up into chain type for visual aid
def damage_wrtatmnum(where,PDBmulti,densmet,epsilon):
	# epsilon sets a threshold such that density values (1+epsilon) and (1-epsilon) 
	# mean of the density values (using the densmet density metric). Note that 
	# epsilon set to 0 means no filtering is applied

	# deterine the start and stop positions of each chain 
	PDBmulti.sort(key=lambda x: x.atomnum)
	chains = []
	startend_atmnums = []
	for atom in PDBmulti:
		if atom.chaintype not in chains:
			chains.append(atom.chaintype)
			chainstart_atmnum = atom.atomnum
			prev_chainend = atom.atomnum-1
			startend_atmnums.append(atom.atomnum-1)
			startend_atmnums.append(atom.atomnum)


	# also add last atom number (end of last chain) and remove first element
	# since this does not correspond to a real atom number (-1..)
	startend_atmnums.append(PDBmulti[len(PDBmulti)-1].atomnum)
	startend_atmnums.pop(0)

	groups = []
	for i in range(0,len(chains)):
		groups.append((str(chains[i]), (startend_atmnums[2*i], startend_atmnums[2*i+1]) ) )

	# determine the list of density metric for each atom
	x = range(0,len(PDBmulti))
	densities = []
	for i in range(0,len(PDBmulti[0].meandensity)):
		if densmet == 'mean':
			densities.append([atom.meandensity[i] for atom in PDBmulti])
		elif densmet == 'median':
			densities.append([atom.mediandensity[i] for atom in PDBmulti])
		elif densmet == 'min':
			densities.append([atom.mindensity[i] for atom in PDBmulti])
		elif densmet == 'max':
			densities.append([atom.maxdensity[i] for atom in PDBmulti])
		else:
			print 'Unrecognised density metric'
			print '---> terminating script...'
			sys.exit()

	# apply determine the mean of the density metric and use this to set a 
	# epsilon filtering as defined above
	y_multi = []
	for i in range(0,len(PDBmulti[0].meandensity)):
		y_single = []
		mean_of_densmetric = np.mean(densities[i])
		for dens in densities[i]:
			if math.fabs(dens - mean_of_densmetric) >= float(epsilon):
				y_single.append(dens)
			else:
				y_single.append(mean_of_densmetric)
		y_multi.append(y_single)

    #-- Plot the results ------------------------------------------------------
	fig = plt.figure(figsize=(30, 20))
	ax = fig.add_subplot(111)

	# Give ourselves a bit more room at the bottom
	plt.subplots_adjust(bottom=0.2)
	#colorlist = ['#7479ce','#c75e5e','#34ceab','#1388a1','#ea85e5','#a17b4d','#a14d6a']
	colorlist = ['#ad1af0','#d11141','#00b159','#00aedb','#f37735','#f01ad1','#ffc425']
	for i in range(0,len(PDBmulti[0].meandensity)):
		ax.plot(x,y_multi[i],color=colorlist[i],linewidth=0.7,label = 'dataset '+str(i+1))

	# Drop the bottom spine by 40 pts
	ax.spines['bottom'].set_position(('outward', 40))

	# Make a second bottom spine in the position of the original bottom spine
	make_second_bottom_spine(label='Chain')

	# Annotate the groups
	for name, xspan in groups:
	    annotate_group(name, xspan)

	plt.xlabel('Atom number',fontsize=18)
	plt.ylabel(str(densmet)+ ' density loss',fontsize=18)
	plt.title('3clc: Atom number vs '+str(densmet)+' density loss',fontsize=20)
	plt.legend(loc='best',fontsize=18)

	#my_yticks = ['No damage','Damage']
	#plt.yticks([0,1], my_yticks)

	plt.savefig(where + str(densmet)+'densloss_wrtatmnum_3clc_epsilon'+str(epsilon)+'.png', bbox_inches="tight")



