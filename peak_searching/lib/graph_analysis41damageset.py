# -*- coding: utf-8 -*-

# Graphical analysis of filtered difference map peaks for 1 damage set
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker
import matplotlib as mpl

def peakfiltering_bystage(datedir,num_peaks_kept_list):
	# this function reads in a list of the number of peaks at each stage
	# of the filtering process and outputs a line graph to show the reduction
	# of number of peaks at each stage
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (8, 4)})
	ax = plt.subplots()

	x = [1,2,3,4]
	y = []
	for num in num_peaks_kept_list:
		y.append(int(num))

	plt.plot(x,y,color='#7479ce',linewidth=3.0)

	my_xticks = ['unfiltered','neg only','distance\nfilter','manual\nsearch']
	plt.xticks(x, my_xticks)

	plt.xlabel('Filtering stage', fontsize=16, color='black')
	plt.ylabel('Peaks remaining', fontsize=16, color='black')
	plt.title('Number of peaks filtered')

	# plt.show()
	plt.savefig(datedir.path + '/peakfiltering_bystage.png', bbox_inches="tight")


def barplot_resfreq(datedir,peakobjlist):
	# this function reads in the list of diff peaks assigned to atoms in the PDB 
	# structure and displays a bar plot of the frequency of each residue type 

	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (8, 4)})
	ax = plt.subplots()


	peakobjlist.sort(key=lambda x: x.basetype)
	print '----------------------------------------'
	print 'Summary of damage sites by residue type:'

	resi_list = [atom.basetype for atom in peakobjlist]
	unique_list = []
	resi_freq_list = []

	for resi in resi_list:
		if resi not in unique_list:
			unique_list.append(resi)
			resfreq  = resi_list.count(resi)
			resi_freq_list.append(resfreq)
			print '\t%s ---> %s' %(str(resi),str(resfreq))
		else:
			pass

	sns.barplot(x=np.array(unique_list),y=np.array(resi_freq_list),palette="BuGn_d")


	plt.xlabel('Residue/base type', fontsize=16, color='black')
	plt.ylabel('Frequency', fontsize=16, color='black')
	plt.title('Frequency per residue/base type')

	sns.despine(bottom=True)

	# plt.show()
	plt.savefig(datedir.path + '/barplot_resfreq.png', bbox_inches="tight")



def barplot_chainfreq(datedir,peakobjlist):
	# this function reads in the list of diff peaks assigned to atoms in the PDB 
	# structure and displays a bar plot of the frequency of each residue type 

	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (8, 4)})
	ax = plt.subplots()


	peakobjlist.sort(key=lambda x: x.chaintype)
	print '----------------------------------------'
	print 'Summary of damage sites by chain type:'

	chain_list = [atom.chaintype for atom in peakobjlist]
	unique_list = []
	chain_freq_list = []

	for chain in chain_list:
		if chain not in unique_list:
			unique_list.append(chain)
			chainfreq  = chain_list.count(chain)
			chain_freq_list.append(chainfreq)
			print '\t%s ---> %s' %(str(chain),str(chainfreq))
		else:
			pass

	sns.barplot(x=np.array(unique_list),y=np.array(chain_freq_list),palette="Paired")


	plt.xlabel('Chain type', fontsize=16, color='black')
	plt.ylabel('Frequency', fontsize=16, color='black')
	plt.title('Frequency per chain type')

	sns.despine(bottom=True)

	# plt.show()
	plt.savefig(datedir.path + '/barplot_chainfreq.png', bbox_inches="tight")




def barplot_atompropfreq(datedir,peakobjlist,atomprop):
	# this function reads in the list of diff peaks assigned to atoms in the PDB 
	# structure and displays a bar plot of the frequency of each atom by atom property
	# 'atomprop' = ('basetype','atomtype','chaintype')
 
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (8, 4)})
	ax = plt.subplots()

	print '----------------------------------------'
	print 'Summary of damage sites by %s:' %(str(atomprop))

	if atomprop == 'basetype':
		peakobjlist.sort(key=lambda x: x.basetype)
		atmprop_list = [atom.basetype for atom in peakobjlist]
	elif atomprop == 'atomtype':
		peakobjlist.sort(key=lambda x: x.atomtype)
		atmprop_list = [atom.atomtype for atom in peakobjlist]
	elif atomprop == 'chaintype':
		peakobjlist.sort(key=lambda x: x.chaintype)
		atmprop_list = [atom.chaintype for atom in peakobjlist]

	unique_list = []
	atmprop_freq_list = []

	for atmprop in atmprop_list:
		if atmprop not in unique_list:
			unique_list.append(atmprop)
			atmpropfreq  = atmprop_list.count(atmprop)
			atmprop_freq_list.append(atmpropfreq)
			print '\t%s ---> %s' %(str(atmprop),str(atmpropfreq))
		else:
			pass

	sns.barplot(x=np.array(unique_list),y=np.array(atmprop_freq_list),palette="Paired")


	plt.xlabel(str(atomprop), fontsize=16, color='black')
	plt.ylabel('Frequency', fontsize=16, color='black')
	plt.title('Frequency per '+str(atomprop))

	sns.despine(bottom=True)

	# plt.show()
	plt.savefig(datedir.path + '/barplot_'+str(atomprop)+'_freq.png', bbox_inches="tight")






## the following two functions override the default behavior or twiny()
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)

def make_spine_invisible(ax, direction):
    if direction in ["right", "left"]:
        ax.yaxis.set_ticks_position(direction)
        ax.yaxis.set_label_position(direction)
    elif direction in ["top", "bottom"]:
        ax.xaxis.set_ticks_position(direction)
        ax.xaxis.set_label_position(direction)
    else:
        raise ValueError("Unknown Direction : %s" % (direction,))

    ax.spines[direction].set_visible(True)



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
def damage_wrtatmnum(datedir,peakobjlist,PDBarray):
	x = range(0,len(PDBarray))

	y = []
	# determine the list of atom numbers associated with damage sites
	damagesite_atmlist = [int(atom.atomnum) for atom in peakobjlist]

	for i in range(0,len(PDBarray)):
		if i in damagesite_atmlist:
			y.append(1)
		else:
			y.append(0)

	peakobjlist.sort(key=lambda x: x.atomnum)
	chains = []
	startend_atmnums = []
	for i in range(0,len(PDBarray)):
		if PDBarray[i][21] not in chains:
			chains.append(PDBarray[i][21])
			chainstart_atmnum = i
			prev_chainend = i-1
			startend_atmnums.append(i-1)
			startend_atmnums.append(i)


	# also add last atom number (end of last chain) and remove first element
	# since this does not correspond to a real atom number (-1..)
	startend_atmnums.append(i)
	startend_atmnums.pop(0)

	groups = []
	for i in range(0,len(chains)):
		groups.append((str(chains[i]), (startend_atmnums[2*i], startend_atmnums[2*i+1]) ) )


    #-- Plot the results ------------------------------------------------------
	fig = plt.figure()
	ax = fig.add_subplot(111)

	# Give ourselves a bit more room at the bottom
	plt.subplots_adjust(bottom=0.2)
	ax.plot(x,y,color='#7479ce',linewidth=0.5)

	# Drop the bottom spine by 40 pts
	ax.spines['bottom'].set_position(('outward', 40))

	# Make a second bottom spine in the position of the original bottom spine
	make_second_bottom_spine(label='Chain')

	# Annotate the groups
	for name, xspan in groups:
	    annotate_group(name, xspan)

	plt.xlabel('Atom number')
	plt.title('Atom number vs damage sites')

	my_yticks = ['No damage','Damage']
	plt.yticks([0,1], my_yticks)

	plt.savefig(datedir.path + '/damage_wrtatmnum.png', bbox_inches="tight")









