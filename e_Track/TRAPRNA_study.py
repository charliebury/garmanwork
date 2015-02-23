# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def TRAPRNA_heterostudy(PDBmulti,datasetnum,atomtypes,resitypes):
	# function to plot violin plots for the TRAP-RNA case study
	# to determine whether there is any heterogeneity between 
	# the spread of protein damage between the RNA-bound and non-RNA
	# protein TRAP macromolecules

	# specify which chains are TRAP RNA bound, not RNA bound 
	# and also RNA
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	RNA = ['W']

	# determine the atom types, residue types to be included
	PDBmulti_included = []
	for atom in PDBmulti:
		if atomtypes != 'all' and resitypes != 'all':
			if (atom.atomtype in atomtypes and atom.basetype in resitypes):
				PDBmulti_included.append(atom)
		elif atomtypes == 'all' and resitypes != 'all':
			if atom.basetype in resitypes:
				PDBmulti_included.append(atom)
		elif atomtypes != 'all' and resitypes == 'all':
			if atom.atomtype in atomtypes:
				PDBmulti_included.append(atom)
		else:
			PDBmulti_included.append(atom)

	# loop to only include atoms in RNA bound TRAP closer than
	# 4.5 Angstroms to RNA. 
	distancelist = []
	PDBmulti_included2 = []
	for atom in PDBmulti_included:
		if atom.chaintype in RNAbound:
			for otheratom in PDBmulti:
				if otheratom.chaintype in RNA:
					distance = np.sqrt(np.square(atom.X_coord - otheratom.X_coord) +
	                           		   np.square(atom.Y_coord - otheratom.Y_coord) + 
	                                   np.square(atom.Z_coord - otheratom.Z_coord))
					distancelist.append(distance)
					if distance < 4.5:
						PDBmulti_included2.append(atom)
						break
	# additional loop to include respective atoms in RNA unbound TRAP
	# that would be closer than 4.5 Angstroms to RNA (if bound)
	for atom in PDBmulti_included:
		if atom.chaintype in nonRNAbound:
			for otheratom in PDBmulti_included2:
				if (atom.atomtype == otheratom.atomtype and
				   atom.residuenum == otheratom.residuenum and
				   atom.basetype == otheratom.basetype):
					PDBmulti_included2.append(atom)
					break


	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (8, 4)})
	fig = plt.figure()

	# determine the dataset number i to be chosen:
	i = datasetnum


	# determine mean density for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.meandensity[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.meandensity[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 1)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="muted", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label             
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('Mean density change', fontsize=12)



	# determine min density change for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.mindensity[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.mindensity[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 2)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="deep", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label            
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('min density change', fontsize=12)



	# determine max density change for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.maxdensity[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.maxdensity[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 3)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="pastel", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label            
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('max density change', fontsize=12)



	# determine Bfactor for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.Bfactor[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.Bfactor[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 4)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="bright", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label            
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('Bfactor', fontsize=12)


	fig.suptitle('TRAPRNA binding/non-binding statistics',
	             fontsize=20)   
	# Save the figure
	fig.savefig('TRAPRNA_meandensVSproteintype.png',bbox_inches='tight')

	return distancelist


def densitychange_v_atomnum(where,PDBmulti):
	# function to plot mean density change and min density change 
	# as function of atom number on same plot

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(2, 1, figsize=(100, 20), sharex=True)

	PDBmulti.sort(key=lambda x: x.atomnum)

	y = [atom.meandensity for atom in PDBmulti]
	x = [atom.atomnum for atom in PDBmulti]

	# add line plot for this dataset 
	axes[0].plot(x, y)
	axes[0].set_ylabel('mean change')

	y = [atom.mindensity for atom in PDBmulti]

	# add line plot for this dataset 
	axes[1].plot(x, y)
	axes[1].set_ylabel('min change')


	plt.xlabel('atom number', fontsize=18)
	f.suptitle('density change vs atom number',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	f.savefig(where+'densitychange_vs_atom_number.png')



def TRAPringcompare_format(PDBmulti):
	# function to compare the electron density loss for each atom in each TRAP
	# ring, to determine whether there is a difference between the atoms on average
	# between the two.

	# specify which chains are TRAP RNA bound, not RNA bound 
	# and also RNA
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	RNA = ['W']

	# determine the corresponding atom in each chain, starting with chain A
	# only considering protein TRAP chains here
	atomsbychain = []
	for atom in PDBmulti:
		if atom.chaintype in nonRNAbound[0]:
			atomperchain = []
			atomperchain.append(atom)
			for otheratom in PDBmulti:
				if otheratom.chaintype not in nonRNAbound[0]:
					if (atom.atomtype == otheratom.atomtype and
					   atom.residuenum == otheratom.residuenum and
					   atom.basetype == otheratom.basetype):
						atomperchain.append(otheratom)
						print otheratom.atomtype
			atomsbychain.append(atomperchain)

	# only want to include atoms which are present in all refined chains
	# since it was observed that some atoms only appeared in particular 
	# chains in the refined model
	atomsbychain_present = []
	for j in range(0,len(atomsbychain)):
		if len(atomsbychain[j]) != (len(nonRNAbound) + len(RNAbound)):
			print 'found!'
		else:
			atomsbychain_present.append(atomsbychain[j])

	return atomsbychain_present


def TRAPringcompare_plot(atomsbychain_present,atomtype,resitype,resinum,densmet):
	# following on from formatting function above, plot min density for corresponding 
	# atom found in all chains against chain number, with a separate subplot for each
	# dataset number

	# determine the desired atom by atom type, residue type and residue number
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	x = range(1,len(atomsbychain_present[0])+1)

	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		if densmet == 'min':
			y = [atom.mindensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'mean':
			y = [atom.meandensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'max':
			y = [atom.maxdensity[d_num-1] for atom in atomsbychain_present[counter]]
		else:
			print 'unrecognised density metric\n---> terminating script'
			sys.exit()

		ax.plot(x, y,linewidth=3)
		plt.axvline(11.5, color='#d64d4d', linestyle='dashed', linewidth=3)
		ax.text(0.05, 0.90, '# '+str(d_num), color='black',transform=ax.transAxes,fontsize=18)
	
	f.text(0.5, 0.04, 'chain number', ha='center',fontsize=18)
	f.text(0.04, 0.5, str(densmet)+' density loss', va='center', rotation='vertical',fontsize=18)

	f.suptitle(str(densmet)+' density change vs chain number: '+\
		       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)

	f.savefig('%s_densitychange_vs_chain_number_%s_%s_%s.png' %(str(densmet),str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))












