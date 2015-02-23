# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from operator import truediv

def metricVSdensity_scatter(where,PDBmulti,atmmetric,densmetric,resitypes,atomtypes,chaintypes):
	# this function plots as scatter plot of density metric (mean,
	# median, min, max) vs other atomic metric (Bfactor, Bdamage etc)
	# for each dataset individually, to determine correlations between
	# the two chosen metrics.
	# resitypes specifies which residues should be included --> if 'all' is
	# specified then all residues are included

	# determine the atoms types to include, specified by residue type, atom type, chain type
	includedatoms1 = []
	for atom in PDBmulti:
		if resitypes == 'all':
			includedatoms1.append(atom)
		else:
			if atom.basetype in resitypes:
				includedatoms1.append(atom)

	includedatoms2 = []
	for atom in includedatoms1:
		if atomtypes == 'all':
			includedatoms2.append(atom)
		else:
			if atom.atomtype in atomtypes:
				includedatoms2.append(atom)

	includedatoms3 = []
	for atom in includedatoms2:
		if chaintypes == 'all':
			includedatoms3.append(atom)
		else:
			if atom.chaintype in chaintypes:
				includedatoms3.append(atom)

	includedatoms_final = includedatoms3

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(1, 1, figsize=(12, 12),)

	# number of datasets present
	num_datasets = len(PDBmulti[0].meandensity)

	# number of atoms present in structure
	num_atoms = len(includedatoms_final)

	multi_densmetric = []
	multi_atmmetric = []
	multi_colorlist = []

	# create list for atomic metric for each dataset individually
	for i in range(0,num_datasets):
		atmmetric_vector = []
		if atmmetric == 'bfactor':
			for atom in includedatoms_final:
				atmmetric_vector.append(atom.Bfactor[i])
		elif atmmetric == 'bdamage':
			for atom in includedatoms_final:
				atmmetric_vector.append(atom.bdam[i])
		elif atmmetric == 'bdamchange':
			for atom in includedatoms_final:
				atmmetric_vector.append(atom.bdamchange[i])
		elif atmmetric == 'Bfactorchange':
			for atom in includedatoms_final:
				atmmetric_vector.append(atom.Bfactorchange[i])
		else:
			print 'Unrecognised atomic metric\n--->terminating script...'
			sys.exit()

		densmetric_vector = []
		if densmetric == 'max':
			for atom in includedatoms_final:
				densmetric_vector.append(atom.maxdensity[i])
		elif densmetric == 'min':
			for atom in includedatoms_final:
				densmetric_vector.append(atom.mindensity[i])
		elif densmetric == 'mean':
			for atom in includedatoms_final:
				densmetric_vector.append(atom.meandensity[i])
		elif densmetric == 'median':
			for atom in includedatoms_final:
				densmetric_vector.append(atom.mediandensity[i])
		else:
			print 'Unrecognised density metric\n--->terminating script...'
			sys.exit()

		multi_densmetric.append(densmetric_vector)
		multi_atmmetric.append(atmmetric_vector)
		multi_colorlist.append([i]*num_atoms)

	cm = plt.cm.get_cmap('hsv')
	sc = plt.scatter(sum(multi_densmetric,[]),sum(multi_atmmetric,[]), s=24, c=sum(multi_colorlist,[]), alpha=0.5,cmap=cm)
	plt.xlabel(str(densmetric)+' density', fontsize=18)
	plt.ylabel(str(atmmetric)+' density', fontsize=18)
	plt.colorbar(sc)
	f.suptitle(str(densmetric)+' density vs '+ str(atmmetric),fontsize=20)
	f.savefig(str(where) + str(densmetric)+'_density_vs_'+ str(atmmetric)+'.png')







