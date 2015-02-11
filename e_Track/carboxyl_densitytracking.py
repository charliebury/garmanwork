# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def carboxyl_lineplot(PDBmulti,where,densmet):
 	# function to plot a line plot to show the dynamics of each carboxyl group atom
 	# with repsect to dataset number

	# find atoms which are part of carboxyl groups
	carboxyl_atms = []
	for atom in PDBmulti:
		if ((atom.basetype in ('ASP') and atom.atomtype in ('CG','OD1','OD2')) or
			atom.basetype in ('GLU') and atom.atomtype in ('CD','OE1','OE2')):
			carboxyl_atms.append(atom)

	fig = plt.figure()
	ax = plt.axes()

	# Make some line plots
	x = range(0,len(PDBmulti[0].meandensity))

	CG_counter,OD1_counter,OD2_counter,CD_counter,OE1_counter,OE2_counter = 0,0,0,0,0,0

	# no show lines for you legend
	plt.plot([], label="ASP: CG", color="#b2ffdc")  
	plt.plot([], label="ASP: OD1", color="#bbaaee")  
	plt.plot([], label="ASP: OD2", color="#ffac93")  
	plt.plot([], label="ASP: CD", color="#d6f285")  
	plt.plot([], label="ASP: OE1", color="#ff80ff")  
	plt.plot([], label="ASP: OE2", color="#970898")  

	for atom in carboxyl_atms:	 
		if atom.atomtype == 'CG': 
			cc = '#b2ffdc'
		elif atom.atomtype == 'OD1':
			cc = '#bbaaee'
		elif atom.atomtype == 'OD2':
			cc = '#ffac93'
		elif atom.atomtype == 'CD':
			cc = '#d6f285'
		elif atom.atomtype == 'OE1':
			cc = '#ff80ff'
		elif atom.atomtype == 'OE2':
			cc = '#970898'

		# Make some line plots using the specified density metric
		if densmet == 'mean':
			ax.plot(x,atom.meandensity,color=cc,linewidth = 3)
		elif densmet == 'median':
			ax.plot(x,atom.mediandensity,color=cc,linewidth = 3)
		elif densmet == 'min':
			ax.plot(x,atom.mindensity,color=cc,linewidth = 3)
		elif densmet == 'max':
			ax.plot(x,atom.maxdensity,color=cc,linewidth = 3)

		plt.legend(loc='best')

		plt.title(str(densmet)+' density loss for carboxyl groups')
		plt.xlabel('Dataset')
		plt.ylabel(str(densmet)+' density loss')
		fig.savefig(str(where)+str(densmet)+'_3clc_carboxyls.png')













