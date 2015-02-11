# -*- coding: utf-8 -*-

from savevariables import retrieve_objectlist
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math


# retrieve the pickled list of objects
#PDBmulti = retrieve_objectlist('602_mergedpeaksearch3clc_data.pkl')

def barplot_resvsdatasets(path,PDBmulti):
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (12, 6)})
	ax = plt.subplots()


	PDBmulti.sort(key=lambda x: x.basetype)
	print '----------------------------------------'
	print 'Summary of damage sites by residue type:'

	unique_list_multi = []
	resi_freq_list_multi = []

	for i in range(0,len(PDBmulti[0].peaknum)):
		print '------------------------'
		print 'Dataset %s:' %(str(i+1))
		unique_list = []
		resi_freq_list = []
		resi_list = [atom.basetype for atom in PDBmulti if math.isnan(atom.peaknum[i]) != True]
		for resi in resi_list:
			if resi not in unique_list:
				unique_list.append(resi)
				resfreq  = resi_list.count(resi)
				resi_freq_list.append(resfreq)
				print '\t%s ---> %s' %(str(resi),str(resfreq))
			else:
				pass
		unique_list_multi.append(unique_list)
		resi_freq_list_multi.append(resi_freq_list)

	# the full number of distinct residue types by the last dataset is:
	full_res_present = unique_list_multi[len(PDBmulti[0].peaknum)-1]

	# next want to 'bulk out' the earlier dataset lists, to include the residue
	# types present in the later datasets but not in the earlier datasets
	for i in range(0,len(PDBmulti[0].peaknum)):
		for resi in full_res_present:
			if resi not in unique_list_multi[i]:
				 unique_list_multi[i].append(resi)
				 resi_freq_list_multi[i].append(0)
		unique_list_multi[i], resi_freq_list_multi[i] = (list(t) for t in zip(*sorted(zip(unique_list_multi[i], resi_freq_list_multi[i]))))


	N = len(unique_list_multi[0])
	ind = np.arange(N)  # the x locations for the groups
	width = 0.1       # the width of the bars

	fig = plt.figure()
	ax = fig.add_subplot(111)

	for i in range(0,len(PDBmulti[0].peaknum)):

		colorlist = ['#008744','#0057e7','#d62d20','#ffa700','#eeeeee']*len(PDBmulti[0].peaknum)

		yvals = resi_freq_list_multi[i]
		rect = ax.bar(ind+width*i, yvals, width, color=colorlist[i],label='dataset '+str(i+1))

	ax.set_ylabel('Frequency')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(tuple(unique_list_multi[0]))
	ax.set_xlabel('Residue/base type')
	ax.legend()
	fig.suptitle('Residue/base type vs frequency')

	def autolabel(rects):
	    for rect in rects:
	        h = rect.get_height()
	        ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
	                ha='center', va='bottom')

	plt.savefig(path+'barplot_resvsdatasets.png', bbox_inches="tight")




def barplot_atompropvsdatasets(path,PDBmulti,atomprop):
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (12, 6)})
	ax = plt.subplots()

	print '----------------------------------------'
	if atomprop == 'residuetype':		
		PDBmulti.sort(key=lambda x: x.basetype)
		print 'Summary of damage sites by residue type:'
	elif atomprop == 'atomtype':
		PDBmulti.sort(key=lambda x: x.atomtype)
		print 'Summary of damage sites by atom type type:'
	elif atomprop == 'chaintype':
		PDBmulti.sort(key=lambda x: x.chaintype)
		print 'Summary of damage sites by chain type:'

	unique_list_multi = []
	resi_freq_list_multi = []

	for i in range(0,len(PDBmulti[0].peaknum)):
		print '------------------------'
		print 'Dataset %s:' %(str(i+1))
		unique_list = []
		resi_freq_list = []

		if atomprop == 'residuetype':		
			resi_list = [atom.basetype for atom in PDBmulti if math.isnan(atom.peaknum[i]) != True]
		if atomprop == 'atomtype':		
			resi_list = [atom.atomtype for atom in PDBmulti if math.isnan(atom.peaknum[i]) != True]
		if atomprop == 'chaintype':		
			resi_list = [atom.chaintype for atom in PDBmulti if math.isnan(atom.peaknum[i]) != True]

		for resi in resi_list:
			if resi not in unique_list:
				unique_list.append(resi)
				resfreq  = resi_list.count(resi)
				resi_freq_list.append(resfreq)
				print '\t%s ---> %s' %(str(resi),str(resfreq))
			else:
				pass
		unique_list_multi.append(unique_list)
		resi_freq_list_multi.append(resi_freq_list)

	# the full number of distinct residue types by the last dataset is:
	full_res_present = unique_list_multi[len(PDBmulti[0].peaknum)-1]

	# next want to 'bulk out' the earlier dataset lists, to include the residue
	# types present in the later datasets but not in the earlier datasets
	for i in range(0,len(PDBmulti[0].peaknum)):
		for resi in full_res_present:
			if resi not in unique_list_multi[i]:
				 unique_list_multi[i].append(resi)
				 resi_freq_list_multi[i].append(0)
		unique_list_multi[i], resi_freq_list_multi[i] = (list(t) for t in zip(*sorted(zip(unique_list_multi[i], resi_freq_list_multi[i]))))


	N = len(unique_list_multi[0])
	ind = np.arange(N)  # the x locations for the groups
	width = 0.1       # the width of the bars

	fig = plt.figure()
	ax = fig.add_subplot(111)

	for i in range(0,len(PDBmulti[0].peaknum)):

		colorlist = ['#008744','#0057e7','#d62d20','#ffa700','#eeeeee']*len(PDBmulti[0].peaknum)

		yvals = resi_freq_list_multi[i]
		rect = ax.bar(ind+width*i, yvals, width, color=colorlist[i],label='dataset '+str(i+1))

	ax.set_ylabel('Frequency')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(tuple(unique_list_multi[0]))


	ax.set_xlabel(str(atomprop))
	ax.legend()
	fig.suptitle(str(atomprop)+ ' vs frequency')

	def autolabel(rects):
	    for rect in rects:
	        h = rect.get_height()
	        ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
	                ha='center', va='bottom')

	plt.savefig(path+'barplot_multidatasets_'+str(atomprop)+'.png', bbox_inches="tight")




def barplot_damagesitesperdataset(path,PDBmulti):
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (12, 6)})
	ax = plt.subplots()

	fig = plt.figure()
	ax = fig.add_subplot(111)

	# determine the number of damage sites located per dataset
	num_damsites = []
	print '--------------------------------------------------'
	print 'Final filtered number of damage sites per dataset:'
	for i in range(0,len(PDBmulti[0].peaknum)):
		counter = 0
		for damsite in PDBmulti:
			if math.isnan(damsite.peaknum[i]) == False:
				counter += 1
		num_damsites.append(counter)
		
		print '\tDataset %s ---> %s' %(str(i+1),str(counter))


	# detemine a colour scheme and x-axis labels for the barplot
	xaxislabels = []
	for i in range(0,len(PDBmulti[0].peaknum)):
		colorlist = ['#008744','#0057e7','#d62d20','#ffa700','#eeeeee']*len(PDBmulti[0].peaknum)
		xaxislabels.append(str(i+1))

	ind = np.arange(len(num_damsites))

	rect = ax.bar(ind, num_damsites, color=colorlist)

	ax.set_ylabel('Frequency')
	ax.set_xticklabels(np.array(xaxislabels))


	ax.set_xlabel('Dataset number')
	ax.legend()
	fig.suptitle('Dataset vs Damage Frequency')

	def autolabel(rects):
	    for rect in rects:
	        h = rect.get_height()
	        ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
	                ha='center', va='bottom')

	plt.savefig(path+'barplot_damageperdataset.png', bbox_inches="tight")



