# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 22:23:09 2015

@author: charlie
"""
from convert_maps2pkldobjectlists import maps2pkldobjs
from savevariables import retrieve_objectlist
from PDBfile_manipulation import multiARRAY_diffatomnumbers
from topdamagehits import topNdamsites_resibarplotter,topNdamsites_chainbarplotter,topNdamsites_printer
import os
from PDBmulti2txt import objlist2txt
from find_metricchange import find_Bchange
from PDBfile_manipulation import PDBtoCLASSARRAY_v2 as pdb2list
from bdam import numsurroundatoms_calculate,bdamage_calculate,numsurroundatms_extract
from bdamchange_vs_atomnum import bdamBfac_change_v_atomnum,bdamchange_v_atomnum


########################################################################
#                      master file for eTrack                          #
########################################################################

def map_processing():

	# location of map and pdb files required
	where = '../TRAPRNAdamage_26feb2015/'
	########################################################################

	# dataset names -> prefix only needed here. eg: '3clcdamage1' for 
	# '3clcdamage.pdb','3clcdamage_atoms.map','3clcdamage_density.map'
	pdbname = ['TRAPRNAdamage2','TRAPRNAdamage3','TRAPRNAdamage4','TRAPRNAdamage5',
			   'TRAPRNAdamage6','TRAPRNAdamage7','TRAPRNAdamage8','TRAPRNAdamage9','TRAPRNAdamage10']

	# want to create the following directories
	where2 = where+'output/'
	if not os.path.exists(where2):
	    os.makedirs(where2)

	where3 = where2+'plots/'
	if not os.path.exists(where3):
	    os.makedirs(where3)

	pkl_filenames = []
	for dataset in pdbname:
		mapfilname1 = dataset+'_atoms.map'
		mapfilname2 = dataset+'_density.map'
		pkl_filename = maps2pkldobjs(where,dataset,mapfilname1,mapfilname2)
		pkl_filenames.append(pkl_filename)

	return pkl_filenames 



def post_processing(pkl_filenames,initialPDBfile):

	print '••••••••••••••••••••••••••••••'
	print 'Reading in initial pdb file...'
	# next read in the pdb structure file:
	# run function to fill PDBarray list with atom objects from structure

	initialPDB = pdb2list(initialPDBfile,[])

	# determine the number of surrounding atoms for each atom in structure
	numsurroundatoms_calculate(initialPDBfile,initialPDB,10)

	# determine Bdamage metric for initial PDB structure
	bdamage_calculate(initialPDB)

	# retrieve object lists of atoms for each damage set
	print '•••••••••••••••••••••••••••••••'
	print 'Reading in damaged pkl files...'
	data_list = []
	for pkl_filename in pkl_filenames:
		print '\nDamage file number %s:' %(len(data_list)+1)
		PDB_ret = retrieve_objectlist(pkl_filename)

		# extract number of surrounding atoms for each atom in later structure
		# from initial structure
		numsurroundatms_extract(initialPDB,PDB_ret)

		# determine Bdamage metric for later PDB structures
		bdamage_calculate(PDB_ret)

		# add new retrieved damage set list to data_list
		data_list.append(PDB_ret)

	# create a list of atom objects with attributes as lists varying over 
	# dose range, only including atoms present in ALL damage datasets
	PDBmulti = multiARRAY_diffatomnumbers(data_list)

	# determine Bfactor change between later datasets and initial 
	find_Bchange(initialPDB,PDBmulti,'Bfactor')

	# determine Bdamage change between later datasets and initial 
	find_Bchange(initialPDB,PDBmulti,'Bdamage')

	print len(initialPDB)

	# write atom numbers and density metrics to a simple text file
	where = '../TRAPRNAdamage_26feb2015/output/'
	    
	objlist2txt(PDBmulti,where)

	# return the list of merged objects PDBmulti
	return PDBmulti



def graphanalysis(PDBmulti,N,densmet):

	# location of output plots (make folder if doesn't exist)
	where = '../TRAPRNAdamage_26feb2015/output/plots/combineddatasets/'
	if not os.path.exists(where):
	    os.makedirs(where)

	# determine top N damage sites:
	topNdamsites_resibarplotter(PDBmulti,N,where,densmet,'NOTnormalised')
	topNdamsites_resibarplotter(PDBmulti,N,where,densmet,'normalised')
	topNdamsites_chainbarplotter(PDBmulti,N,where,densmet)
	topNdamsites_printer(PDBmulti,N,where,densmet)

	# plot graphs of Bfactor change and Bdamage change against atom number
	bdamchange_v_atomnum(where,PDBmulti)
	bdamBfac_change_v_atomnum(where,PDBmulti)




# a = ['13830_TRAPRNAdamage2_data.pkl','13830_TRAPRNAdamage3_data.pkl',
# 	 '13830_TRAPRNAdamage4_data.pkl','13830_TRAPRNAdamage5_data.pkl',
# 	 '13830_TRAPRNAdamage6_data.pkl','13830_TRAPRNAdamage7_data.pkl',
# 	 '13830_TRAPRNAdamage8_data.pkl','13830_TRAPRNAdamage9_data.pkl',
# 	 '13830_TRAPRNAdamage10_data.pkl']






