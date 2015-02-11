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
########################################################################
#             master file for Oli's protein datasets                   #
########################################################################

def map_processing():

	# location of map and pdb files required
	where = '../2BN3damage_11feb2014/'
	########################################################################

	# dataset names -> prefix only needed here. eg: '3clcdamage1' for 
	# '3clcdamage.pdb','3clcdamage_atoms.map','3clcdamage_density.map'
	pdbname = ['2BN3damage2','2BN3damage3','2BN3damage4','2BN3damage5',
			   '2BN3damage6','2BN3damage7','2BN3damage8','2BN3damage9','2BN3damage10']

	# want to create the following directories
	where2 = '../2BN3damage_11feb2014/output/'
	if not os.path.exists(where2):
	    os.makedirs(where2)

	where3 = '../2BN3damage_11feb2014/output/plots/'
	if not os.path.exists(where3):
	    os.makedirs(where3)

	pkl_filenames = []
	for dataset in pdbname:
		mapfilname1 = dataset+'_atoms.map'
		mapfilname2 = dataset+'_density.map'
		pkl_filename = maps2pkldobjs(where,dataset,mapfilname1,mapfilname2)
		pkl_filenames.append(pkl_filename)

	return pkl_filenames 



def post_processing(pkl_filenames):

	# retrieve object lists of atoms for each damage set
	data_list = []

	for pkl_filename in pkl_filenames:
		PDB_ret = retrieve_objectlist(pkl_filename)
		data_list.append(PDB_ret)

	# create a list of atom objects with attributes as lists varying over 
	# dose range, only including atoms present in ALL damage datasets
	PDBmulti = multiARRAY_diffatomnumbers(data_list)

	# write atom numbers and density metrics to a simple text file
	where = '../2BN3damage_11feb2014/output/'
	    
	objlist2txt(PDBmulti,where)

	# return the list of merged objects PDBmulti
	return PDBmulti


def graphanalysis(PDBmulti,N,densmet):

	# location of output plots (make folder if doesn't exist)
	where = '../2BN3damage_11feb2014/output/plots/combineddatasets/'
	if not os.path.exists(where):
	    os.makedirs(where)

	# determine top N damage sites:
	topNdamsites_resibarplotter(PDBmulti,N,where,densmet,'NOTnormalised')
	topNdamsites_resibarplotter(PDBmulti,N,where,densmet,'normalised')
	topNdamsites_chainbarplotter(PDBmulti,N,where,densmet)
	topNdamsites_printer(PDBmulti,N,where,densmet)








