# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 22:23:09 2015

@author: charlie
"""
from convert_maps2pkldobjectlists import maps2pkldobjs

########################################################################
#             master file for Oli's protein datasets                   #
########################################################################

def map_processing():

	# location of map and pdb files required
	where = './Oli_proteinfiles/'
	########################################################################

	# the individual dataset inputs begin here
	############--------1---------##################
	# input file info here for xfel dataset
	pdbname = 'xfel'
	mapfilname1 = 'xfel_atoms.map'
	mapfilname2 = 'xfel_density.map'
	pkl_filename1 = maps2pkldobjs(where,pdbname,mapfilname1,mapfilname2)

	############--------2---------##################
	# input file info here for sync 4.1A dataset
	pdbname = 'sync41A'
	mapfilname1 = 'sync41A_atoms.map'
	mapfilname2 = 'sync41A_density.map'
	pkl_filename2 = maps2pkldobjs(where,pdbname,mapfilname1,mapfilname2)

	############--------3---------##################
	# input file info here for sync 3.6A dataset
	pdbname = 'sync36A'
	mapfilname1 = 'sync36A_atoms.map'
	mapfilname2 = 'sync36A_density.map'
	pkl_filename3 = maps2pkldobjs(where,pdbname,mapfilname1,mapfilname2)



def post_processing():

	pkl_filename1 = '10723_xfel_data.pkl'
	pkl_filename2 = '10788_sync41A_data.pkl'
	pkl_filename3 = '10784_sync36A_data.pkl'

	########################################################################
	# retrieve object lists of atoms for each damage set
	from savevariables import retrieve_objectlist

	PDB_xfel_ret = retrieve_objectlist(pkl_filename1)
	PDB_sync41A_ret = retrieve_objectlist(pkl_filename2)
	PDB_sync36A_ret = retrieve_objectlist(pkl_filename3)


	data_list = [PDB_xfel_ret,PDB_sync41A_ret,PDB_sync36A_ret]


	########################################################################
	# create a list of atom objects with attributes as lists varying over 
	# dose range
	from PDBfile_manipulation import multiARRAY_diffatomnumbers
	PDBmulti = multiARRAY_diffatomnumbers(data_list)