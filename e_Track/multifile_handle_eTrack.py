# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 22:23:09 2015

@author: charlie
"""
from CHECK_readinatommap import maps2atmdensity
from savevariables import save_objectlist, retrieve_objectlist

############--------1---------##################
# input file info here
where = './3clc_files/'
pdbname = '3clcdamage2'
mapfilname1 = '3clcdamage2_atoms.map'
mapfilname2 = '3clcdamage2_density.map'
PDB_dam2 = maps2atmdensity(where + pdbname,where + mapfilname1,'atom_map',where + mapfilname2,'density_map')

# pickle the big list of atom objects
pkl_filename2 = save_objectlist(PDB_dam2,pdbname)
PDB_dam2 = []

############--------2---------##################
# input file info here
pdbname = '3clcdamage3'
mapfilname1 = '3clcdamage3_atoms.map'
mapfilname2 = '3clcdamage3_density.map'
PDB_dam3 = maps2atmdensity(where + pdbname,where + mapfilname1,'atom_map',where + mapfilname2,'density_map')

# pickle the big list of atom objects
pkl_filename3 = save_objectlist(PDB_dam3,pdbname)
PDB_dam3 = []

############--------3---------##################
# input file info here
pdbname = '3clcdamage4'
mapfilname1 = '3clcdamage4_atoms.map'
mapfilname2 = '3clcdamage4_density.map'
PDB_dam4 = maps2atmdensity(where + pdbname,where + mapfilname1,'atom_map',where + mapfilname2,'density_map')

# pickle the big list of atom objects
pkl_filename4 = save_objectlist(PDB_dam4,pdbname)
PDB_dam4 = []

############--------4---------##################
# input file info here
pdbname = '3clcdamage5'
mapfilname1 = '3clcdamage5_atoms.map'
mapfilname2 = '3clcdamage5_density.map'
PDB_dam5 = maps2atmdensity(where + pdbname,where + mapfilname1,'atom_map',where + mapfilname2,'density_map')

# pickle the big list of atom objects
pkl_filename5 = save_objectlist(PDB_dam5,pdbname)
PDB_dam5 = []

############--------5---------##################
# input file info here
pdbname = '3clcdamage6'
mapfilname1 = '3clcdamage6_atoms.map'
mapfilname2 = '3clcdamage6_density.map'
PDB_dam6 = maps2atmdensity(where + pdbname,where + mapfilname1,'atom_map',where + mapfilname2,'density_map')

# pickle the big list of atom objects
pkl_filename6 = save_objectlist(PDB_dam6,pdbname)
PDB_dam6 = []

############--------6---------##################
# input file info here
pdbname = '3clcdamage7'
mapfilname1 = '3clcdamage7_atoms.map'
mapfilname2 = '3clcdamage7_density.map'
PDB_dam7 = maps2atmdensity(where + pdbname,where + mapfilname1,'atom_map',where + mapfilname2,'density_map')

# pickle the big list of atom objects
pkl_filename7 = save_objectlist(PDB_dam7,pdbname)
PDB_dam7 = []

############--------7---------##################
# input file info here
pdbname = '3clcdamage8'
mapfilname1 = '3clcdamage8_atoms.map'
mapfilname2 = '3clcdamage8_density.map'
PDB_dam8 = maps2atmdensity(where + pdbname,where + mapfilname1,'atom_map',where + mapfilname2,'density_map')

# pickle the big list of atom objects
pkl_filename8 = save_objectlist(PDB_dam8,pdbname)
PDB_dam8 = []


###############################################################################
from PDBfile_manipulation import PDBtoCLASSARRAY
PDBarray_empty = []
pdbfilename = '3clcdamage1'
# run function to fill PDBinitial list with atom objects from initial structure
PDBinitial = PDBtoCLASSARRAY(where + pdbfilename+'.pdb',PDBarray_empty)
###############################################################################

pkl_filename2 = '3927_3clcdamage2_data.pkl'
pkl_filename3 = '3927_3clcdamage3_data.pkl'
pkl_filename4 = '3927_3clcdamage4_data.pkl'
pkl_filename5 = '3927_3clcdamage5_data.pkl'
pkl_filename6 = '3927_3clcdamage6_data.pkl'
pkl_filename7 = '3927_3clcdamage7_data.pkl'
pkl_filename8 = '3927_3clcdamage8_data.pkl'

# retrieve object lists of atoms for each damage set
PDB_dam2_ret = retrieve_objectlist(pkl_filename2)
PDB_dam3_ret = retrieve_objectlist(pkl_filename3)
PDB_dam4_ret = retrieve_objectlist(pkl_filename4)
PDB_dam5_ret = retrieve_objectlist(pkl_filename5)
PDB_dam6_ret = retrieve_objectlist(pkl_filename6)
PDB_dam7_ret = retrieve_objectlist(pkl_filename7)
PDB_dam8_ret = retrieve_objectlist(pkl_filename8)


#from bdam import Bfactor_change
## calculate Bfactor change between initial data set and damaged datasets
data_list = [PDB_dam2_ret,PDB_dam3_ret,PDB_dam4_ret,PDB_dam5_ret,PDB_dam6_ret,PDB_dam7_ret,PDB_dam8_ret]
#for dataset in data_list:
#    Bfactor_change(PDBinitial,dataset)

###############################################################################
# create a list of atom objects with attributes as lists varying over dose range
from PDBfile_manipulation import multiARRAY_PDBonly
PDBmulti = multiARRAY_PDBonly(data_list)