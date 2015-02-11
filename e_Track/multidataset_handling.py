# -*- coding: utf-8 -*-
"""
Created on Sun Jan  4 01:53:56 2015
@author: charlie
"""
from masterscript import master

# This is the script where all the input data for multiple datasets is specified

# small class to hold info on each dataset:
# OPTIONS:
# Do you want to threshold electron density around each atom in structure 
# by VDW radius (type: 'vdw') or a global (non-element-specific) threshold 
# value chosen by the user (select value - see below) or a Bfactor weighted VDW radius 
# threshold per atom (type: 'bfac_vdw')?'

# If you have selected a global threshold level, then select your threshold
# level (in Angstroms)! Note that too low a threshold <1 may result in no 
# electron density being associated for a given atom in structure, resulting 
# in program breakdown. Consider grid spacing in electron density .map file 
# before choosing a suitable threshold). Threshold will be read as string and 
# converted to float later

class dataset_info:
    def __init__(self,mapname="",pdbname="",PDBDENSITY=[],thres_type="",voxel_list=[],disjoint_vxlperatom="",voxel_list_filtered=[]):
        self.mapname = mapname
        self.pdbname = pdbname
        self.PDBDENSITY = PDBDENSITY
        self.thres_type = thres_type
        self.voxel_list = voxel_list
        self.disjoint_vxlperatom = disjoint_vxlperatom
        self.voxel_list_filtered = voxel_list_filtered

    def get_pdbdens(self):
        self.PDBDENSITY = master(self)


############--------1---------##################
data2 = dataset_info('./map_files/3clcdiffmap6.map','./pdb_files/3clcdamage2',[],'vdw',[],'yes',[])
data2.get_pdbdens()

############--------2---------##################
data3 = dataset_info('./map_files/3clcdiffmap7.map','./pdb_files/3clcdamage3',[],'vdw',[],'yes',[])
data3.get_pdbdens()

############--------3---------##################
data4 = dataset_info('./map_files/3clcdiffmap8.map','./pdb_files/3clcdamage4',[],'vdw',[],'yes',[])
data4.get_pdbdens()

############--------4---------##################
data5 = dataset_info('./map_files/3clcdiffmap9.map','./pdb_files/3clcdamage5',[],'vdw',[],'yes',[])
data5.get_pdbdens()

############--------5---------##################
data6 = dataset_info('./map_files/3clcdiffmap10.map','./pdb_files/3clcdamage6',[],'vdw',[],'yes',[])
data6.get_pdbdens()

############--------6---------##################
data7 = dataset_info('./map_files/3clcdiffmap11.map','./pdb_files/3clcdamage7',[],'vdw',[],'yes',[])
data7.get_pdbdens()

############--------7---------##################
data8 = dataset_info('./map_files/3clcdiffmap12.map','./pdb_files/3clcdamage8',[],'vdw',[],'yes',[])
data8.get_pdbdens()





###############################################################################
from PDBfile_manipulation import PDBtoCLASSARRAY
PDBarray_empty = []
pdbfilename = './pdb_files/3clcdamage1'
# run function to fill PDBinitial list with atom objects from initial structure
PDBinitial = PDBtoCLASSARRAY(pdbfilename+'.pdb',PDBarray_empty)
# calculate Bdamage style metric for current pdb structure
from bdam import bdam_calculate
threshold = 14
bdam_calculate(PDBinitial,threshold) 
###############################################################################

from bdam import bdam_change,Bfactor_change
# calculate Bfactor change between initial data set and damaged datasets
data_list = [data2,data3,data4,data5,data6,data7,data8]
for data in data_list:
    bdam_change(PDBinitial,data.PDBDENSITY)
    Bfactor_change(PDBinitial,data.PDBDENSITY)

###############################################################################
# create a list of atom objects with attributes as lists varying over dose range
from PDBfile_manipulation import multiARRAY
PDBmulti = multiARRAY(data_list)