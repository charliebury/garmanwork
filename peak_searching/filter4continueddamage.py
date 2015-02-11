# -*- coding: utf-8 -*-
"""
Created on Wed Oct 08 14:56:38 2014

@author: lina2532
"""

#8/10/2014: created new python script converting previous rotation1 matlab 
#script into python. Script filters distance.txt difference map peak files 
#further, by ensuring that if damage is observed in an early dataset, it
#must also be observed in a later dataset to be kept in the distance.txt file

# 2/02/2015: rewrote most of this script to read in pickled lists of objects
# (one for each damage data set) and process to create a new list of objects
# containing damage sites which once are present, remain present


#!/usr/bin/env python
import os
from math import exp, expm1, sqrt, pow
import datetime
import numpy as np
from lib.savevariables import retrieve_objectlist,save_objectlist
from lib.multidatasets_plots import barplot_atompropvsdatasets,barplot_damagesitesperdataset
from lib.glossy_output import *
import sys
from lib.write_mergedoutputfile import mergedoutput_writer


# A class for difference map peaks
class peakobj_multi:
    def __init__(self,peaknum=[],mindist=[],atomnum=[],residuenum=0,
                 atomtype="",basetype="",chaintype="",peakvalue=[],
                 X_coord=[],Y_coord=[],Z_coord=[],Bfactor=[],Occupancy=[]):

        self.peaknum = peaknum
        self.atomnum = atomnum
        self.residuenum = residuenum
        self.atomtype = atomtype
        self.basetype = basetype 
        self.chaintype = chaintype
        self.X_coord = X_coord
        self.Y_coord = Y_coord
        self.Z_coord = Z_coord
        self.Bfactor = Bfactor 
        self.Occupancy = Occupancy 
        self.mindist = mindist
        self.peakvalue = peakvalue   
    
# new title page for new script:
title_page2()
    
# MAKE A NEW FOLDER TO SAVE NEW FILES TO
print '\nWhat do you want to call your folder (date will be automatically included)?'
foldername = raw_input("Please enter something:\n")    
    
date =  datetime.datetime.now().strftime("%Y%m%d_%H%M") 
    
path = 'outputfiles_mergeddatasets/'+str(foldername)+str('_')+str(date);
if not os.path.exists(path):
    os.makedirs(path)

# determine how many damage sets are present
print 'How many damage sets are present'
num_damsets = raw_input('Enter a number please...\n')

filler_section()
print 'Reading inputfile4merging.txt to locate dataset names...'
# read in input file and convert pkl files for each dataset into multi-dim
# list PDBmulti, containing list of damage site objects for each dataset
inputfile = open('inputfile4merging.txt','r')
inputfilelines = inputfile.readlines()
PDBmulti = []
for line in inputfilelines:
    if 'DATA' in line.split()[0]:
        print line
        PDBretrieved = retrieve_objectlist(line.split()[1])
        PDBmulti.append(PDBretrieved)

# check that same number of damage sets read in as specified by user:
if len(PDBmulti) != int(num_damsets):
    print 'Inconsistency in number of datasets found in inputfile4merging.txt'
    print '---> not %s' %(str(num_damsets))
    sys.exit()

# want to do a check for atoms which have more than 1 diff peak associated with them
filler_section()
print 'For atoms in structure with multiple damage sites located per dataset, locating'+\
      'closest and removing all others...'
datacount = 0
for dataset in PDBmulti:
    datacount += 1
    counter = -1
    for atom in dataset:
        counter += 1
        for otheratom in dataset:
            if (atom.atomtype == otheratom.atomtype and atom.residuenum == otheratom.residuenum and
                atom.basetype == otheratom.basetype and atom.chaintype == otheratom.chaintype):
                if float(atom.mindist) > float(otheratom.mindist):
                    print 'Dataset %s: found atom with multiple diff peaks associated to it --> removing all but closest' %(str(datacount))
                    print str(atom.atomnum) + ' ' + str(atom.atomtype) + ' ' +\
                          str(atom.residuenum) + ' ' + str(atom.basetype) + ' ' + str(atom.chaintype)
                    dataset.pop(counter)
                    break
print '---> success!'

# next read through the lists of objects and create a new list of same form, only containing atoms which are present
# in all later lists
filler_section()
print 'Filtering to remove damage sites that once appear, are not found in ALL later datasets...'
PDBmulti_filtered = []

for i in range(0,len(PDBmulti)):
    PDB_filtered = []
    for atom in PDBmulti[i]:
        counter = 0
        for laterlist in PDBmulti[i:len(PDBmulti)]: 
            for otheratom in laterlist:
                if (atom.atomtype == otheratom.atomtype and atom.residuenum == otheratom.residuenum and
                    atom.basetype == otheratom.basetype and atom.chaintype == otheratom.chaintype):
                    counter += 1
        if counter == (len(PDBmulti) - i):
            PDB_filtered.append(atom)
        else:
            print 'Damage @ %s %s %s %s %s : dataset %s not found in all later datasets --> removing'\
                   %(str(atom.atomnum),str(atom.atomtype),str(atom.residuenum),
                     str(atom.basetype),str(atom.chaintype),str(i+1))
    PDBmulti_filtered.append(PDB_filtered)

damage_sites = []
for atom in PDBmulti_filtered[len(PDBmulti_filtered)-1]:

    peaknum_vec = []
    mindist_vec = []
    atomnum_vec = []
    peakvalue_vec = []
    X_coord_vec = []
    Y_coord_vec = []
    Z_coord_vec = []
    Bfactor_vec = []
    Occupancy_vec = []

    for i in range(0,len(PDBmulti_filtered)):
        counter = 0
        for otheratom in PDBmulti_filtered[i]:

            if (atom.atomtype == otheratom.atomtype and atom.residuenum == otheratom.residuenum and
                atom.basetype == otheratom.basetype and atom.chaintype == otheratom.chaintype):

                peaknum_vec.append(otheratom.peaknum)
                mindist_vec.append(otheratom.mindist)
                atomnum_vec.append(otheratom.atomnum)
                peakvalue_vec.append(otheratom.peakvalue)
                X_coord_vec.append(otheratom.X_coord)
                Y_coord_vec.append(otheratom.Y_coord)
                Z_coord_vec.append(otheratom.Z_coord)
                Bfactor_vec.append(otheratom.Bfactor)
                Occupancy_vec.append(otheratom.Occupancy)
                counter += 1 
                break

        if counter == 0:
            peaknum_vec.append(float('nan'))
            mindist_vec.append(float('nan'))
            atomnum_vec.append(float('nan'))
            peakvalue_vec.append(float('nan'))
            X_coord_vec.append(float('nan'))
            Y_coord_vec.append(float('nan'))
            Z_coord_vec.append(float('nan'))
            Bfactor_vec.append(float('nan'))
            Occupancy_vec.append(float('nan'))

    y = peakobj_multi(peaknum_vec,mindist_vec,atomnum_vec,atom.residuenum,atom.atomtype,atom.basetype,
                      atom.chaintype,peakvalue_vec,X_coord_vec,Y_coord_vec,Z_coord_vec,Bfactor_vec,Occupancy_vec)

    damage_sites.append(y)

# write an output txt file containing all the sites of damage
mergedoutput_writer('./'+path+'/',damage_sites,num_damsets)

# save the new list of objects using cPickle
save_objectlist(damage_sites,'mergedpeaksearch3clc')

# output some barplots over the datasets, to show different residue, atom and chain types
barplot_atompropvsdatasets('./'+path+'/',damage_sites,'residuetype')
barplot_atompropvsdatasets('./'+path+'/',damage_sites,'atomtype')
barplot_atompropvsdatasets('./'+path+'/',damage_sites,'chaintype')

barplot_damagesitesperdataset('./'+path+'/',damage_sites)
















    
    
    
    
    
    
    
    
    
    
    
    
    
    