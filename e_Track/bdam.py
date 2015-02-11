# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 17:08:27 2014

@author: charlie
"""
import numpy as np
from classholder_v2 import StructurePDB
import sys

##function determines for each atom in structure the number of neighbouring atoms within 
#a threshold (defined above) for all atoms. For each atom, number of contacts added
#as class attribute for atom
def bdam_calculate(PDBarray,threshold):
    
    ####    #####   ####    ####
    #unessential loading bar add-in
    total = len(PDBarray)
    point = total / 100
    increment = total / 20
    ####    #####   ####    ####

    print 'Calculating contact number for atoms in structure...'
    counter_bar = 0
    for atom in PDBarray:
        counter_bar += 1
        
        ####    #####   ####    ####
        #unessential loading bar add-in
        if(counter_bar % (5 * point) == 0):
                sys.stdout.write("\r[" + "=" * (counter_bar / increment) +  " " * ((total - counter_bar)/ increment) + "]" +  str(counter_bar / point) + "%")
                sys.stdout.flush()        
        ####    ####    ####    ####        
        
        num_contacts = 0
        for otheratoms in PDBarray:
            distance = np.sqrt(np.square(atom.X_coord - otheratoms.X_coord) + np.square(atom.Y_coord - otheratoms.Y_coord) + np.square(atom.Z_coord - otheratoms.Z_coord))
            if distance != 0 and distance < threshold:
                num_contacts += 1
        atom.numsurroundatoms = num_contacts
    
    ##Next calculate Bdamage style metric for each atom, to save bdam attribute for each atom
    for atom in PDBarray:
        print 'Calculating bdam style metric for atoms in structure...'
        simpacking_bfactors = []
        for otheratoms in PDBarray:
            if round(atom.numsurroundatoms/10) == round(otheratoms.numsurroundatoms/10):
                simpacking_bfactors.append(float(otheratoms.Bfactor))
        bdam = float(atom.Bfactor)/(np.mean(simpacking_bfactors))
        
        atom.bdam = bdam




##function to calculate for damaged structure the change in bdamage from the initial structure
def bdam_change(initialpdb,laterpdb):
    #first ensure both PDB lists sorted by atom number
    initialpdb.sort(key=lambda x: x.atomnum)
    laterpdb.sort(key=lambda x: x.atomnum)
    
    for atom in laterpdb:
        atom.bdamchange = float(atom.bdam) -float(initialpdb[int(atom.atomnum)-1].bdam)
        
        
        
        
##function to calculate for damaged structure the change in Bfactor from the initial structure
def Bfactor_change(initialpdb,laterpdb):
    #first ensure both PDB lists sorted by atom number
    initialpdb.sort(key=lambda x: x.atomnum)
    laterpdb.sort(key=lambda x: x.atomnum)
    
    for atom in laterpdb:
        atom.Bfactorchange = float(atom.Bfactor) -float(initialpdb[int(atom.atomnum)-1].Bfactor)
        

