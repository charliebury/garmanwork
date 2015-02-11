# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 02:39:27 2014

@author: charlie
"""
from StructurePDB import StructurePDB

def PDBtoCLASSARRAY(pdbfilename):
    PDBarray = []
    pdbin = open(str(pdbfilename), "r")
    lines = pdbin.readlines()
    for line in lines:
        if 'ATOM' in line.split()[0]:
            print line
            y = StructurePDB()
            y.atomnum = line.split()[1]
            y.atomtype = line.split()[2]
            y.basetype = line.split()[3]
            y.chaintype = line.split()[4]
            y.residuenum = line.split()[5]
            y.X_coord = line.split()[6]
            y.Y_coord = line.split()[7]
            y.Z_coord = line.split()[8]
            print y.atomnum
            PDBarray.append(y)
        else: 
            pass
    pdbin.close()
    return PDBarray


def PDBtoCLASSARRAY_fixedlinewidth(pdbfilename):
    #this variant of the above utilises the fixed line width of the PDB file
    PDBarray = []
    pdbin = open(str(pdbfilename), "r")
    lines = pdbin.readlines()
    for line in lines:
        if 'ATOM' in line[0:6]:
            print line
            y = StructurePDB()
            y.atomnum = line[6:11]
            y.atomtype = line[12:16]
            y.basetype = line[17:20]
            y.chaintype = line[21]
            y.residuenum = line[22:26]
            y.X_coord = line[30:38]
            y.Y_coord = line[38:46]
            y.Z_coord = line[46:54]
            PDBarray.append(y)
        else: 
            pass
    pdbin.close()
    return PDBarray