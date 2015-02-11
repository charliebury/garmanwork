# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 01:48:33 2014

@author: charlie
"""

#A class for structure PDB file
class StructurePDB:
    def __init__(self,atomnum=0,residuenum=0,atomtype="",basetype="",chaintype="",X_coord=0,Y_coord=0,Z_coord=0,Bfactor=0,Occupancy=0):
    #def __init__(self) :   
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