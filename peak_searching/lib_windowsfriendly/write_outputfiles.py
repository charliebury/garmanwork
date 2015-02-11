# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:42:12 2014

@author: charlie
"""

from math import sqrt, pow
from file_len_function import file_len

##############################################################################################################################################################
##############################################################################################################################################################
#SECTION: Write new .txt file containing processed information - listing specific atom details from the structure associated with each selected specific damage sites
##############################################################################################################################################################
##############################################################################################################################################################

# A class for damage sites associated with atoms in structure
class damagedatoms:
    def __init__(self,atomnum="",residuenum="",atomtype="",basetype="",chaintype="",
        X_coord="",Y_coord="",Z_coord="",atomidentifier="",mindist="",peaknum="",
        peakvalue="",Bfactor="",Occupancy=""):
            
        self.atomnum = atomnum
        self.residuenum = residuenum
        self.atomtype = atomtype
        self.basetype = basetype 
        self.chaintype = chaintype
        self.X_coord = X_coord
        self.Y_coord = Y_coord
        self.Z_coord = Z_coord
        self.atomidentifier = atomidentifier
        self.mindist = mindist
        self.peaknum = peaknum
        self.peakvalue = peakvalue
        self.Bfactor = Bfactor
        self.Occupancy = Occupancy


def write_outputfiles(datedir,PDBarray):

    #read through most recent .pdb peak file generated
    diffpeakpdb = str("edittednegpeaksonly_thresholddistance"+datedir.date+".pdb")
    diffpdbin = open(datedir.path + '/' + str(diffpeakpdb), "r")
    #write a .txt file here to complete filtering script
    distfile = open(datedir.path + '/' + "distancefile_"+datedir.date+".txt", "w")
    
    distfile.write('peaknumber' + "\t" + 'distmetric' + "\t" + 'atomnum' + "\t" + 'residuenum' + "\t" + 'atomtype' + "\t" + 'basetype' + "\t" + 'chaintype' + "\t" + 'peakvalue' + "\t" + 'X_coord' + "\t" + 'Y_coord' + "\t" + 'Z_coord' + "\t" +"\n")
    diffpeaklines = diffpdbin.readlines()
    
    #read through most recent .ha file generated 
    hafile = str("edittednegpeaksonly_thresholddistance"+datedir.date+".ha")
    hafileread = open(datedir.path + '/' + str(hafile), "r")
    
    #convert .ha file to array form for easier manipulation
    haarray = []
    for line in hafileread:
        haarray.append(line)
    hafileread.close()
    
    counter = 0
    damagedatoms_list = []
    for line in diffpeaklines[4:file_len(datedir.path + '/' + str(diffpeakpdb))-1]:
        counter += 1
        
        #specify the x,y,z cartesian coordinates for each filtered damage site in 
        #the newest .ha file
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
    
        proDNAcoords = open(datedir.path+'/'+"pdbcoords.txt", "r")
        proDNAlines = proDNAcoords.readlines()
        Dvector = []
        
        #search through the pdb structure file and assign the difference map peak 
        #to the closest atom of the structure
        for pdbline in proDNAlines:
            xi = float(pdbline[30:38].strip())
            yi = float(pdbline[38:46].strip())
            zi = float(pdbline[46:54].strip())
            Di = sqrt( pow(x-xi,2) + pow(y-yi,2) + pow(z-zi,2))
            Dvector.append(Di)
    
        mindist = min(Dvector)
        minIndex = Dvector.index(min(Dvector))
        minPDBatom = minIndex + 1
        
        peakvalue = str(haarray[counter + 1].split()[5])
    
        y = damagedatoms(str(PDBarray[minIndex][6:11]).strip(),str(PDBarray[minIndex][22:26]).strip(),
                         str(PDBarray[minIndex][12:16]).strip(),str(PDBarray[minIndex][17:20]).strip(),
                         str(PDBarray[minIndex][21]).strip(),str(PDBarray[minIndex][30:38]).strip(),
                         str(PDBarray[minIndex][38:46]).strip(),str(PDBarray[minIndex][46:54]).strip(),
                         str(PDBarray[minIndex][76:78]).strip(),mindist,counter+1,peakvalue)


        proDNAcoords.close()

        distfile.write(str(y.peaknum) + "\t" + str(y.mindist) + "\t" + str(y.atomnum) +\
                       "\t" + str(y.residuenum) + "\t" + str(y.atomtype) + "\t" +\
                       str(y.basetype) + "\t" + str(y.chaintype) + "\t" + str(y.peakvalue) +\
                       "\t" + str(y.X_coord) + "\t" + str(y.Y_coord) + "\t" + str(y.Z_coord) + "\n")
    
        damagedatoms_list.append(y)

    diffpdbin.close()
    distfile.close()

    return damagedatoms_list


    ###############################################################################