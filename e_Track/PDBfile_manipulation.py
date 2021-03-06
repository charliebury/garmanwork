# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 23:39:24 2014
@author: charlie
"""
import sys
from classholder_v2 import StructurePDB,singlePDB,multiPDB
import numpy as np
from operator import sub,truediv
from random import randint
from progbar import progress


###############################################################################
def PDBtoCLASSARRAY(pdbfilename,PDBarray):
    #this function inputs a pdb file name and returns an array of pdb objects, 
    #organised following the StructurePDB class 
    
    pdbin = open(str(pdbfilename), "r")
    lines = pdbin.readlines()
    print 'Reading PDB file and converting to array of objects'
    for line in lines:
        if 'ATOM' in line.split()[0]:
            y = singlePDB(StructurePDB)
            y.atomnum = int(line.split()[1])
            y.atomtype = line.split()[2]
            y.basetype = line.split()[3]
            y.chaintype = line.split()[4]
            y.residuenum = int(line.split()[5])
            y.X_coord = float(line.split()[6])
            y.Y_coord = float(line.split()[7])
            y.Z_coord = float(line.split()[8])
            
            #it has been observed that in some cases the occupancy and Bfactor pdb
            #columns are merged (due to the Bfactor value being over 100). To address
            #this problem the following splitting is used:

            #Case1: no merging present:
            try:
                y.atomidentifier = line.split()[11]
                y.Occupancy = line.split()[9]
                y.Bfactor = line.split()[10]

            #Case2: merging present:    
            except IndexError:
                #print 'Occupancy and Bfactor columns merged in .pdb file -> splitting them now'
                y.atomidentifier = line.split()[10]
                y.Occupancy = str(line.split()[9])[0:4]
                y.Bfactor = str(line.split()[9])[4:len(line.split()[9])]
                #print 'Occupancy: ' +str(y.Occupancy) + '; Bfactor: ' + str(y.Bfactor) 
            
            PDBarray.append(y)
        else: 
            pass
    pdbin.close()
    return PDBarray
###############################################################################



###############################################################################
def PDBtoCLASSARRAY_v2(pdbfilename,PDBarray):
    # this function inputs a pdb file name and returns an array of pdb objects, 
    # organised following the StructurePDB class. It is the same as the above 
    # function (PDBtoCLASSARRAY) but directly uses the PDB official line space
    # formatting
    
    pdbin = open(str(pdbfilename), "r")
    lines = pdbin.readlines()
    print 'Reading PDB file and converting to array of objects'
    for line in lines:
        if ('ATOM' in str(line[0:6])) or ('HETATM' in str(line[0:6])):
            y = singlePDB(StructurePDB)
            y.atomnum = int(line[6:11].strip())
            y.atomtype = str(line[12:16].strip())
            y.basetype = str(line[17:20].strip())                       
            y.chaintype = str(line[21])                     
            y.residuenum = int(line[22:26].strip())             
            y.X_coord = float(line[30:38].strip())                       
            y.Y_coord = float(line[38:46].strip())                                             
            y.Z_coord = float(line[46:54].strip())                                                 
            y.Occupancy = str(line[54:60].strip())                                                    
            y.Bfactor = str(line[60:66].strip())
            y.atomidentifier = str(line[76:78].strip())  
            PDBarray.append(y)
        else: 
            pass
    pdbin.close()
    return PDBarray
###############################################################################



###############################################################################
def multiARRAY(data_list):
    #this function inputs a list of lists of PDB atom objects (see StructurePDB class)
    #and formats as an object of the class 'multiPDB'
     
    #first check that each PDBarray contains the same number of atoms (consistency check)
    if len(data_list) > 1:
        print 'Multiple datasets detected...'
        for dataset in data_list:
            if len(dataset.PDBDENSITY) != len(data_list[0].PDBDENSITY):
                print 'Not all PDB structures have same number of atoms -> stopping script'
                sys.exit()
    elif len(data_list) == 1:
        print 'Single dataset detected...'
        
    PDBdoses = []
    for atom in data_list[0].PDBDENSITY:
        y = multiPDB()
        y.atomnum = atom.atomnum
        y.residuenum = atom.residuenum
        y.atomtype = atom.atomtype
        y.basetype = atom.basetype 
        y.chaintype = atom.chaintype
        y.X_coord = atom.X_coord
        y.Y_coord = atom.Y_coord
        y.Z_coord = atom.Z_coord
        y.atomidentifier = atom.atomidentifier
        y.numsurroundatoms = atom.numsurroundatoms

        Bfactor_comb = []
        Occupancy_comb = []
        meandensity_comb = []
        maxdensity_comb = []
        mindensity_comb = []
        mediandensity_comb = []
        bdam_comb = []
        bdamchange_comb = []
        Bfactorchange_comb = []
        
        for dataset in data_list:
            Bfactor_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].Bfactor)
            Occupancy_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].Occupancy)
            meandensity_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].meandensity)
            maxdensity_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].maxdensity)
            mindensity_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].mindensity)
            mediandensity_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].mediandensity)
            bdam_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].bdam)
            bdamchange_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].bdamchange)
            Bfactorchange_comb.append(dataset.PDBDENSITY[int(atom.atomnum)-1].Bfactorchange)
            
        y.Bfactor = Bfactor_comb 
        y.Occupancy = Occupancy_comb
        y.meandensity = meandensity_comb
        y.maxdensity = maxdensity_comb
        y.mindensity = mindensity_comb
        y.mediandensity = mediandensity_comb
        y.bdam = bdam_comb
        y.bdamchange = bdamchange_comb
        y.Bfactorchange = Bfactorchange_comb
           
        PDBdoses.append(y)
    
    return PDBdoses
###############################################################################


###############################################################################
def multiARRAY_PDBonly(data_list):
    # this function inputs a list of lists of PDB atom objects (see 
    # StructurePDB class)and formats as an object of the class 'multiPDB'. 
    # Same as function above, but takes in list of atom lists (NOT 
    # datainfo objects as above)
     
    # first check that each PDBarray contains the same number of atoms 
    # (consistency check)
    if len(data_list) > 1:
        print 'Multiple datasets detected...'
        for dataset in data_list:
            if len(dataset) != len(data_list[0]):
                print 'Not all PDB structures have same number of atoms'\
                '-> stopping script'
                sys.exit()
    elif len(data_list) == 1:
        print 'Single dataset detected...'
        
    PDBdoses = []
    for atom in data_list[0]:
        y = multiPDB()
        y.atomnum = atom.atomnum
        y.residuenum = atom.residuenum
        y.atomtype = atom.atomtype
        y.basetype = atom.basetype 
        y.chaintype = atom.chaintype
        y.X_coord = atom.X_coord
        y.Y_coord = atom.Y_coord
        y.Z_coord = atom.Z_coord
        y.atomidentifier = atom.atomidentifier
        y.numsurroundatoms = atom.numsurroundatoms

        Bfactor_comb = []
        Occupancy_comb = []
        meandensity_comb = []
        maxdensity_comb = []
        mindensity_comb = []
        mediandensity_comb = []
        bdam_comb = []
        bdamchange_comb = []
        Bfactorchange_comb = []
        
        # for each damaged pdb file, find the atom above (from its type,
        # chain type, residue num etc, and then append the Bfactor etc.
        # to this atom). Quite slow, but ensures atoms are matched between
        # datasets if orders of atoms in the list happens to be switched
        for dataset in data_list:
            for otheratom in dataset:
                if (atom.residuenum == otheratom.residuenum and
                    atom.basetype == otheratom.basetype and
                    atom.chaintype == otheratom.chaintype and
                    atom.atomtype == otheratom.atomtype):
            
                    Bfactor_comb.append(otheratom.Bfactor)
                    Occupancy_comb.append(otheratom.Occupancy)
                    meandensity_comb.append(otheratom.meandensity)
                    maxdensity_comb.append(otheratom.maxdensity)
                    mindensity_comb.append(otheratom.mindensity)
                    mediandensity_comb.append(otheratom.mediandensity)
                    bdam_comb.append(otheratom.bdam)
                    bdamchange_comb.append(otheratom.bdamchange)
                    Bfactorchange_comb.append(otheratom.Bfactorchange)
        
        # to check that every dataset has contributed to the list as 
        # intended. If not, do not append this atom to the list. This 
        # means that solvent atoms which may be present in only a subset
        # of the datasets will be removed 
        if len(Bfactor_comb) != len(data_list):
            print 'Atom detected not present in all datasets'
            continue
        
        # if atom present in all datasets, give the atom its dose-dependent
        # attributes as specified below
        y.Bfactor = Bfactor_comb 
        y.Occupancy = Occupancy_comb
        y.meandensity = meandensity_comb
        y.maxdensity = maxdensity_comb
        y.mindensity = mindensity_comb
        y.mediandensity = mediandensity_comb
        y.bdam = bdam_comb
        y.bdamchange = bdamchange_comb
        y.Bfactorchange = Bfactorchange_comb
           
        # add this atom to the list   
        PDBdoses.append(y)
    
    return PDBdoses
###############################################################################



###############################################################################
def multiARRAY_diffatomnumbers(data_list):
    # this function inputs a list of lists of PDB atom objects (see StructurePDB class)
    # and formats as an object of the class 'multiPDB'. It is a variant of the function 
    # above which can also cope with structures containing different numbers of atoms 
    # (say if solvent molecules/ligands are included in a subset of the structures). 
    # In this case, the smallest common substructure between all structures will be used
     
    #first check that each PDBarray contains the same number of atoms (consistency check)
    if len(data_list) > 1:
        print 'Multiple datasets detected...'
        for dataset in data_list:
            if len(dataset) != len(data_list[0]):
                print 'Not all PDB structures have same number of atoms!'\
                ' Will only include atoms common to ALL structures...'
    elif len(data_list) == 1:
        print 'Single dataset detected...'
    
    # find the reference point atom in each dataset that will be used to
    # normalise the density metrics later on
    print '------------------------------------------------'
    print 'Locating reference point atom in all datasets...'

    # generate a random reference atom in the first dataset and check that
    # present in all later datasets
    found = 'NO'
    while found == 'NO':
        randnum = randint(0,len(data_list[0]))
        randatm = data_list[0][randnum]

        counter = 0
        for dataset in data_list:
            for atom in dataset:
                if (atom.atomtype == randatm.atomtype and
                   atom.basetype == randatm.basetype and
                   atom.chaintype == randatm.chaintype and
                   atom.residuenum == randatm.residuenum):
                    counter += 1

        if counter == len(data_list):
            found = 'YES'
            print 'Found random atom present in all datasets to use as reference'+\
                  ' atom:\n---> %s %s %s %s' %(str(randatm.atomtype),
                    str(randatm.basetype),str(randatm.chaintype),str(randatm.residuenum))


    refpoint_mean,refpoint_median,refpoint_min,refpoint_max = [],[],[],[]
    counter = 0
    for dataset in data_list:
        counter += 1
        for atom in dataset:       
            if (atom.atomtype == randatm.atomtype and
                atom.basetype == randatm.basetype and
                atom.chaintype == randatm.chaintype and
                atom.residuenum == randatm.residuenum):
                    
                refpoint_mean.append(atom.meandensity)
                refpoint_median.append(atom.mediandensity)
                refpoint_min.append(atom.mindensity)
                refpoint_max.append(atom.maxdensity)
                print 'Reference atom found in dataset %s' %str(counter)
                print 'Atom number: %s' %str(atom.atomnum)
    
    # check that referenced atom found in all datasets
    if len(refpoint_mean) != len(data_list):        
        print 'Cannot find reference point atom in all datasets'
        print '---> terminating script...'
        sys.exit()
    else:
        print '---> success!'
   
    PDBdoses = []
    notincludedatmcounter = 0
    
    print '------------------------------------------------'
    print 'Locating common atoms to ALL datasets...'
    print 'SUMMARY:'

    i = 0
    num_atoms = len(data_list[0])

    for atom in data_list[0]:

        # unessential loading bar add-in
        i += 1
        progress(i, num_atoms, suffix='')

        atm_counter = 1
        
        Bfactor_comb = [atom.Bfactor]
        Occupancy_comb = [atom.Occupancy]
        meandensity_comb = [atom.meandensity]
        maxdensity_comb = [atom.maxdensity]
        mindensity_comb = [atom.mindensity]
        mediandensity_comb = [atom.mediandensity]
        bdam_comb = [atom.bdam]
        numvoxels_comb = [atom.numvoxels]
        
        # list of index of same atom in each later dataset
        indexindataset = []

        # check whether atom in all datasets:
        for dataset in data_list[1:]:
            k = -1        
            for otheratom in dataset: 
                k += 1  
                if (atom.residuenum == otheratom.residuenum and 
                    atom.atomtype == otheratom.atomtype and 
                    atom.basetype == otheratom.basetype and 
                    atom.chaintype == otheratom.chaintype):
                    
                    atm_counter += 1       
                    Bfactor_comb.append(otheratom.Bfactor)
                    Occupancy_comb.append(otheratom.Occupancy)
                    meandensity_comb.append(otheratom.meandensity)
                    maxdensity_comb.append(otheratom.maxdensity)
                    mindensity_comb.append(otheratom.mindensity)
                    mediandensity_comb.append(otheratom.mediandensity)
                    bdam_comb.append(otheratom.bdam)
                    numvoxels_comb.append(otheratom.numvoxels)
                    indexindataset.append(k)
                    break

        # remove this located atom from the later dataset now that it
        # has been located --> to make loop faster 
        for j in range(1,len(indexindataset)+1):
            if indexindataset[j-1] != -1:
                data_list[j].pop(indexindataset[j-1])
                        
        if atm_counter != len(data_list):
            print 'Atom not found in all datasets!'
            print 'res num: %s, atm type: %s, res type: %s, chain: %s '\
            %(atom.residuenum,atom.atomtype,atom.basetype,atom.chaintype)
            print '---> not including atom...'
            notincludedatmcounter += 1
            continue
      
        else:                 
            y = multiPDB()
            y.atomnum = atom.atomnum
            y.residuenum = atom.residuenum
            y.atomtype = atom.atomtype
            y.basetype = atom.basetype 
            y.chaintype = atom.chaintype
            y.X_coord = atom.X_coord
            y.Y_coord = atom.Y_coord
            y.Z_coord = atom.Z_coord
            y.atomidentifier = atom.atomidentifier
            y.numsurroundatoms = atom.numsurroundatoms
            
            # calculate the normalised mean,median,max,min for each dataset
            y.meandensity_norm = map(truediv, meandensity_comb, refpoint_mean)
            y.mediandensity_norm = map(truediv, mediandensity_comb, refpoint_median)
            y.mindensity_norm = map(truediv, mindensity_comb, refpoint_min)
            y.maxdensity_norm = map(truediv, maxdensity_comb, refpoint_max)

            y.Bfactor = Bfactor_comb 
            y.Occupancy = Occupancy_comb
            y.meandensity = meandensity_comb
            y.maxdensity = maxdensity_comb
            y.mindensity = mindensity_comb
            y.mediandensity = mediandensity_comb
            y.bdam = bdam_comb
            y.numvoxels = numvoxels_comb
               
            PDBdoses.append(y)
    print '\n------------------------------------------------'    
    print 'Number of atoms removed since not in all datasets: %s' %str(notincludedatmcounter)
    print '---> Finished!'
    return PDBdoses
###############################################################################



###############################################################################
def convertPDBobject_toPDBline_fn(element,Occupancy):
    #script to convert atom information (in class format) to 'ATOM' line format for
    #pdb output files. Note that the following PDB ATOM line convention is used:
    
    #FIELD      #COLUMNS        DATA TYPE       CONTENTS                            
    #--------------------------------------------------------------------------------
    #FIELD1      1 -  6        Record name     "ATOM  "                                            
    #FIELD2      7 - 11        Integer         Atom serial number.                   
    #FIELD3     13 - 16        Atom            Atom name.                            
    #FIELD4     17             Character       Alternate location indicator.         
    #FIELD5     18 - 20        Residue name    Residue name.                         
    #FIELD6     22             Character       Chain identifier.                     
    #FIELD7     23 - 26        Integer         Residue sequence number.              
    #FIELD8     27             AChar           Code for insertion of residues.       
    #FIELD9     31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
    #FIELD10    39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
    #FIELD11    47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
    #FIELD12    55 - 60        Real(6.2)       Occupancy.                            
    #FIELD13    61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
    #FIELD14    73 - 76        LString(4)      Segment identifier, left-justified.   
    #FIELD15    77 - 78        LString(2)      Element symbol, right-justified.      
    #FIELD16    79 - 80        LString(2)      Charge on the atom.   
    
    #NOTE THAT OCCUPANCY IS NOT SPECIFIC TO THE 'ELEMENT' CHOSEN - HENCE IT CAN 
    #BE ASSIGNED TO ANOTHER METRIC IF DESIRED
    
    #for 'atom' of PDBarray (ie a specific atom of the structure):
    FIELD1 = "ATOM  " #has length 6
    FIELD2 = " "*(5-len(str(element.atomnum))) + str(element.atomnum) #has length 5
    BREAK1 = " "*1 #has length 2
    FIELD3 = str(element.atomtype) + " "*(4-len(str(element.atomtype))) #has length 4
    FIELD4 = " " #has length 1
    FIELD5 = str(element.basetype).rjust(3) #has length 3
    BREAK2 = " " # has length 1
    FIELD6 = str(element.chaintype)
    FIELD7 = " "*(4-len(str(element.residuenum))) + str(element.residuenum) #has length 4
    FIELD8 = " " #has length 1
    BREAK3 = " "*3 #has length 3
    FIELD9 = "{0:8.3f}".format(element.X_coord) #has length 8
    FIELD10 = "{0:8.3f}".format(element.Y_coord) #has length 8
    FIELD11 = "{0:8.3f}".format(element.Z_coord) #has length 8
    FIELD12 = "{0:6.2f}".format(float(Occupancy)) #has length 6
    FIELD13 = "{0:6.2f}".format(float(element.Bfactor)) #has length 6
    BREAK4 = " "*6 # has length 6
    FIELD14 = " "*4 #has length 4
    FIELD15 = str(element.atomidentifier).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16
    
    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line
###############################################################################
    

###############################################################################
def convertPDBobject_toPDBline_fn_Bfactorsetter(element,Bfactor):
    #script to convert atom information (in class format) to 'ATOM' line format for
    #pdb output files. Note that the following PDB ATOM line convention is used:
        
    #NOTE THAT Bfactor IS NOT SPECIFIC TO THE 'ELEMENT' CHOSEN - HENCE IT CAN 
    #BE ASSIGNED TO ANOTHER METRIC IF DESIRED
    
    #for 'atom' of PDBarray (ie a specific atom of the structure):
    FIELD1 = "ATOM  " #has length 6
    FIELD2 = " "*(5-len(str(element.atomnum))) + str(element.atomnum) #has length 5
    BREAK1 = " "*1 #has length 2
    FIELD3 = str(element.atomtype) + " "*(4-len(str(element.atomtype))) #has length 4
    FIELD4 = " " #has length 1
    FIELD5 = str(element.basetype).rjust(3) #has length 3
    BREAK2 = " " # has length 1
    FIELD6 = str(element.chaintype)
    FIELD7 = " "*(4-len(str(element.residuenum))) + str(element.residuenum) #has length 4
    FIELD8 = " " #has length 1
    BREAK3 = " "*3 #has length 3
    FIELD9 = "{0:8.3f}".format(element.X_coord) #has length 8
    FIELD10 = "{0:8.3f}".format(element.Y_coord) #has length 8
    FIELD11 = "{0:8.3f}".format(element.Z_coord) #has length 8
    FIELD12 = "{0:6.2f}".format(float(element.Occupancy)) #has length 6
    FIELD13 = "{0:6.2f}".format(float(Bfactor)) #has length 6
    BREAK4 = " "*6 # has length 6
    FIELD14 = " "*4 #has length 4
    FIELD15 = str(element.atomidentifier).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16
    
    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line
###############################################################################
