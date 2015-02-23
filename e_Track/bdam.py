# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 17:08:27 2014

@author: charlie
"""
import numpy as np
from classholder_v2 import StructurePDB
import sys
import os
from colrowsec_to_cartesian import findCartesianTranslationVectors
from PDBfile_manipulation import PDBtoCLASSARRAY_v2 as pdb2list

def pdbCUR_symgen(inputpdbfile,outputpdbfile,space_group):
    # function to run CCP4 program pdbCUR to generate all symmetry related atoms of 
    # original structure 
    print '---------------------------------------------------------------------------'
    print 'Determining symmetrically related atoms to original structure, using pdbCUR'

    input1 = "/Applications/ccp4-6.4.0/bin/pdbcur "+\
             "XYZIN %s " %(str(inputpdbfile))+\
             "XYZOUT %s " %(str(outputpdbfile))+\
             "SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib "

    # specify to remove hydrogen atoms and pick on the most probably conformation
    # (for two conformations with 0.5 occupancy, the first - A is chosen and occupancy
    # set to 1.00). Also remove all anisou info from file - since it is not needed for 
    # current analysis
    input2 = "delhydrogen\n"+\
             "mostprob\n"+\
             "noanisou\n"+\
             "symm %s\n" %(str(space_group))+\
             "GENUNIT\n"+\
             "mkchainIDs\n"+\
             "END"

    # run input2 in pdbcur
    textinput = open('inputfile.txt','w')
    textinput.write(input2)
    textinput.close()
    os.system(input1 +' < inputfile.txt')
    os.remove('inputfile.txt')

    print '\n---> pdbcur run complete...'
    print '--------------------------'
    print 'Summary:'
    print 'Input pdb file: %s' %(str(inputpdbfile))
    print 'Output pdb file: %s' %(str(outputpdbfile))

    # determine initial number of atoms in input pdb file
    pdbin = open(str(inputpdbfile),'r')
    pdbinlines = pdbin.readlines()
    counter = 0
    for line in pdbinlines:
        if 'ATOM' in line[0:5]:
            counter += 1
    pdbin.close()
    print 'number of atoms in input pdb file: %s' %(str(counter))

    # determine number of atoms in output pdb file
    pdbout = open(str(outputpdbfile),'r')
    pdboutlines = pdbout.readlines()
    counter = 0
    for line in pdboutlines:
        if 'ATOM' in line[0:5]:
            counter += 1
    pdbout.close()
    print 'number of atoms in output pdb file: %s' %(str(counter))
    print '--------------------------'



def translate26cells(inputpdbfile,outpdbfile):
    # function to create an extended pdb file containing all atoms translated from
    # the original structure into the adjacent 26 unit cells surrounding the unit
    # cell the structure is situated in
    print '------------------------------------------------'
    print 'Determining related atoms to original structure,\n'+\
          'translated into adjacent 26 unit cells'

    pdbin = open(str(inputpdbfile),'r')
    pdbout = open(str(outpdbfile),'w')

    for line in pdbin.readlines():
        if 'CRYST1' in line[0:6]:
            a = float(line.split()[1])
            b = float(line.split()[2])
            c = float(line.split()[3])
            alpha = float(line.split()[4])
            beta = float(line.split()[5])
            gamma = float(line.split()[6])
            pdbout.write(line)

            avec,bvec,cvec = findCartesianTranslationVectors(a,b,c,alpha,beta,gamma)

        if ('ATOM' in line[0:5] or 'HETATM' in line[0:6]):

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            for i in ([-1,0,1]):
                for j in ([-1,0,1]):
                    for k in ([-1,0,1]):
                        xyz_new = np.array([x,y,z]) + i*np.array(avec) + j*np.array(bvec) + k*np.array(cvec)
                        x_new = xyz_new[0]
                        y_new = xyz_new[1]
                        z_new = xyz_new[2]

                        pdbout.write(line[0:30])
                        pdbout.write("{0:8.3f}".format(x_new))
                        pdbout.write("{0:8.3f}".format(y_new))
                        pdbout.write("{0:8.3f}".format(z_new))
                        pdbout.write(line[54:80])
                        pdbout.write('\n')
    pdbin.close()
    pdbout.close()
    print '---> success!'



def restrict14A(PDBarrayin,expanded_inputfile,outputpdbfile):
    # function to restrict the pdb file containing neighbouring 
    # atoms of original pdb structure to 14 Angstroms around
    # original structure
    print '--------------------------------------------------------------'
    print 'Restricting atoms to within 14 Angstroms of original structure'

    # find max and min x values in original PDB file
    PDBarrayin.sort(key=lambda x: (x.X_coord))
    minX = PDBarrayin[0].X_coord
    maxX = PDBarrayin[len(PDBarrayin)-1].X_coord

    # find max and min y values in original PDB file
    PDBarrayin.sort(key=lambda x: (x.Y_coord))
    minY = PDBarrayin[0].Y_coord
    maxY = PDBarrayin[len(PDBarrayin)-1].Y_coord

    # find max and min z values in original PDB file
    PDBarrayin.sort(key=lambda x: (x.Z_coord))
    minZ = PDBarrayin[0].Z_coord
    maxZ = PDBarrayin[len(PDBarrayin)-1].Z_coord

    pdbin = open(str(expanded_inputfile),'r')
    pdbout = open(str(outputpdbfile),'w')
    for line in pdbin.readlines():
        if ('ATOM' in line[0:5] or 'HETATM' in line[0:6]):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if x - maxX < 14 and minX - x < 14:
                if y - maxY < 14 and minY - y < 14:
                    if z - maxZ < 14 and minZ - z < 14:
                        pdbout.write(line)
        else:
            pdbout.write(line)

    pdbin.close()
    pdbout.close()
    print '---> success!'



def numsurroundatoms_calculate(initialPDBfile,PDBarray,threshold):
    # function determines for each atom in structure the number of neighbouring atoms within 
    # a threshold (defined above) for all atoms. For each atom, number of contacts added
    # as class attribute for atom
    print '----------------------------------------------------'
    print 'Calculating contact number for atoms in structure...'

    # determine the correct extended pdb file, with atoms present up to 1 unit cell 
    # away from the original structure
    inputpdbfile1 = initialPDBfile

    # determine the space group for the input pdb file:
    pdbin = open(str(inputpdbfile1),'r')
    for line in pdbin.readlines():
        if 'CRYST1' in line[0:6]:
            space_group = line[55:66]
    pdbin.close

    # run the above functions to (a) determine the symmetrically related
    # atoms to the original structure, (b) translate to determine the 
    # location of all atoms within the adjacent 26 unit cells to the 
    # original structure, and (c) to restrict to atoms only within 14 
    # Angstroms of the original structure.
    outputpdbfile1 = initialPDBfile[:-4]+'_pdbCURsymgenOUT.pdb'
    pdbCUR_symgen(inputpdbfile1,outputpdbfile1,space_group)
    outputpdbfile2 = initialPDBfile[:-4]+'_translate26cells.pdb'
    translate26cells(outputpdbfile1,outputpdbfile2)
    extended_pdbfile = initialPDBfile[:-4]+'_restrict14A.pdb'
    restrict14A(PDBarray,outputpdbfile2,extended_pdbfile)

    counter = 0
    for atom in PDBarray:
        print counter
        counter += 1
        num_contacts = 0

        # read through extended pdb file and determine how many atoms are neighbours 
        # (within the specified threshold distance) 
        pdbin = open(extended_pdbfile,'r')
        for line in pdbin.readlines():
            if ('ATOM' in line[0:5] or 'HETATM' in line[0:6]):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # determine whether the atom specified in this line has coordinates
                # close to the current atom in PDBarray
                if abs(x - atom.X_coord) < threshold or abs(atom.X_coord - x) < threshold:
                    continue
                if abs(y - atom.Y_coord) < threshold or abs(atom.Y_coord - y) < threshold:
                    continue
                if abs(z - atom.Z_coord) < threshold or abs(atom.Z_coord - z) < threshold:
                    continue

                distance = np.sqrt(np.square(atom.X_coord - x) +
                                   np.square(atom.Y_coord - y) + 
                                   np.square(atom.Z_coord - z))

                if distance != 0 and distance < threshold:
                    num_contacts += 1
        # add the number of surrounding atoms as an attribute to the atom object            
        atom.numsurroundatoms = num_contacts
    print '---> success!'



def bdamage_calculate(PDBarray):
    # function to calculate Bdamage style metric for each atom, to save bdam 
    # attribute for each atom
    print '-------------------------------------------------------'
    print 'Calculating bdam style metric for atoms in structure...'
    for atom in PDBarray:
        simpacking_bfactors = []
        for otheratoms in PDBarray:
            if round(atom.numsurroundatoms/10) == round(otheratoms.numsurroundatoms/10):
                simpacking_bfactors.append(float(otheratoms.Bfactor))
        bdam = float(atom.Bfactor)/(np.mean(simpacking_bfactors))
        
        atom.bdam = bdam
    print '---> success!'



def numsurroundatms_extract(initialPDBarray,laterPDBarray):
    # function to extract numbers for surrounding atoms from intial pdb structure
    # and extend these values to the same atoms in later pdb structures (for the 
    # same damage series)

    # loop through the later dataset and assign the corresponding num of 
    # neighbouring atoms from the same atom in the initial dataset
    for atom in laterPDBarray:
        for otheratom in initialPDBarray:
            if (atom.atomtype == otheratom.atomtype and
               atom.basetype == otheratom.basetype and
               atom.residuenum == otheratom.residuenum and
               atom.chaintype == otheratom.chaintype):
                atom.numsurroundatoms = otheratom.numsurroundatoms





        

