# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 20:04:30 2015

@author: charlie
"""
from PDBfile_manipulation import PDBtoCLASSARRAY_v2 as pdb2list
from PDBfile_manipulation import convertPDBobject_toPDBline_fn_Bfactorsetter as switchBfactor
import sys

# the name of the initial pdb structure (lowest damage)
pdbname1 = '/Users/charlie/DPhil/YEAR1/NOV/PDB_submission_frompaper/phenixrefined_models/dataset1/ROTATION1_GARMAN2014_refine_41'

# the name of the damaged pdb structure
pdbname2 = '/Users/charlie/DPhil/YEAR1/JAN/sd_finder/pdb_files/3clcdamage5'

# the name of the file to be created
pdbname3 = '/Users/charlie/Desktop/3clcdamage1withdamage5Bfactor'
pdbout = open(str(pdbname3)+'.pdb','w')
pdbin = open(str(pdbname1)+'.pdb','r')

PDBarray1 = pdb2list(pdbname1+'.pdb',[])
PDBarray2 = pdb2list(pdbname2+'.pdb',[])


for line in pdbin.readlines():
     if ('ATOM' not in str(line[0:6])) and ('HETATM' not in str(line[0:6])):
         pdbout.write(line)
     else:
         break
     
counter = -1     
for atom in PDBarray1: 
    counter += 1
    # check that the atom is the same in both files and add line if so   
    try:
        PDBarray2[counter]
    except IndexError:
        break
        
    if atom.atomnum == PDBarray2[counter].atomnum and atom.atomtype == PDBarray2[counter].atomtype:
         editted_line = switchBfactor(atom,PDBarray2[counter].Bfactor)
         pdbout.write(editted_line+'\n')
         print editted_line
         counter + 1
    else:
         print 'discrepancy between two atoms which should be the same identity. Showing pdb lines...'
         print 'Initial dataset: '
         print switchBfactor(atom,atom.Bfactor)
         print 'Later dataset: '
         print switchBfactor(PDBarray2[counter],PDBarray2[counter].Bfactor)
         response = raw_input('Want to include? (y/n) --> will use initial dataset Bfactor...')
         if response in ('y'):
             pdbout.write(switchBfactor(atom,atom.Bfactor)+'\n')
         elif response in ('n'):
             pass
         else:
             print 'Response not recognised... --> terminating script'
             sys.exit()

pdbout.write('END')

pdbin.close()
pdbout.close()     
