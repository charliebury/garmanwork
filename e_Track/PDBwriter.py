# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 23:04:31 2015

@author: charlie
"""
from PDBfile_manipulation import convertPDBobject_toPDBline_fn
from datetime import datetime as dt


def densPDB_write(dataset_info,PDBarray):
    # WRITE A NEW PDB FILE FOR THE MEAN,MEDIAN,MAX,MIN 
    # LOCALISED ELECTRON DENSITY TO EACH ATOM - WRITTEN IN THE OCCUPANCY 
    # COLUMN FOR EACH ATOM:
    PDBout_mean = open(dataset_info.pdbname + '_mean' + '.pdb','w')
    PDBout_median = open(dataset_info.pdbname + '_median' + '.pdb','w')
    PDBout_max = open(dataset_info.pdbname + '_max' + '.pdb','w')
    PDBout_min = open(dataset_info.pdbname + '_min' + '.pdb','w')
    
    # first read through input pdb file and write all remark lines to 
    # output file
    pdbin = open(dataset_info.pdbname + '.pdb', "r")
    for line in pdbin.readlines():
            if line.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
                PDBout_mean.write(line)
                PDBout_median.write(line)
                PDBout_max.write(line)
                PDBout_min.write(line)
    pdbin.close()
    
    # now for each atom in pdb file, write the new pdb line with density 
    # metrics added
    for atom in PDBarray:
        mean_line = convertPDBobject_toPDBline_fn(atom,atom.meandensity)    
        median_line = convertPDBobject_toPDBline_fn(atom,atom.mediandensity)    
        max_line = convertPDBobject_toPDBline_fn(atom,atom.maxdensity)    
        min_line = convertPDBobject_toPDBline_fn(atom,atom.mindensity)
        PDBout_mean.write(mean_line + '\n')
        PDBout_median.write(median_line + '\n')
        PDBout_max.write(max_line + '\n')
        PDBout_min.write(min_line + '\n')    
        
    PDBout_mean.write('END')
    PDBout_median.write('END')
    PDBout_max.write('END')
    PDBout_min.write('END')
    PDBout_mean.close()
    PDBout_median.close()
    PDBout_max.close()
    PDBout_min.close() 
    
    
    
def denssummary_densorder(dataset_info,PDBarray):    
    # write an output file containing the mean/median/max/min density 
    # for each atom in file
    outputfile = str(dataset_info.pdbname) + '_edensperatom.txt'
    atom_edensity = open(str(outputfile),'w')
    
    # write some header information for output file:
    datenow = dt.now().strftime("%Y%m%d_%H%M%S")
    atom_edensity.write('#Datestamp: ' + str(datenow) + '\n')
    atom_edensity.write(
        '#Input structure file: ' + str(dataset_info.pdbname) + 
        '.pdb' + '\n')
    
    # to sort the PDBarray list of atom objects by meandensity (in order 
    # from most positive to most negative), use the following:
    PDBarray.sort(key=lambda x: x.meandensity, reverse=True)
    print 'Sorting by mean density:'
    for i in range(0,len(PDBarray)):   
        
        print 'Mean density: %s; Residue: %s; Chain: %s; Atom: %s'\
        %(str(PDBarray[i].meandensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype))
          
        atom_edensity.write(
        'Mean density: %s; Residue: %s; Chain: %s; Atom: %s \n'\
        %(str(PDBarray[i].meandensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype)))
          
    atom_edensity.write('...\n')
    
    # to sort the PDBarray list of atom objects by maxdensity (in order 
    # from most positive to most negative), use the following:
    PDBarray.sort(key=lambda x: x.maxdensity, reverse=True)
    print 'Sorting by max density:'
    for i in range(0,len(PDBarray)): 
        
        print 'Max density: %s; Residue: %s; Chain: %s; Atom: %s'\
        %(str(PDBarray[i].maxdensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype))
          
        atom_edensity.write(
        'Max density: %s; Residue: %s; Chain: %s; Atom: %s \n'\
        %(str(PDBarray[i].maxdensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype)))
          
    atom_edensity.write('...\n')
    
    # to sort the PDBarray list of atom objects by mindensity (in order 
    # from most positive to most negative), use the following:
    PDBarray.sort(key=lambda x: x.mindensity, reverse=True)
    print 'Sorting by min density:'
    for i in range(0,len(PDBarray)):    

        print 'Min density: %s; Residue: %s; Chain: %s; Atom: %s'\
        %(str(PDBarray[i].mindensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype))
          
        atom_edensity.write(
        'Min density: %s; Residue: %s; Chain: %s; Atom: %s \n'\
        %(str(PDBarray[i].mindensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype)))

    atom_edensity.write('...\n')
    
    #to sort the PDBarray list of atom objects by mediandensity (in order 
    # from most positive to most negative), use the following:
    PDBarray.sort(key=lambda x: x.mediandensity, reverse=True)
    print 'Sorting by median density:'
    for i in range(0,len(PDBarray)):    

        print 'Median density: %s; Residue: %s; Chain: %s; Atom: %s'\
        %(str(PDBarray[i].mediandensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype))
          
        atom_edensity.write(
        'Median density: %s; Residue: %s; Chain: %s; Atom: %s \n'\
        %(str(PDBarray[i].mediandensity),str(PDBarray[i].basetype),
          str(PDBarray[i].chaintype),str(PDBarray[i].atomtype)))

    atom_edensity.write('...\n')
    
    atom_edensity.close()


    
def denssummary_atomorder(dataset_info,PDBarray):
    # Want to write a .txt file with mean, median,max,min density for 
    # each atom in structure, ordered by atom number in structure:
    summaryfilename = str(dataset_info.pdbname) + 'summary_edensity.txt'
    summary_edensityfile = open(str(summaryfilename),'w')
    
    # write some header information for output file:
    datenow = dt.now().strftime("%Y%m%d_%H%M%S")
    summary_edensityfile.write('#Datestamp: ' + str(datenow) + '\n')
    summary_edensityfile.write(
        '#Input structure file: ' + str(dataset_info.pdbname) +
        '.pdb' + '\n')    
    summary_edensityfile.write('Atom Number \t Mean \t Median \t Min \t Max \n')
    
    PDBarray.sort(key=lambda x: x.atomnum)
    for i in range(0,len(PDBarray)):    
        summary_edensityfile.write(
            str(PDBarray[i].atomnum) + '\t' + str(PDBarray[i].meandensity) + 
            '\t' + str(PDBarray[i].mediandensity) + '\t' + str(PDBarray[i].mindensity) + 
            '\t' + str(PDBarray[i].maxdensity) + '\n')
    summary_edensityfile.close()    
    
    
