# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 11:28:50 2015

@author: charlie
"""
from delete_listindices import multi_delete 
from map2voxelclasslist_nocartconversion import densmap2class
import sys   
from PDBfile_manipulation import PDBtoCLASSARRAY_v2 as pdb2list
from densityanalysisplots import edens_scatter
from res_formatter import densper_resatom_NOresidueclass,densper_res
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import izip as zip, count
from densandatommaps2vxllist import combinevxl_atmanddensvals,combinevxl_atmanddensvals_gen
from vxlsperatm_analysisplots import vxlsperatm_hist,vxlsperatm_kde
import math


###----------------###----------------###----------------###----------------###
###############################################################################      
# A class for .map file voxel
class voxel_density:
    def __init__(self,density=0,atmnum=0):
        self.density = density
        self.atmnum = atmnum

# A class for collecting voxels per atom
class vxls_of_atm:
    def __init__(self,vxls=[],atmnum = 0):
        self.vxls = vxls
        self.atmnum = atmnum
###############################################################################
###----------------###----------------###----------------###----------------###

def maps2atmdensity(where,pdbname,mapfilname1,maptype1,mapfilname2,maptype2):
    # title
    print '\n'
    print '================================================'
    print '------------------------------------------------'
    print '|||              eTrack run                  |||'
    print '------------------------------------------------'
    print '================================================'
    print '\n'
    sys.stdout.flush()
    
    print 'Reading in pdb file...'
    sys.stdout.flush()
    # next read in the pdb structure file:
    # run function to fill PDBarray list with atom objects from structure
    PDBarray = pdb2list(where+pdbname+'.pdb',[])
    print '---> success'
    sys.stdout.flush()
    
    # want to make sure array of structure atoms ordered by atomnumber
    # before reading through them
    PDBarray.sort(key=lambda x: x.atomnum)
       
    # need to get VDW radius for each atom:
    for atom in PDBarray:
        atom.VDW_get()  
    
    # find number of atoms in structure
    num_atoms = len(PDBarray)
    
    # read in the atom map
    print '------------------------------------------------'
    print 'Atom map read in: '
    sys.stdout.flush()
    atmmap,atom_indices = densmap2class(where+mapfilname1,maptype1,[])
    
    # find atom numbers present in list (repeated atom numbers removed)
    seen = set()
    seen_add = seen.add
    uniq_atms = [x for x in atmmap.vxls_val if not (x in seen or seen_add(x))] 
    
    # find set of atoms numbers not present (i.e atoms not assigned to voxels)
    Atms_notpres = set(range(1,num_atoms+1)) - set(uniq_atms)
    print 'Number of atoms not assigned to voxels: %s' %str(len(Atms_notpres))
    sys.stdout.flush()
    
    # read in the density map
    print '------------------------------------------------'
    print 'Density map read in: '
    sys.stdout.flush()
    densmap = densmap2class(where+mapfilname2,maptype2,atom_indices)  
    
    
    print '------------------------------------------------'
    print 'Checking that maps have same dimensions and sampling properties...' 
    sys.stdout.flush()
    # Checks that the maps have the same dimensions, grid sampling etc.
    # Note that for the cell dimensions, 3dp is take only since unpredictable
    # fluctuations between the .mtz (and thus .map) and pdb recorded
    # cell dimensions have been observed in trial datasets
    if ('%.3f' %atmmap.celldim_a != '%.3f' %densmap.celldim_a or
    '%.3f' %atmmap.celldim_b != '%.3f' %densmap.celldim_b or
    '%.3f' %atmmap.celldim_c != '%.3f' %densmap.celldim_c or
    atmmap.celldim_alpha != densmap.celldim_alpha or
    atmmap.celldim_beta != densmap.celldim_beta or 
    atmmap.celldim_gamma != densmap.celldim_gamma or
    atmmap.fast_axis != densmap.fast_axis or
    atmmap.med_axis != densmap.med_axis or
    atmmap.slow_axis != densmap.slow_axis or
    atmmap.gridsamp1 != densmap.gridsamp1 or
    atmmap.gridsamp2 != densmap.gridsamp2 or 
    atmmap.gridsamp3 != densmap.gridsamp3 or
    atmmap.nx != densmap.nx or
    atmmap.ny != densmap.ny or
    atmmap.nz != densmap.nz or
    atmmap.start1 != densmap.start1 or
    atmmap.start2 != densmap.start2 or
    atmmap.start3 != densmap.start3 or
    atmmap.type != densmap.type):
        print 'Incompatible map properties --> terminating script'
        sys.exit()
    else:
        print '---> success: The atom and density map are of compatible format!'
    

    print '------------------------------------------------'      
    print 'Total number of voxels assigned to atoms: %s' %str(len(atmmap.vxls_val))
    sys.stdout.flush()


    # create list of voxel objects in class voxel_density 
    print '------------------------------------------------'
    print 'Combining voxel density and atom values...'
    sys.stdout.flush()
    vxl_list = combinevxl_atmanddensvals(atmmap.vxls_val,densmap.vxls_val)
    
    # delete atmmap and densmap now to save memory
    densmap =[]
    atmmap = []
    
    # next assign voxels to atom objects in PDBarray list
    print '\n------------------------------------------------'
    print 'Assigning voxels to corresponding atoms...'
    sys.stdout.flush()
    # first order voxels by atom number
    vxl_list.sort(key=lambda x: x.atmnum)
    
    # initialise list of voxels for each atom before can append voxels. Note 
    # that voxels are saved in list attribute of local vxls_of_atom class objects
    PDBarray.sort(key=lambda x: x.atomnum)   
    vxlsperatom = []
    pdbappend = vxlsperatom.append

    for i in range(1,PDBarray[len(PDBarray)-1].atomnum + 1):
        counter = 0
        for atom in PDBarray:
            if i == atom.atomnum:
                y = vxls_of_atm([],i)
                counter += 1
        if counter == 0:
            y = vxls_of_atm(['atom not present'],i)
        pdbappend(y)

    # unessential loading bar add-in (part1)
    total = len(vxl_list)
    point = total / 100
    increment = total / 100

    counter = -1  
    for vxl in vxl_list:
        counter += 1
        (vxlsperatom[vxl.atmnum-1].vxls).append(vxl)

        # check here
        if vxlsperatom[vxl.atmnum-1].vxls[0] == 'atom not present':
            print 'something has gone wrong...'
            sys.exit()

        # unessential loading bar add-in (part2)
        if(counter % (1 * point) == 0):
            sys.stdout.write(
                    "\r[" + "=" * (counter / increment) +
                    " " * ((total - counter)/ increment) + 
                    "]" +  str(counter / point) + "%")
            sys.stdout.flush()

    # can remove vxl_list variable now to save memory
    vxl_list = []    
        
    # histogram plot of number of voxels per atom
    vxlsperatm_hist(pdbname,where,vxlsperatom)
    vxlsperatm_kde(pdbname,where,vxlsperatom)    
    
    # determine min,max,mean,median density per atom
    print '\n------------------------------------------------'
    print 'Calculating electron density statistics per atom...'
    sys.stdout.flush()
    for atom in PDBarray:
        if len(vxlsperatom[atom.atomnum-1].vxls) != 0:
            voxelsofatom = [vxl.density for vxl in vxlsperatom[atom.atomnum-1].vxls]
            atom.meandensity = np.mean(voxelsofatom)
            atom.mediandensity = np.median(voxelsofatom)
            atom.mindensity = min(voxelsofatom)
            atom.maxdensity = max(voxelsofatom)
            # also determine number of voxels assigned per atom
            atom.numvoxels = len(voxelsofatom)
    
    # check that the last step has worked
    # (may not be necessary to have but good to check when testing!)
    for atom in PDBarray:
        for vxl in vxlsperatom[atom.atomnum-1].vxls:
            if atom.atomnum != vxl.atmnum:
                print 'error!'
     
    print '---> success'
    sys.stdout.flush()

    # delete the vxlsperatom list now to save memory
    del vxlsperatom
      
    print 'Plotting scatter plots for electron density statistics...'
    #plot a scatter plot of mean vs max
    var = ['mean','max']
    edens_scatter(where,var,PDBarray,pdbname)
    #plot a scatter plot of mean vs median
    var = ['mean','median']
    edens_scatter(where,var,PDBarray,pdbname)
    #plot a scatter plot of mean vs min
    var = ['mean','min']
    edens_scatter(where,var,PDBarray,pdbname)
    #plot a scatter plot of min vs max
    var = ['min','max']
    edens_scatter(where,var,PDBarray,pdbname)
    
    #perform residue analysis for datatset, outputting boxplots for each atom specific
    #to each residue, and also a combined boxplot across all residues in structures
    toplot = 'y'
    densmet = 'mean'
    residueArray = densper_resatom_NOresidueclass(where,PDBarray,toplot,densmet,pdbname)
    
    minresnum = 0
    sideormain = ['sidechain','mainchain']
    densper_res(where,residueArray,minresnum,sideormain,'mean',pdbname)

    # remove residueArray now to save memory 
    residueArray = []

    return PDBarray

