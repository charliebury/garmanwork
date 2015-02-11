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
from densandatommaps2vxllist import combinevxl_atmanddensvals

###----------------###----------------###----------------###----------------###
###############################################################################      
# A class for .map file voxel
class voxel_density:
    def __init__(self,density=0,atmnum=0):
        self.density = density
        self.atmnum = atmnum
###############################################################################
###----------------###----------------###----------------###----------------###


def maps2atmdensity(pdbname,mapfilname1,maptype1,mapfilname2,maptype2):

    print '------------------------------------------------'
    print 'Reading in pdb file...'   
    # next read in the pdb structure file:
    # run function to fill PDBarray list with atom objects from structure
    PDBarray = pdb2list(pdbname+'.pdb',[])
    print '---> success'
    
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
    atmmap = densmap2class(mapfilname1,maptype1)
    
    # find atom numbers present in list (repeated atom numbers removed)
    seen = set()
    seen_add = seen.add
    uniq_atms = [x for x in atmmap.vxls_val if not (x in seen or seen_add(x))] 
    
    # find set of atoms numbers not present (i.e atoms not assigned to voxels)
    Atms_notpres = set(range(1,num_atoms+1)) - set(uniq_atms)
    print 'Number of atoms not assigned to voxels: %s' %str(len(Atms_notpres))
    
    # read in the density map
    print '------------------------------------------------'
    print 'Density map read in: '
    densmap = densmap2class(mapfilname2,maptype2)  
    
    print '------------------------------------------------'
    print 'Checking that maps have same number of voxels...'
    # check that the maps have the same length 1D list of voxels
    if len(atmmap.vxls_val) != len(densmap.vxls_val):
        print 'Two maps have different number of voxels --> cannot handle'
        print 'Terminating script...'
        sys.exit()
    print '---> success'
    
    
    print '------------------------------------------------'
    print 'Checking that maps have same dimensions and sampling properties...' 
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
    

    # remove voxels not assigned to atoms
    print '------------------------------------------------'      
    print 'Total number of voxels: %s' %str(len(atmmap.vxls_val)) 
    # find list indices corresponding to non-atom voxels (zero values)
    
    nonatm_indices = [i for i, j in zip(count(),atmmap.vxls_val) if j == 0]    
    print 'Number of voxels not assigned to atoms: %s' %str(len(nonatm_indices)) 
   
    print '---> removing these voxels'
    filteredlist_atms = multi_delete(atmmap.vxls_val,nonatm_indices)
    filteredlist_dens = multi_delete(densmap.vxls_val,nonatm_indices)
    print '    ---> success!'

    # delete atmmap and densmap now to save memory
    densmap =[]
    atmmap = []
    
    
    # create list of voxel objects in class voxel_density 
    print '------------------------------------------------'
    print 'Combining voxel density and atom values...'
    vxl_list = combinevxl_atmanddensvals(filteredlist_atms,filteredlist_dens)
    
    # next assign voxels to atom objects in PDBarray list
    print '------------------------------------------------'
    print 'Assigning voxels to corresponding atoms...'
    # first order voxels by atom number
    vxl_list.sort(key=lambda x: x.atmnum)
    
    #initialise list of voxels for each atom before can append voxels
    for atom in PDBarray:
        atom.mapvoxels = []
    
    # unessential loading bar add-in (part1)
    total = len(vxl_list)
    point = total / 100
    increment = total / 100
    
    counter = -1  
    for vxl in vxl_list:
        counter += 1
        (PDBarray[vxl.atmnum-1].mapvoxels).append(vxl)
    
    # can remove vxl_list_filt variable now to save memory
    vxl_list = []    
    
    # unessential loading bar add-in (part2)
    if(counter % (1 * point) == 0):
        sys.stdout.write(
                "\r[" + "=" * (counter / increment) +
                " " * ((total - counter)/ increment) + 
                "]" +  str(counter / point) + "%")
        sys.stdout.flush()
    
    # check that the last step has worked
    # (may not be necessary to have but good to check when testing!)
    for atom in PDBarray:
        for vxl in atom.mapvoxels:
            if atom.atomnum != vxl.atmnum:
                print 'error!'
      
    print '---> success'   
    
    # histogram plot of number of voxels per atom
    print '------------------------------------------------'
    print 'Plotting histogram of number of voxels per atom...' 
    sns.set_palette("deep", desat=.6)
    sns.set_context(rc={"figure.figsize": (10, 6)})
    fig = plt.figure()
    datax = [len(atom.mapvoxels) for atom in PDBarray]
    plt.hist(datax, 300, histtype="stepfilled", alpha=.7);
    plt.xlabel('Voxels per atom')
    plt.ylabel('Frequency')
    plt.title('Histrogram of voxels per atom')
    
    
    # determine min,max,mean,median density per atom
    print '------------------------------------------------'
    print 'Calculating electron density statistics per atom...'
    for atom in PDBarray:
        if len(atom.mapvoxels) != 0:
            atom.meandensity = np.mean([vxl.density for vxl in atom.mapvoxels])
            atom.mediandensity = np.median([vxl.density for vxl in atom.mapvoxels])
            atom.mindensity = min([vxl.density for vxl in atom.mapvoxels])
            atom.maxdensity = np.max([vxl.density for vxl in atom.mapvoxels])
    
    print 'Plotting scatter plots for electron density statistics...'
    #plot a scatter plot of mean vs max
    var = ['mean','max']
    edens_scatter(var,PDBarray,pdbname)
    #plot a scatter plot of mean vs median
    var = ['mean','median']
    edens_scatter(var,PDBarray,pdbname)
    #plot a scatter plot of mean vs min
    var = ['mean','min']
    edens_scatter(var,PDBarray,pdbname)
    #plot a scatter plot of min vs max
    var = ['min','max']
    edens_scatter(var,PDBarray,pdbname)
    
    #perform residue analysis for datatset, outputting boxplots for each atom specific
    #to each residue, and also a combined boxplot across all residues in structures
    toplot = 'y'
    densmet = 'min'
    residueArray = densper_resatom_NOresidueclass(PDBarray,toplot,densmet)
    
    minresnum = 0
    sideormain = ['sidechain','mainchain']
    densper_res(residueArray,minresnum,sideormain)

    # remove residueArray now to save memory 
    residueArray = []

    return PDBarray

