# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 11:28:50 2015

@author: charlie
"""

import sys   
from PDBfile_manipulation import PDBtoCLASSARRAY_v2 as pdb2list
from densityanalysisplots import edens_scatter
from res_formatter import densper_resatom_NOresidueclass,densper_res
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from map2voxelclasslist_nocartconversion import densmap2class_readheader,densmap2class_readvoxels,densmap2class_consistencycheck

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
    atmmap, densitystart_atom = densmap2class_readheader(mapfilname1)
    
    # read in the density map
    print '------------------------------------------------'
    print 'Density map read in: '
    densmap, densitystart_dens = densmap2class_readheader(mapfilname2)      
    
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
    
    
    print '------------------------------------------------'
    print 'Checking that maps have same number of voxels...'
    # check that both maps have the same starting position
    if densitystart_atom != densitystart_dens:
        print 'Error --> maps not compatibile with this script, since'\
        + 'non identical voxel start location'
        print '---> terminating script...'
        sys.exit()
        
    # check that both maps have same file size
    filesize_atommap = os.path.getsize(mapfilname1)
    filesize_densmap = os.path.getsize(mapfilname2)
    if filesize_atommap != filesize_densmap:
        print 'Error --> maps not compatibile with this script, since'\
        + 'non identical map file size'
        print '---> terminating script...'
        sys.exit()
    print '---> success'

    print 'Reading through maps to find voxel information...'
    # now read through data array in both map files and output a 
    vxl_list_filt = densmap2class_readvoxels(mapfilname1,mapfilname2,densitystart_atom)
    print '---> success'
    # perform consistency check with header of both files    
    densmap2class_consistencycheck(vxl_list_filt,
                                   atmmap.maxdensity,atmmap.mindensity,
                                   densmap.maxdensity,densmap.mindensity)    
    
    print 'Total number of voxels associated with atoms: %s' %str(len(vxl_list_filt)) 
    # delete atmmap and densmap now to save memory
    densmap =[]
    atmmap = []
    
    # find atom numbers present in list (repeated atom numbers removed)
    seen = set()
    seen_add = seen.add
    atm_vxls = [vxl.atmnum for vxl in vxl_list_filt]
    uniq_atms = [x for x in atm_vxls if not (x in seen or seen_add(x))] 
    
    # find set of atoms numbers not present (i.e atoms not assigned to voxels)
    Atms_notpres = set(range(1,num_atoms+1)) - set(uniq_atms)
    print 'Number of atoms not assigned to voxels: %s' %str(len(Atms_notpres))
    
#    # remove voxels not assigned to atoms
#    print '------------------------------------------------'
#    vxl_list_filt = [vxl for vxl in vxl_list if vxl.atmnum != 0]
#    print 'Total number of voxels: %s' %str(len(vxl_list)) 
#    print 'Number of voxels not assigned to atoms: %s' %str(len(vxl_list)-len(vxl_list_filt)) 
#    print '---> removing these voxels'
#    
#    # next delete variable vxl_list to save memory
#    vxl_list = []
    
    # next assign voxels to atom objects in PDBarray list
    print '------------------------------------------------'
    print 'Assigning voxels to corresponding atoms...'
    # first order voxels by atom number
    vxl_list_filt.sort(key=lambda x: x.atmnum)
    
    #initialise list of voxels for each atom before can append voxels
    for atom in PDBarray:
        atom.mapvoxels = []
    
    # unessential loading bar add-in (part1)
    total = len(vxl_list_filt)
    point = total / 100
    increment = total / 100
    
    counter = -1  
    for vxl in vxl_list_filt:
        counter += 1
        (PDBarray[vxl.atmnum-1].mapvoxels).append(vxl)
    
    # can remove vxl_list_filt variable now to save memory
    vxl_list_filt = []    
    
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

