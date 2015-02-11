# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 20:35:38 2015

@author: charlie
"""
import emailme
from assign_map2atoms_v2 import densperatom
from densdist_perresidue import densper_resatom,densper_res
from densityanalysisplots import edens_scatter
from bdam import bdam_calculate
from map2voxelclasslist import densmap2class

def master(dataset_info):
    ###############################################################################
    print 'Converting density map...'
    dataset_info.voxel_list = densmap2class(dataset_info.mapname)
    ###############################################################################
    
    
    ###############################################################################
    #work out start time for script (to be used at end to determine duration of script)
    #starttime = emailme.scriptstart()

    print '\nCalculating density per atom...'
    PDBarray_full = densperatom(dataset_info)

    #upon successful end of script, email csbury@me.com to confirm job success
    #emailme.scriptdone(starttime)
    ###############################################################################
    
    
    ###############################################################################
    #Analysis of summary_edensity.txt file generated after running 
    #'findaverageelectrondensityaroundPDBatoms.py' script
    ###############################################################################
    #plot a scatter plot of mean vs max
    var = ['mean','max']
    edens_scatter(var,PDBarray_full,dataset_info.pdbname)
    #plot a scatter plot of mean vs median
    var = ['mean','median']
    edens_scatter(var,PDBarray_full,dataset_info.pdbname)
    #plot a scatter plot of mean vs min
    var = ['mean','min']
    edens_scatter(var,PDBarray_full,dataset_info.pdbname)
    #plot a scatter plot of min vs max
    var = ['min','max']
    edens_scatter(var,PDBarray_full,dataset_info.pdbname)
    ###############################################################################
    

    ###############################################################################
    #perform residue analysis for datatset, outputting boxplots for each atom specific
    #to each residue, and also a combined boxplot across all residues in structures
    toplot = 'y'
    densmet = 'mean'
    residueArray = densper_resatom(PDBarray_full,toplot,densmet)
    
    minresnum = 0
    sideormain = 'sidechain'
    densper_res(residueArray,minresnum,sideormain)
    ###############################################################################
    
       
#    ###############################################################################
#    #calculate Bdamage style metric for current pdb structure
#    threshold = 14
#    bdam_calculate(PDBarray_full,threshold)
#    ###############################################################################
    
    return PDBarray_full