# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 22:49:37 2014
@author: charlie
"""
from format_by_residue import createresiduelist
import sys
from matplotlib import pyplot as plt 
import numpy as np
from boxplotter import residue_boxplotter,residue_violinplotter
import plotly.plotly as py
from plotly.graph_objs import *


def densper_resatom(PDBarray,toplot,densmet):
    # 'PDBarray' is list of atoms of structure
    # 'toplot' takes values 'y' and 'n' to determine whether individual 
    # boxplots plotted
    # 'densmet' takes values 'mean','median','max','min' to determine 
    # metric of electron density
    
    # first want to read in a list of objects grouped by residue type
    residueArray = createresiduelist(PDBarray)
    
    # next, for each specific residue type, want to calculate average 
    # distribution in electron density over atoms 
    for index in range(0,24):
        if residueArray[index].quantity > 0:
            print residueArray[index].name\
            + ' --------------------------------------------> '\
            + str(residueArray[index].quantity)
            
            densitylist2 = []
            atomtypelist2 = []
            
            for atoms in residueArray[index].res_list:
                print '\tChain: ' + str(atoms[0].chaintype)\
                + '; resi num: ' + str(atoms[0].residuenum)\
                + '; num atoms: ' + str(len(atoms)) 
                
                densitylist1 = []
                atomtypelist1 = []
                for atom in atoms:
                    
                    if densmet in ('mean'):
                        metric = atom.meandensity
                    elif densmet in ('median'):
                        metric = atom.mediandensity
                    elif densmet in ('min'):
                        metric = atom.mindensity
                    elif densmet in ('max'):
                        metric = atom.maxdensity
                    else:
                        print 'Unknown electron density metric specified'
                        sys.exit()
                    
                    print '\t' + str(atom.atomtype) +'\t'\
                    +str(atom.atomnum)+ ' ---------> ' + str(densmet)\
                    + ' density: ' + str.format('{0:.10f}',metric)
                    
                    densitylist1.append(metric)
                    atomtypelist1.append(atom.atomtype)
                    
                densitylist2.append(densitylist1)
                atomtypelist2.append(atomtypelist1)
            
            # determine the size of the full amino/nucleic acid -> do 
            # not want to include partial residues in distribution of 
            # electron density per atom type for each residue to remove 
            # potential bias towards more frequently represented atom types
            newdensitylist = []        
            completeresiduesize = max(len(densitylist2[i]) for i in range(0,len(densitylist2)))
            counter = 0        
            for residue in densitylist2:
                if len(residue) == completeresiduesize:
                    counter += 1
                    newdensitylist.append(residue)
                    
            print '------------SUMMARY---------------'
            print 'Total number of residues of type '\
            + residueArray[index].name + ': ' + str(len(densitylist2))
            print 'Maximum size of residue calculated to be: '\
            + str(completeresiduesize)      
            print 'Residues kept for use: ' + str(counter)
            print 'Residues removed since not complete: '\
            + str(len(densitylist2)-counter)
           
            # next need to find transpose of array of residue number 
            # against atom type, in order to collect atoms of each type 
            # and their corresponding density for the given residue            
            densitylist_atomorder = zip(*newdensitylist)
            
            # additionally need to find identities of atoms in order for
            # complete residues of given type        
            for residue in atomtypelist2:
                if len(residue) == completeresiduesize:
                    atoms_in_residue = residue
                    break
            
            # call the boxplotting function to plot boxplot of atom type
            # against electron density change for given residue type
            if toplot in ('y','Y','YES','Yes','yes'):
                print 'Plotting now...'
                residue_violinplotter(densitylist_atomorder,
                                      atoms_in_residue,
                                      residueArray[index].name,
                                      residueArray[index].quantity)
    
    return residueArray
    
    
    
from res_formatter import res2atomsbytype

def densper_resatom_NOresidueclass(PDBarray,toplot,densmet):
    # 'PDBarray' is list of atoms of structure
    # 'toplot' takes values 'y' and 'n' to determine whether individual 
    # boxplots plotted
    # 'densmet' takes values 'mean','median','max','min' to determine 
    # metric of electron density
        
    # first ensure PDB list ordered by atom number
    PDBarray.sort(key=lambda x: (x.basetype,x.atomnum))
    
    aminoacids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                  'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                  'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']                  
    DNAbases = ['DA','DC','DG','DT']   
    
    restype_list  = aminoacids + DNAbases
    

    residueArray = []
    for res in restype_list:
        #create an object for the residue type (in residuetype class)
        res_obj = res2atomsbytype(PDBarray,res)
                 
        print '------------SUMMARY---------------'
        print 'Total number of residues of type %s: %s'\
        %(res_obj.name,str(res_obj.frequency))
        print 'Maximum size of residue calculated to be: %s'\
        %str(res_obj.atm_names) 
      
        # need to merge 2d list res_obj of atoms by atom type into 1d list
        merged_atms = sum(res_obj.atms_bytype,[])
        
        # determine the density metric to use and create list of densities
        # from merged_atms
        if densmet in ('mean'):
            merged_dens = [atom.meandensity for atom in merged_atms]
        elif densmet in ('median'):
            merged_dens = [atom.mediandensity for atom in merged_atms]
        elif densmet in ('min'):
            merged_dens = [atom.mindensity for atom in merged_atms]
        elif densmet in ('max'):
            merged_dens = [atom.maxdensity for atom in merged_atms]
        else:
            print 'Unknown electron density metric specified'
            sys.exit()         
          
      
        # call the boxplotting function to plot boxplot of atom type
        # against electron density change for given residue type
        if toplot in ('y','Y','YES','Yes','yes'):
            print 'Plotting now...'
            residue_violinplotter(merged_dens,
                                  res_obj.atm_names,
                                  res_obj.name,
                                  res_obj.frequency)
        
        residueArray.append(res_obj)
    return residueArray
    
    
    
def densper_res(residueArray,minresnum,sideormain):  
    # plots a boxplot for each residue detailing the electron density 
    # distn. 
    # 'minresnum' is the threshold for the min number of residues of a 
    # given type that need to be present in structure to be included in 
    # the end plot
    # 'sideormain' specifies whether 'sidechain' only, 'mainchain' only, 
    # or ['sidechain','mainchain'] are selected
    residue_densities = []
    residue_label = []

    # The next section of code plots a boxplot for each residue/base 
    # type in structure    
    for index in range(0,24):
        if residueArray[index].quantity > minresnum:
            
            print str(residueArray[index].name) + ' -----> '\
            + str(residueArray[index].quantity)

            densitylist2 = []
            for atoms in residueArray[index].res_list:
                densitylist1 = []
                for atom in atoms:
                    if atom.side_or_main() in (str(sideormain)):
                        
                        densitymetric = (1+atom.mindensity)*(1+atom.mediandensity)
                        #densitymetric2 = (1+atom.meandensity)*(atom.maxdensity-atom.mindensity)
                        #densitymetric3 = atom.mindensity
                        #densitymetric4 = atom.mediandensity
                        
                        densitylist1.append(densitymetric)
                densitylist2.append(densitylist1)
            
            # convert list of atom density lists for each residue 
            # representative to 1d merged list for further use
            mergedlist = []
            for element in densitylist2:
                mergedlist = mergedlist + element
            residue_densities.append(mergedlist)
            residue_label.append(residueArray[index].name)
            
    print 'Plotting all present residues now......'
    residue_violinplotter(residue_densities,
                          residue_label,
                          'All residues',
                          str(sideormain))
    



def densper_sidemain(residueArray,minresnum):
    # the following plots a plotly boxplot of sidechain and mainchain 
    # density changes for each residue type individually to identify 
    # specific damage heterogeneity between residue types (in the form 
    # of significant side chain electron density disorder).
    # 'minresnum' is the threshold for the min number of residues of a 
    # given type that need to be present in structure to be included in 
    # the end plot

    residue_densities_main = []
    residue_densities_side = []
    residue_label = []
       
    # The next section of code plots a boxplot for each residue/base 
    # type in structure    
    for index in range(0,24):
        if residueArray[index].quantity > minresnum:
            print str(residueArray[index].name) + ' -----> '\
            + str(residueArray[index].quantity)
            
            densitylist2_side = []
            densitylist2_main = []
            for atoms in residueArray[index].res_list:
                densitylist1_side = []
                densitylist1_main = []
                for atom in atoms:
                    densitymetric = (1+atom.mindensity)*(1+atom.mediandensity)   
                    
                    if atom.side_or_main() in ('sidechain'):
                        densitylist1_side.append(densitymetric)
                    elif atom.side_or_main() in ('mainchain'):
                        densitylist1_main.append(densitymetric)

                densitylist2_side.append(densitylist1_side)
                densitylist2_main.append(densitylist1_main)
                            
            # convert list of atom density lists for each residue 
            # representative to 1d merged list for further use
            mergedlist_main = []
            for element in densitylist2_main:
                mergedlist_main = mergedlist_main + element
            residue_densities_main.append(mergedlist_main)

            mergedlist_side = []
            for element in densitylist2_side:
                mergedlist_side = mergedlist_side + element
            residue_densities_side.append(mergedlist_side)
            
            residue_label.append(residueArray[index].name)
            
    py.sign_in("charlieb333", "11090mkre9")

    x_side = []
    y_side = []
    counter = 0
    for res in residue_densities_side:
        y_side = y_side + res
        x_side = x_side + [str(residue_label[counter])]*len(res)
        counter += 1
        
    x_main = []
    y_main = []
    counter = 0
    for res in residue_densities_main:
        y_main = y_main + res
        x_main = x_main + [str(residue_label[counter])]*len(res)
        counter += 1
        
    trace1 = Box(y=y_main,x=x_main,name='main',marker=Marker(color='#3D9970'))
    trace2 = Box(y=y_side,x=x_side,name='side',marker=Marker(color='#FF4136'))

    data = Data([trace1, trace2])
    
    layout = Layout(yaxis=YAxis(title='density loss',zeroline=True),boxmode='group')
    fig = Figure(data=data, layout=layout)
    plot_url = py.plot(fig, filename='box-grouped')
                    
            
            
            
    
    
         
            