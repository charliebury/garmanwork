# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 15:44:34 2014

@author: lina2532
"""
from matplotlib import pyplot as plt    
import sys

###############################################################################
#following scatterplot generating code allows a check between the correlation 
#between the mean and median calculated electron density for each atom in structure
def electrondensity_scatter(var):
    
    filenamein = 'summary_edensity.txt'
    edensityfile = open(filenamein,'r')
        
    var1_list = []
    var2_list = []
    
    #first determine column numbers for var1 and var2 variables (such that var = [var1,var2]):
    col = []
    
    for i in range(0,2):
        if var[i] in ('mean','Mean'):
            col.append(1)
        elif var[i] in ('median','Median'):
            col.append(2)
        elif var[i] in ('min','Min'):
            col.append(3)
        elif var[i] in ('max','Max'):
            col.append(4)
        else:
            print 'Unrecognised variable name, stopping program!'
            sys.exit()

    for line in edensityfile.readlines():
        if line[0] in ('#','.','A'):
            continue
        else:
            var1_list.append(float(line.split()[col[0]]))
            var2_list.append(float(line.split()[col[1]]))
        
        #check two generated lists of same length:
        if len(var1_list) != len(var2_list):
            print 'Error: lists not same length for scatter plot'
            sys.exit()   
            
    edensityfile.close()

    scatter1 = plt.figure()
    plt.scatter(var1_list,var2_list)
    scatter1.suptitle(str(var[0]) + ' vs ' + str(var[1]) + ' electron density',fontsize=20)
    plt.xlabel(str(var[0]) + ' density',fontsize=18)
    plt.ylabel(str(var[1]) + ' density',fontsize=16)

    figname = str(var[0])+'_vs_'+str(var[1])+'.png'
    scatter1.savefig(figname)
###############################################################################
    
    
    
###############################################################################
# following scatterplot generating code is updated version (2Jan2015) of plotting 
# function above - using the PDBarray list of atom objects
def edens_scatter(where,var,PDBarray,pdbname):
    # var of form ['mean','median'] to choose two metrics of electron density per 
    # atom to plot against each other in a scatter plot
    
    valsperparam = []
    for param in var:   
        if param in ('mean','Mean'):
            valperatom = [atom.meandensity for atom in PDBarray]
        elif param in ('median','Median'):
            valperatom = [atom.mediandensity for atom in PDBarray]
        elif param in ('min','Min'):
            valperatom = [atom.mindensity for atom in PDBarray]
        elif param in ('max','Max'):
            valperatom = [atom.maxdensity for atom in PDBarray]
        else:
            print 'Unrecognised variable name, stopping program!'
            sys.exit()
        valsperparam.append(valperatom)
        
    #check two generated lists of same length:
    if len(valsperparam[0]) != len(valsperparam[1]):
        print 'Error: lists not same length for scatter plot'
        sys.exit()   
            
    scatter1 = plt.figure()
    plt.scatter(valsperparam[0],valsperparam[1])
    
    scatter1.suptitle(
        str(var[0]) + ' vs ' + str(var[1]) + ' electron density',
        fontsize=20)
        
    plt.xlabel(
        str(var[0]) + ' density',
        fontsize=18)
        
    plt.ylabel(
        str(var[1]) + ' density',
        fontsize=16)

    figname = str(where)+'output/plots/' + str(pdbname) + '_' + str(var[0])+'_vs_'+str(var[1])+'.jpg'
    scatter1.savefig(figname)
###############################################################################





###############################################################################
def Npercentmostdamaged(var,PDBarray,percent):
    densityfile = open('atom_edensity.txt','r')
    counter = 0
    num_atoms = len(PDBarray)
    res_list = []
    for line in densityfile.readlines():
        if str(var) == line.split()[0]:
            counter += 1
            if counter > num_atoms*(float(100-percent)/100):
                res_list.append(str(line.split()[4])[0:len(line.split()[4])-1])
    densityfile.close()
    
    aminoacids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    DNAbases = ['DA','DC','DG','DT']   
    restype_list  = aminoacids + DNAbases
    
    res_list_counter = [0]*len(restype_list)
    counter = 0
    for res in res_list:
        index = restype_list.index(res)
        newlist = []
        for i in range(0,index):
            newlist.append(res_list_counter[i])
        newlist.append(res_list_counter[index]+1)
        for i in range(index+1,len(restype_list)):
            newlist.append(res_list_counter[i])
        res_list_counter = newlist

    #want fractions of 100 for pie chart: 
    res_list_counter_normalised = [(float(x)*float(100)) / float(sum(res_list_counter)) for x in res_list_counter]
    explode = [0.2]*len(restype_list)
    pie1 = plt.figure()
    plt.pie(res_list_counter_normalised, explode = explode, labels=restype_list,
                autopct='%1.1f%%', shadow=True, startangle=90)


    pie1.suptitle(str(percent)+'%-tile ' +str(var) + ' density' ,fontsize=20)
    figname = str(percent)+'%' +str(var) + ' density_piechart'+'.png'
    pie1.savefig(figname)
###############################################################################
    


###############################################################################
def bfac_scatter(pdbfilename,PDBarray,mainorside,restypes,atomtypes,densmet,bmet):
    ##simply plots scatter plot of two numerical attributes of atom objects of StructurePDB class
    #'PDBarray' is list of atoms
    #'mainorside' specifies 'mainchain', 'sidechain'
    #'restypes' specifies 'GLU','MET','HIS','DA',... etc
    #'atomtypes' specifies atom types eg: 'C','O','S',...etc
    #'densmet' specifies electron density metric per atom ('mean','min',...etc)
    #'bmet' specifies the Bfactor associated metric per atom ('Bfactor','Bdamage',...etc)
    
    X = []
    Y = []
    
    for atom in PDBarray:
        if atom.side_or_main() in (mainorside) and atom.basetype in (restypes) and atom.atomidentifier in (atomtypes):
            if densmet in ('mean','Mean'):
                X.append(atom.meandensity) 
            elif densmet in ('median','Median'):
                X.append(atom.mediandensity) 
            elif densmet in ('min','Min'):
                X.append(atom.mindensity) 
            elif densmet in ('max','Max'):
                X.append(atom.maxdensity) 
            else:
                print 'Unrecognised variable name, stopping program!'
                sys.exit()
            
            if bmet in ('Bfactor'):
                Y.append(atom.Bfactor) 
            elif bmet in ('Bfactorchange'):
                Y.append(atom.Bfactorchange) 
            elif bmet in ('Bdamage'):
                Y.append(atom.bdam) 
            elif bmet in ('Bdamagechange'):
                Y.append(atom.bdamchange) 
            else:
                print 'Unrecognised variable name, stopping program!'
                sys.exit()
               
    scatter1 = plt.figure()
    plt.scatter(X,Y)    
    scatter1.suptitle(str(densmet) + ' electron density' + ' vs ' + str(bmet),fontsize=20)
    plt.xlabel(str(densmet) + ' electron density',fontsize=18)
    plt.ylabel(str(bmet),fontsize=16)
    
    figname = str(pdbfilename) + '_' + str(densmet) + ' electron density' + ' vs ' + str(bmet) + '.jpg'
    scatter1.savefig(figname)
############################################################################### 