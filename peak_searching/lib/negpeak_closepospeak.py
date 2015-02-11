# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 16:22:58 2014

@author: charlie
"""
from math import sqrt, pow

######################################################################################################################################################
######################################################################################################################################################
#SECTION: FIND DIFFERENCE MAP PEAKS WITH CLOSE POSITIVE PEAKS TO FLAG UP FOR MANUAL SEARCHING
######################################################################################################################################################

######################################################################################################################################################
#Update:25nov2014: if negative peak is associated with a side chain atom of 
#structure, it is suggested that the nearest positive peak to this could be
#a conformational change (different rotamer) for the side chain - as such it is 
#now flagged up appropriately to the command line
##########################################################################################################################################################

def negpeak_closepospeak(filelist,datedir,PDBarray):

    #read in the negative damage peaks .pdb file into script
    pdbdamagefile_original = filelist[1] 
    pdbdamagefile = str("negpeaksonly_thresholddistance"+datedir.date+".pdb")
    pdbdamagein = open(datedir.path + '/' + str(pdbdamagefile), "r" )    
    pdbdamagelines = pdbdamagein.readlines()
    
    minpospeakdist = []
    
    print 'Searching for potential conformation side chain changes...'
    print 'For negative peaks with positive peaks closer than preset distance (2.5 Angstrom chosen here), damage locations are flagged'
    sidechainlist = []
    
    for line in pdbdamagelines:
        if 'ATOM' not in line[0:5]:
            pass
        else:           
            #for a given line in the .pdb difference map file assign x,y,z 
            #cartesian coordinates for the difference map peak
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            ######################################################################
            #update:25nov2014
            proDNAcoords = open(datedir.path+'/'+"pdbcoords.txt", "r")
            proDNAlines = proDNAcoords.readlines()
            Dvector = []
            #read through the structure pdb file and calculate the closest 
            #atom location to the difference map peak site using a basic
            #L2-norm calculation
            for prodnaline in proDNAlines:
                xi = float(prodnaline[30:38].strip())
                yi = float(prodnaline[38:46].strip())
                zi = float(prodnaline[46:54].strip())
                Di = sqrt( pow(x-xi,2) + pow(y-yi,2) + pow(z-zi,2))
                Dvector.append(Di)
            mindist = min(Dvector)
            minIndex = Dvector.index(min(Dvector))        
            proDNAcoords.close()
            
            if str(PDBarray[minIndex][17:20].strip()) in ('DT','DA','DC','DG'):
                #print str(PDBarray[minIndex].split()[3]) + ' is DNA'
                if str(PDBarray[minIndex][12:16].strip()) not in ("P","OP1","O5'","C5'","C4'","C3'","O3'","C2'","C1'","O4'","OP2"):
                    #print 'base of DNA detected'
                    sidechainindicator = 'baseside'
                else:
                    sidechainindicator = 'no'
                    pass
            else:
                #print 'is protein'
                if str(PDBarray[minIndex][12:16].strip()) not in ('N','CA','C','O'):
                    #print 'is a side chain atom of residue'
                    sidechainindicator = 'proteinside'

                else:
                    sidechainindicator = 'no'
                    pass
            ######################################################################
            
            ######################################################################
            pdbdamagein_withpospeaks = open(str(pdbdamagefile_original), "r" )    
            pdbdamagelines_withpospeaks = pdbdamagein_withpospeaks.readlines()
            Dvector = []
            #read through the original damage peak pdb file (containing positive peaks too) and calculate the closest 
            #positive peak location to the difference map peak site using a basic
            #L2-norm calculation
            for pdbdamageline_withpospeaks in pdbdamagelines_withpospeaks:
                if pdbdamageline_withpospeaks.split()[0] == 'ATOM': 
                    if float(pdbdamageline_withpospeaks[54:60].strip()) > 0:
                        xi = float(pdbdamageline_withpospeaks[30:38].strip())
                        yi = float(pdbdamageline_withpospeaks[38:46].strip())
                        zi = float(pdbdamageline_withpospeaks[46:54].strip())
                        Di = sqrt( pow(x-xi,2) + pow(y-yi,2) + pow(z-zi,2))

                    elif float(pdbdamageline_withpospeaks[54:60].strip()) < 0:
                        Di = float('Inf')
                    Dvector.append(Di)
                else:
                    pass
            ######################################################################
                
            mindist = min(Dvector)
            pdbdamagein_withpospeaks.close()           
            minpospeakdist.append(mindist)
            
            ######################################################################
            #update: 25nov2014: if nearest negative peak is side chain and closest 
            #positive peak is close by (at a threshold distance specified here)
            sidechainlist.append(sidechainindicator)
            if sidechainindicator in ('proteinside','baseside') and mindist < 2.5:
                print 'Potential conformational change detected: '+ str(PDBarray[minIndex][12:16]).strip() +\
                      ' '+ str(PDBarray[minIndex][17:20]).strip() + ' with closest positive peak at: ' +\
                      str(mindist) + ' Angstrom'
                      
            
    pdbdamagein.close()
    return minpospeakdist,sidechainlist
