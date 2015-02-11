# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 00:17:18 2014

@author: charlie
"""
from math import sqrt, pow
from file_len_function import file_len
import sys

##############################################################################################################################################################
#SECTION: REMOVE PEAKS FURTHER THAN A SPECIFIED DISTANCE THRESHOLD FROM ATOMS IN THE STRUCTURE PDB FILE
##############################################################################################################################################################

def threshold_peak_remover_VDVradius(threshold,datedir,PDBarray):    
    #give option to skip distance threshold stage - if omitted, later functions 
    #will not input correct difference map files. This question allows, easy exiting
    #from script if desired.
    print 'Applying a distance threshold to difference map peaks from protein/nucleic acid structure, as specified in input.txt'

    inclusionlist = []
    
    negpdb = open(datedir.path+'/'+"negpeaksonly"+datedir.date+".pdb", "r")
    negha = open(datedir.path+'/'+"negpeaksonly"+datedir.date+".ha", "r")

    pdblines = negpdb.readlines()
    halines = negha.readlines()

    #write new .ha and .pdb file containing negative peaks only, below the given 
    #threshold from the pdb structure being investigated
    negpdbdistthreshold = open(datedir.path+'/'+"negpeaksonly_thresholddistance"+datedir.date+".pdb", "w")
    neghadistthreshold = open(datedir.path+'/'+"negpeaksonly_thresholddistance"+datedir.date+".ha", "w")
    
    counter = 0;
    
    for line in pdblines:
        #write preamble lines to .pdb files so files can be read in COOT/PYMOL
        if 'ATOM' not in line[0:5]:
            negpdbdistthreshold.write(line)
        else:           
            counter = counter + 1           
            #for a given line in the .pdb difference map file assign x,y,z 
            #cartesian coordinates for the difference map peak
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            proDNAcoords = open(datedir.path+'/'+"pdbcoords.txt", "r")
            proDNAlines = proDNAcoords.readlines()
            Dvector = []
            atomnum_list = []

            #read through the structure pdb file and calculate the closest 
            #atom location to the difference map peak site using a basic
            #L2-norm calculation
            for prodnaline in proDNAlines:
                xi = float(prodnaline[30:38].strip())
                yi = float(prodnaline[38:46].strip())
                zi = float(prodnaline[46:54].strip())
                Di = sqrt( pow(x-xi,2) + pow(y-yi,2) + pow(z-zi,2))
                Dvector.append(Di)
                atomnum_list.append(int(prodnaline[6:11].strip()))

            mindist = min(Dvector)
            minIndex = Dvector.index(min(Dvector))        
            proDNAcoords.close()
            
            ################################################################################
            #24NOV2014 UPDATE OF 'threshold_peak_remover': FIND THE ATOM TYPE 
            #IDENTIY OF THE NEAREST ATOM TO PEAK, AND SET THRESHOLD CORRESPONDINGLY
            # 2FEB2015: UPDATE --> just read atom-identifier segment of PDB atom line to get
            # correct atom code
            nearestatomtype = PDBarray[minIndex][12:16].strip()
            nearestresidue = PDBarray[minIndex][17:20].strip()
            nearestchain = PDBarray[minIndex][21].strip()
            atomidentifier= PDBarray[minIndex][76:78].strip()
            
            #read through VDV file (Van der Waals radius list for elements) and pick correct element 
            VDVradiusinfo = open("./inputfiles/VDVradiusfile.txt", "r")
            VDVlines = VDVradiusinfo.readlines()
            indicator = 0

            for VDVline in VDVlines:
                if atomidentifier == VDVline.split()[1]:
                    peak_threshold = VDVline.split()[5]
                    indicator +=1
                    break

            if indicator == 0:
                sys.exit("Nearest atom type to peak not recognised")            
            ################################################################################
            
            #keep any difference map peaks closer than the minimum threshold
            #0.01 factor included since VDV radii given in pm not Angstrom in inputted txt file
            if mindist <= float(peak_threshold)*(0.01):
                print str(counter) + ': Nearest: '+atomidentifier+' '+nearestresidue+' '+ nearestchain +\
                      ' @ ' + str(mindist)[0:4] + ' Angstrom ---> peak kept!'
                negpdbdistthreshold.write(line)
                inclusionlist.append(counter)
            else:
                print str(counter) + ': Nearest: '+atomidentifier+' '+nearestresidue+' '+ nearestchain +\
                      ' @ ' + str(mindist)[0:4] + ' Angstrom ---> peak removed!'         
    #this code writes the new .ha difference peak file accordingly using the above information
    counter2 = 0;
    for line in halines[0:2]:
        neghadistthreshold.write(line)
    for line in halines[2:file_len(datedir.path+'/'+"negpeaksonly"+datedir.date+".ha")]:
        counter2 = counter2 + 1
        if counter2 in inclusionlist:
            neghadistthreshold.write(line)

    negpdb.close()
    negha.close()
    negpdbdistthreshold.close()
    neghadistthreshold.close()
    
    #helpful commandline statement of how many peaks have been kept in file
    #(gives an initial idea of the extent of specific damage, and whether the
    #chosen distance threshold value was suitable)
    num_kept_peaks = file_len(datedir.path+'/'+"negpeaksonly_thresholddistance"+datedir.date+".ha") - 2
    print 'Number of peaks kept in file is', num_kept_peaks
    return num_kept_peaks
###############################################################################