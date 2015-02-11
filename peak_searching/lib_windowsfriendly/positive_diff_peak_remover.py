# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 23:35:57 2014

@author: charlie
"""
import sys

######################################################################################################################################################
######################################################################################################################################################
#SECTION: REMOVE POSITIVE DIFFERENCE MAP PEAKS
######################################################################################################################################################
######################################################################################################################################################

def positive_diff_peak_remover(filelist,datedir):

    pdbfilename = filelist[1]
    hafilename = filelist[2]
    ################################################################################
    #input file names for .pdb and .ha files (outputted by FFT) from which the
    #positive peaks should be removed automatically
    
    pdbin = open(str(pdbfilename), "r" )
    hain = open(str(hafilename), "r" )
    
    pdblines = pdbin.readlines()
    halines = hain.readlines()
    
    #writes a new .pdb and .ha file containing ONLY negative difference map peaks
    #and places it in the directory created above
    newfilepdb = open(datedir.path +'/'+"negpeaksonly"+datedir.date+".pdb", "w")
    newfileha = open(datedir.path +'/'+"negpeaksonly"+datedir.date+".ha", "w")
    
    #reads through the .pdb and .ha files and removes positive peaks, whilst
    #keeping the necessary preamble at the top of each file type to allow them to
    #be read
    counter_pdb = 0
    for line in pdblines:
        if 'ATOM' not in line[0:5]:
            newfilepdb.write(line)
        elif 'ATOM' in line[0:5]:
            if float(line[54:60]) <= 0:
                newfilepdb.write(line)
                counter_pdb += 1
            else:
                pass

    print '#'+ ' '*5 + 'Type' + ' '*10 + 'Sigma' + ' '*11 + 'Status'

    counter = 0
    line_counter = 0
    for line in halines:
        if 'ATOM' not in line[0:4]:
            newfileha.write(line)
        elif 'ATOM' in line[0:4]:
            line_counter += 1
            if float(line.split()[5]) <= 0:
                counter+=1
                print str(line_counter) + ': Negative   @   ' + str(float(line.split()[5])) + '   ---->   kept'
                newfileha.write(line)
            else:
                print str(line_counter) + ':  Positive   @    ' + str(float(line.split()[5])) + '   ---->   removed'

    # as a check, determine whether diff peak pdb and ha files contain same number of neg peaks
    if counter != counter_pdb:
        print 'Incompatible .pdb and .ha files.. different numbers of negative peaks detected'
        print '---> Terminating script...'
        sys.exit()
                
    print 'Breakdown summary: ' + str(counter) + ' negative peaks found in current file.'   
    pdbin.close()
    hain.close()
    newfilepdb.close()
    newfileha.close()
    return counter
    ###############################################################################