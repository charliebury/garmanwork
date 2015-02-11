# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:01:11 2014

@author: charlie
"""
from file_len_function import file_len

##############################################################################################################################################################
#SECTION: MANUAL SEARCHING THROUGH REMAINING DIFFERENCE MAP PEAKS TO REMOVE UNWANTED PEAKS NOT CONSIDERED TO BE ASSOCIATED WITH SPECIFIC DAMAGE
##############################################################################################################################################################
##############################################################################################################################################################

def manual_peak_searching(datedir,minpospeakdist):

    ###############################################################################
    #read in filtered .ha and .pdb files generated above
    print colored('Code to select wanted peaks from .ha file (from FTT over whole .pdb file)','cyan') + colored(' and then to change these difference map peak values to the correct values as displayed by COOT (from FTT over only asymmetric unit)','magenta')
    
    hafilename = str("negpeaksonly_thresholddistance"+datedir.date+".ha")
    pdbfilename = str("negpeaksonly_thresholddistance"+datedir.date+".pdb")
    insha = open(datedir.path + '/' + str(hafilename), "r" )
    inspdb = open(datedir.path + '/' + str(pdbfilename), "r" )
    ###############################################################################
    
    ###############################################################################
    #read .ha and .pdb lines into arrays to allow easier manipulation later
    arrayha = []
    arraypdb = []
    
    for line in insha:
        arrayha.append( line )
    insha.close()
    
    for line in inspdb:
        arraypdb.append( line )
    inspdb.close()
    ###############################################################################
    
    ###############################################################################
    #write new .ha and .pdb files containing only manually selected peaks
    print 'Code to delete unwanted entries in the difference map peak .ha file.'
    newfileha = open(datedir.path + '/' + "editted"+str(hafilename), "w")
    newfilepdb = open(datedir.path + '/' + "editted"+str(pdbfilename), "w")
    
    #write preamble lines to new .ha file
    array2 = arrayha[0:2]
    for element in array2:
        newfileha.write(str(element))
    
    #write preamble lines to new .pdb file
    arraypdb2 = arraypdb[0:4]
    for element in arraypdb2:
        newfilepdb.write(str(element))
        
    #ask if manual peak searching is intended now - if no script will terminate early; 
    #useful since manual searching can be time consuming   
    print 'Do manual peak analysis now? (y/n)'
    print 'For datasets with lots of diff peaks present, this can take a long time...'
    response1 = raw_input("Please enter something: ")
    
    arrayha_coords = arrayha[2:file_len(datedir.path + '/' + str(hafilename))]
    arraypdb_coords = arraypdb[4:file_len(datedir.path + '/' + str(pdbfilename))-1]
    
    counter1 = 0
    counter2 = 0
    if response1 in ('y','Y','yes','YES','Yes'):
        
        #for each peak in .ha file, state properties (cartesian location, 
        #difference map peak value etc) to command line. NEED to use COOT here to 
        #move between difference maps in most recent .pdb difference map file above
        #('search difference map peaks' tool and click enter to move to next peak).
        #Ordering of difference map peaks in COOT is required to be same as in .pdb 
        #and .ha file - FFT CCP4 program and above filtering has been observed to 
        #keep this ordering in previously tested datasets.
        for element in arraypdb_coords:
            counter1 += 1
            print str(counter1)+': Diff peak #' + str(element[22:26]).strip() +' @ '+\
                  'x = ' + str(element[30:38]) +\
                  ', y = ' + str(element[38:46]) +\
                  ', z = '+ str(element[46:54]) +\
                  ', with rmsd peak of ' + str(element[54:60])

            print 'Nearest positive peak at: '+ str(minpospeakdist[counter1-1])[0:4]+ ' Angstrom'           
            print 'Keep this peak?'
            response2 = raw_input("Please enter something: ")
            
            if response2 in ('y','Y','yes','YES','Yes'):
                counter2 += 1
                
                newfileha.write(str(element).rstrip() + "\t" + str(counter1) + "\n")
                newfilepdb.write(str(arraypdb_coords[counter1-1]))
                
                print 'Peak kept in file!'
            else:
                print 'Peak removed from file!'
        newfilepdb.write("END")
        
    elif response1 in ('n','N','no','NO','No'):
        #if no manual searching selected, then no changes made to the difference peak files - renamed 'editted' to follow next parts of code
        for element in arrayha_coords:
            counter1 += 1
            counter2 += 1
            newfileha.write(str(element).rstrip() + "\t" + str(counter1) + "\n")
            newfilepdb.write(str(arraypdb_coords[counter1-1]))
        print 'All peaks chosen to be kept in file!'
        newfilepdb.write("END")
            
    print 'Peak analysis complete! Number of peaks remaining is: '+ str(counter2)
    
    newfileha.close()
    newfilepdb.close()
    
    return counter2,response1
    ###############################################################################