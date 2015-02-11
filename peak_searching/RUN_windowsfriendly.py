# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 17:23:19 2014

@author: charlie
"""
###############################################################################
#!/usr/bin/env python
from math import sqrt, pow
from lib.file_len_function import file_len
from lib.dated_directory_maker_function import dated_directory_maker
from lib.termcolor import colored
from lib.positive_diff_peak_remover import positive_diff_peak_remover
from lib.threshold_peak_remover import threshold_peak_remover
from lib.manual_peak_searching import manual_peak_searching
from lib.write_outputfiles import write_outputfiles
from lib.threshold_peak_remover_VDVradius import threshold_peak_remover_VDVradius
import sys
import os
from lib.negpeak_closepospeak import negpeak_closepospeak
from lib.PDBtoCLASSARRAY import *
from lib.glossy_output import *
from lib.graph_analysis41damageset import peakfiltering_bystage,barplot_atompropfreq,damage_wrtatmnum 
from lib.savevariables import save_objectlist


######################################################################################################################################################
#SECTION: INITIAL SET UP FOR CODE, INCLUDING DIRECTORY MANAGEMENT AND READING input.txt FILE
######################################################################################################################################################

# Title screen here
title_page()

##MAKE A NEW DIRECTORY TO SAVE NEW FILES TO DURING CODE
inputfile = './inputfiles/inputfile.txt'
datedir = dated_directory_maker(inputfile)

#want to make a log file containing info of outputs of difference map peak processing
logfile = open(datedir.path +'/'+"logfile"+datedir.date+".txt", "w")
logfile.write(datedir.date)

#Want to read a input txt file into script which contains all the key files names to be read
inputfilein = open(str(inputfile), "r")
inputlines = inputfilein.readlines()

counter = 0;
numinputvariables = 5;
for line in inputlines:
    if (('#' or '_') == str(line[0])):
        pass
    if 'datasettitle' in str(line.split()[0]):
        runtitle = str(line.split()[1])
    elif 'threshold' in str(line.split()[0]): 
        counter += 1
        print colored('Damage searching threshold chosen to be ' + str(line.split()[1]) +' Angstrom','yellow')
        logfile.write('\nDamage searching threshold chosen to be ' + str(line.split()[1]) +' Angstrom')
        threshold = str(line.split()[1])
    elif 'pdb_file' in str(line.split()[0]):
        counter += 1
        print colored('Structure pdb entry chosen to be ' + str(line.split()[1]),'magenta')
        logfile.write('\nStructure pdb entry chosen to be ' + str(line.split()[1]))
        proDNApdb = str(line.split()[1])
    elif 'diffpeak_pdb' in str(line.split()[0]):
        counter += 1
        print colored('pdb file containing difference map peaks chosen to be ' + str(line.split()[1]),'cyan')
        logfile.write('\npdb file containing difference map peaks chosen to be ' + str(line.split()[1]))        
        pdbfilename = str(line.split()[1])
    elif 'diffpeak_ha' in str(line.split()[0]):
        counter += 1
        print colored('ha file containing difference map peaks chosen to be ' + str(line.split()[1]),'blue')
        logfile.write('\nha file containing difference map peaks chosen to be ' + str(line.split()[1]))        
        hafilename = str(line.split()[1])
    elif 'end' in str(line.split()[0]) and counter == numinputvariables:
        print 'End of input file'
        break
    
inputs_where = './inputfiles/'   
filelist = [inputs_where+proDNApdb,inputs_where+pdbfilename,inputs_where+hafilename]

#new section here 
filler_section()

###############################################################################
#specify here the PDB structure used (not the pdb file with damage sites)
pdbin = open(filelist[0], "r")
#write pdb atom coordinate lines to a txt file for easy access later
pdbcoords = open(datedir.path + '/' + "pdbcoords.txt", "w")
lines = pdbin.readlines()
PDBarray=[]

#detemine which lines of .pdb file are atom coords and these are written to
#pdbcoords file 
for line in lines:
    if 'ATOM' in line.split()[0]:
        pdbcoords.write(line)
        PDBarray.append(str(line))
    else: 
        pass
pdbin.close()
pdbcoords.close()
###############################################################################

# create a list to detail the number of difference map peaks filtered at each stage of analysis
num_peaks_kept_list = []

# first determine the initial number of unfiltered difference map peaks present in .pdb diff peak file
counter = 0
for line in open(str(filelist[1]), "r" ).readlines():
    if 'ATOM' in line[0:4]:
        counter += 1
print 'Initial number of peaks before processing is: ' + str(counter)
num_peaks_kept_list.append(counter)

# trigger response from user before any processing begins:
print 'Ready to begin filtering process?\nWill first filter to remove positive peaks'
welcome = raw_input('press enter..\n')


###############################################################################################################################################################
##UPDATE:26NOV2014: ADDED NEW FUNCTION HERE TO READ PDB STRUCTURE TO A CLASS FOR EASY ACCESSIBILITY
##this inclusion can replace the above code in which the pdb structure is read line by line into an array for use later
#PDBarray2 = PDBtoCLASSARRAY_fixedlinewidth(proDNApdb)
##############################################################################################################################################################

#new section here 
filler_section()

######################################################################################################################################################
#SECTION: REMOVE POSITIVE DIFFERENCE MAP PEAKS
num_neg_peaks = positive_diff_peak_remover(filelist,datedir)
logfile.write('\nNumber of neg peaks after positive_diff_peak_remover function: ' + str(num_neg_peaks))
num_peaks_kept_list.append(num_neg_peaks)
######################################################################################################################################################


#new section here 
filler_section()

##############################################################################################################################################################
#SECTION: REMOVE PEAKS FURTHER THAN A SPECIFIED DISTANCE THRESHOLD FROM ATOMS IN THE STRUCTURE PDB FILE
#num_kept_peaks = threshold_peak_remover(threshold,datedir)
print 'Do you want to threshold with preset threshold specified in input file'+colored('(type:1)','yellow')+', or use Van der Waals atomic radii to threshold '+colored('(type:2)','magenta')+'?'
thresholdtype = raw_input("Please enter something: (1/2) ")

if thresholdtype == '1':
    num_kept_peaks = threshold_peak_remover(threshold,datedir)
    logfile.write('\nInput.txt file global threshold used')
elif thresholdtype == '2':
    logfile.write('\nAtom VDV-radius-specific threshold used')
    num_kept_peaks = threshold_peak_remover_VDVradius(threshold,datedir,PDBarray)
else:
    print colored('ERROR CHOOSING THRESHOLDING TYPE','red')
logfile.write('\nNumber of neg peaks after threshold_peak_remover function: ' + str(num_kept_peaks))
num_peaks_kept_list.append(num_kept_peaks)
##############################################################################################################################################################


#new section here 
filler_section()

##############################################################################################################################################################
#SECTION: LOCATE CLOSEST POSITIVE PEAKS TO EACH REMAINING NEGATIVE PEAK AND OUTPUT LIST OF LOCATIONS -> TO BE USED AS TOOL TO ADD MANUAL FILTERING OF REMAINING PEAKS
minpospeakdist,sidechainlist = negpeak_closepospeak(filelist,datedir,PDBarray)
##############################################################################################################################################################


#new section here 
filler_section()

##############################################################################################################################################################
#SECTION: MANUAL SEARCHING THROUGH REMAINING DIFFERENCE MAP PEAKS TO REMOVE UNWANTED PEAKS NOT CONSIDERED TO BE ASSOCIATED WITH SPECIFIC DAMAGE
peaks_remaining,response1 = manual_peak_searching(datedir,minpospeakdist)
if response1 in ('y','Y','yes','YES','Yes'):
    logfile.write('\nManual searching chosen')
    logfile.write('\nNumber of neg peaks after manual_peak_searching function: ' + str(peaks_remaining))
else:
    logfile.write('\nManual searching not chosen')
##############################################################################################################################################################


#new section here 
filler_section()


##############################################################################################################################################################
#SECTION: Write new .txt file containing processed information - listing specific atom details from the structure associated with each selected specific damage sites
peakobjlist = write_outputfiles(datedir,PDBarray)
logfile.write('\nScript finished and ' + str(peaks_remaining) +' sites of specific damage located')
num_peaks_kept_list.append(peaks_remaining)
##############################################################################################################################################################

logfile.close()


##############################################################################################################################################################
# SECTION: Output some graphs to detail the reduction in peaks during the filtering process,
# and also statistics on damage sites per residue type
peakfiltering_bystage(datedir,num_peaks_kept_list)
barplot_atompropfreq(datedir,peakobjlist,'basetype')
barplot_atompropfreq(datedir,peakobjlist,'atomtype')
barplot_atompropfreq(datedir,peakobjlist,'chaintype')
damage_wrtatmnum(datedir,peakobjlist,PDBarray)
##############################################################################################################################################################



##############################################################################################################################################################
# SECTION: Use cPickle to save the list of objects 'peakobjlist' to a file
save_objectlist(peakobjlist,runtitle+'_'+proDNApdb)
##############################################################################################################################################################

