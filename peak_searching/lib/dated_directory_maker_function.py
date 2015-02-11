# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 22:11:09 2014

@author: charlie
"""
import datetime
import os

###############################################################################
#MAKE A NEW DIRECTORY TO SAVE NEW FILES TO DURING CODE
#OUTPUT CURRENT DATE AND PATH TO DIRECTORY IN WHICH FILES WILL BE SAVED FOR RUN

class directorymaker:
    pass

def dated_directory_maker(inputfile):
    result = directorymaker()
    #Want to read a input txt file into script which contains all the key files names to be read
    inputlines = open(str(inputfile), "r").readlines() 
    for line in inputlines:
        if 'datasettitle' in str(line.split()[0]):
            datasettitle = str(line.split()[1])
        else:
            break   
    result.date =  datetime.datetime.now().strftime("%Y%m%d_%H%M%S") 
    result.path = 'outputfiles/'+str(datasettitle)+str('_')+str(result.date)

    os.makedirs(result.path)
    return result
###############################################################################

