# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 21:49:44 2014

@author: charlie
"""

###############################################################################
#a simple function 'file_len' is created here to calculate the number of lines
#in a file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
###############################################################################
    