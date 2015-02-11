# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 13:13:30 2015

@author: lina2532
"""
import sys 


###----------------###----------------###----------------###----------------###
###############################################################################      
# A class for .map file voxel
class voxel_density:
    def __init__(self,atmnum=0,density=0): 
        self.atmnum = atmnum
        self.density = density
###############################################################################
###----------------###----------------###----------------###----------------###



def combinevxl_atmanddensvals(atmmap_list,densmap_list):
    # create list of voxel objects in class voxel_density 
    
    vxl_list = []
    appendvxl = vxl_list.append
    
    atmmap_len = len(atmmap_list)
    # unessential loading bar add-in (part1)
    point = atmmap_len / 100
    increment = atmmap_len / 100
    
    for i in range(0,atmmap_len):
        appendvxl(voxel_density(atmmap_list[i],densmap_list[i]))

        # unessential loading bar add-in (part2)
        if(i % (1 * point) == 0):
            sys.stdout.write(
                    "\r[" + "=" * (i / increment) +
                    " " * ((atmmap_len - i)/ increment) + 
                    "]" +  str(i / point) + "%")
            sys.stdout.flush()
        
    return vxl_list
            
            
            
def combinevxl_atmanddensvals_gen(atmmap_list,densmap_list):
    # create list of voxel objects in class voxel_density 
        
    atmmap_len = len(atmmap_list)
    
    for i in range(0,atmmap_len):
         yield voxel_density(atmmap_list[i],densmap_list[i])

