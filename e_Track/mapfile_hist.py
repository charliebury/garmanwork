# -*- coding: utf-8 -*-

from map2voxelclasslist_nocartconversion import densmap2class_readheader,densmap2class_readvoxels
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math
import struct
from classholder_v2 import electron_map_info
import os
import sys

def go():
	mapfilename = '../2BN3damage_11feb2014/2BN3damage2_density.map'
	maptype = 'density_map'
	atom_indices = []

	rho,densitystart = densmap2class_readheader(mapfilename)

	vxl_list = densmap2class_readvoxels(mapfilename,densitystart)

	rho.vxls_val = vxl_list

	return rho

def maphist():

	rho = go()

	# histogram plot of number of voxels per atom
	print '\n------------------------------------------------'
	print 'Plotting histogram of density value per voxel...' 
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (10, 6)})
	fig = plt.figure()

	datax = [rho.vxls_val[i] for i in range(0,len(rho.vxls_val))]

	plt.hist(datax, 300, histtype="stepfilled", alpha=.7);
	plt.xlabel('Density per voxel')
	plt.ylabel('Frequency')
	plt.title('Histrogram of density per voxel')
	fig.savefig('denspervxl_hist.png')
    

def densmap2class_MOD(mapfilename,maptype,atom_indices):
    
    # define 'rho' electron map object
    rho = electron_map_info()

    # open electron density .map file here 
    binarymapfile = open(mapfilename,'rb')
    
    # start adding header information into electron_map_info class format. 
    # Note the unpacking of a struct for each byte, read as a long 'l'
    rho.nx = struct.unpack('=l',binarymapfile.read(4))[0]
    rho.ny = struct.unpack('=l',binarymapfile.read(4))[0]
    rho.nz = struct.unpack('=l',binarymapfile.read(4))[0]
    rho.type = struct.unpack('=l',binarymapfile.read(4))[0]
    print 'Num. Col, Row, Sec: '
    print '%s %s %s' %(rho.nx,rho.ny,rho.nz)
    
    rho.start1 = struct.unpack('=l',binarymapfile.read(4))[0] 
    rho.start2 = struct.unpack('=l',binarymapfile.read(4))[0] 
    rho.start3 = struct.unpack('=l',binarymapfile.read(4))[0] 
    print 'Start positions: '
    print '%s %s %s' %(rho.start1,rho.start2,rho.start3)

    rho.gridsamp1 = struct.unpack('=l',binarymapfile.read(4))[0] 
    rho.gridsamp2 = struct.unpack('=l',binarymapfile.read(4))[0] 
    rho.gridsamp3 = struct.unpack('=l',binarymapfile.read(4))[0]
    print 'Grid sampling:'
    print '%s %s %s' %(rho.gridsamp1,rho.gridsamp2,rho.gridsamp3)

    # for cell dimensions, stored in header file as float not long 
    # integer so must account for this
    rho.celldim_a = struct.unpack('f',binarymapfile.read(4))[0]
    rho.celldim_b = struct.unpack('f',binarymapfile.read(4))[0]
    rho.celldim_c = struct.unpack('f',binarymapfile.read(4))[0]
    rho.celldim_alpha = struct.unpack('f',binarymapfile.read(4))[0]
    rho.celldim_beta = struct.unpack('f',binarymapfile.read(4))[0]
    rho.celldim_gamma = struct.unpack('f',binarymapfile.read(4))[0]
    print 'Cell dimensions:'
    print '%s %s %s' %(rho.celldim_a,rho.celldim_b,rho.celldim_c)
    print '%s %s %s' %(rho.celldim_alpha,rho.celldim_beta,rho.celldim_gamma)


    rho.fast_axis = struct.unpack('=l',binarymapfile.read(4))[0] 
    rho.med_axis = struct.unpack('=l',binarymapfile.read(4))[0] 
    rho.slow_axis = struct.unpack('=l',binarymapfile.read(4))[0] 
    print 'Fast,med,slow axes: '
    print '%s %s %s' %(rho.fast_axis,rho.med_axis,rho.slow_axis)

    rho.mindensity = struct.unpack('f',binarymapfile.read(4))[0] 
    rho.maxdensity = struct.unpack('f',binarymapfile.read(4))[0] 
    rho.meandensity = struct.unpack('f',binarymapfile.read(4))[0]
    print 'Density values: '
    print '%s %s %s' %(rho.mindensity,rho.maxdensity,rho.meandensity)

    # next find .map file size, to calculate the last nx*ny*nz bytes of 
    # file (corresponding to the position of the 3D electron density 
    # array). Note factor of 4 is included since 4-byte floats used for 
    # electron density array values.
    filesize = os.path.getsize(mapfilename)
    densitystart = filesize - 4*(rho.nx*rho.ny*rho.nz)
    
    # next seek start of electron density data
    binarymapfile.seek(densitystart,0)   
    
    # if electron density written in floats (which is to be expected 
    # from FFT-CCP4 outputted .map file of electron density)
    if rho.type is 2:  
        struct_fmt = '=f4'
        struct_len = struct.calcsize(struct_fmt)
        density = []
        appenddens = density.append
        
        if maptype in ('atom_map'):
            atom_indices = []
            nonatom_indices = []
            appendindex = atom_indices.append
            appendnonindex = nonatom_indices.append
            counter = -1
            while True:
                data = binarymapfile.read(struct_len)
                if not data: break
                s = struct.unpack(struct_fmt,data)[0]
                counter += 1
                if int(s) == 0:
                	appendnonindex(counter)
                else:    
                    appenddens(s)
                    appendindex(counter)
        
        # efficient way to read through density map file using indices of atoms
        # from atom map file above
        elif maptype in ('density_map'):
            for i in range(0,len(atom_indices)):
                if i != 0:
                    binarymapfile.read(struct_len*(atom_indices[i]-atom_indices[i-1] - 1))
                else:
                    binarymapfile.read(struct_len*(atom_indices[0]))
                    
                data = binarymapfile.read(struct_len)
                s = struct.unpack(struct_fmt,data)[0]
                appenddens(s)
                
        else:
            print 'Unknown map type --> terminating script'
            sys.exit()
                    
    binarymapfile.close()
    
    # as a check that file has been read correctly, check that the min 
    # and max electron density stated in .map file header correspond to 
    # calculated min and max here
    # Note for the case of atom map, the min voxel val will be 0 (no atom present)
    # --> since these voxels are removed during the filtering of the map, only 
    # the map value is tested.
    # For density map, cannot currently perform a check, since there is no 
    # guarantee that the max and min voxel values may be non atom voxels and
    # thus removed
    if maptype in ('atom_map'):        
        if max(density) == rho.maxdensity:
            print 'calculated max voxel value match value stated in file header'
        else:
            print 'calculated max voxel value does NOT match value stated in file header'
            print 'have now calculated max voxel value to be: %s'\
            %str(max(density))
            sys.exit()
    
    # if each voxel value is an atom number, then want to convert to integer
    if maptype in ('atom_map'):
        density_final = [int(dens)/100 for dens in density]
    elif maptype in ('density_map'):
        density_final = density
    else:
        print 'Unknown map type --> terminating script'
        sys.exit()

    rho.vxls_val = density_final
    
    if maptype in ('atom_map'):
        return rho,atom_indices,nonatom_indices
    else:
        return rho


def go2():
	mapfilename1 = '../2BN3damage_12feb2014/2BN3damage10_density.map'
	mapfilename2 = '../2BN3damage_12feb2014/2BN3damage10_atoms.map'
	maptype1 = 'density_map'
	maptype2 = 'atom_map'

	# first want to determine the location of the atoms and nonatom voxels in atom map
	atom_indices = []
	rho1,atom_indices,nonatom_indices = densmap2class_MOD(mapfilename2,maptype2,atom_indices)

	# next for the density map, can determine the density associated with atoms only
	rho2 = densmap2class_MOD(mapfilename1,maptype1,atom_indices)

	# next for the density map, can determine the density associated with non-atom voxels only
	rho3 = densmap2class_MOD(mapfilename1,maptype1,nonatom_indices)

	# can now determine list of atom and non-atom voxel densities
	atm_vxlvals = rho2.vxls_val
	nonatm_vxlvals = rho3.vxls_val

	# histogram plot of density per atom and per non-atom
	print '\n------------------------------------------------'
	print 'Plotting histogram of density value per voxel...'

	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (10, 6)})
	fig = plt.figure()

	data1 = atm_vxlvals
	data2 = nonatm_vxlvals

	sns.kdeplot(np.array(data1), color="#6495ED", alpha=.5,shade=True)
	sns.kdeplot(np.array(data2), color="#F08080", alpha=.5,shade=True)

	plt.xlabel('Density per voxel')
	plt.ylabel('Frequency')
	plt.title('Histrogram of density per voxel')
	fig.savefig('denspervxl_hist.png')


