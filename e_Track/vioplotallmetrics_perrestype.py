# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 02:35:33 2015

@author: charlie
"""

##########--------------##########--------------##########--------------##########
##########--------------##########--------------##########--------------##########


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_palette("deep", desat=.6)
sns.set_context(rc={"figure.figsize": (8, 4)})
np.random.seed(9221999)

sns.set(rc={"figure.figsize": (6, 6)})

#create full list of residues in structure:
PDBmulti.sort(key=lambda x: x.basetype)
res_list = []
for atom in PDBmulti:
    if atom.basetype not in (res_list) and atom.basetype != 'CA':
        res_list.append(atom.basetype)
        

###########
mainorside = ['sidechain','mainchain']
#create boxplots for each dataset number for residues in res_list above
counter = 0
for res_type in res_list:
    
    # Create a figure instance
    fig = plt.figure()

    dens_multilist_min,dens_multilist_max,dens_multilist_mean,dens_multilist_median = [],[],[],[]    
    for ctype in mainorside:
        for i in range(0,len(PDBmulti[0].Bfactor)):
            dens_list_min,dens_list_max,dens_list_mean,dens_list_median = [],[],[],[]
            for atom in PDBmulti:   
                if atom.basetype in (str(res_type)) and atom.side_or_main() in (ctype):
                    dens_list_min.append(atom.mindensity[i])
                    dens_list_max.append(atom.maxdensity[i])
                    dens_list_mean.append(atom.meandensity[i])
                    dens_list_median.append(atom.mediandensity[i])

            dens_multilist_min.append(dens_list_min)
            dens_multilist_max.append(dens_list_max)
            dens_multilist_mean.append(dens_list_mean)
            dens_multilist_median.append(dens_list_median)

        counter += 1
    
    # Create an axes instance
    ax = fig.add_subplot(221)    
    sns.violinplot(dens_multilist_min, color="coolwarm_r", lw=2);
    X_labels = ['xfel\nside','Sync 4.1A\nside','Sync 3.6A\nside',
                'xfel\nmain','Sync 4.1A\nmain','Sync 3.6A\nmain']
    ## Custom x-axis labels
    ax.set_xticklabels(X_labels)  
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()   
    ## Custom title,x-label,y-label    
    ax.set_title(str(res_type)+':min',
                 fontsize=20)                
    plt.xlabel('Dataset', fontsize=18)
    plt.ylabel('Electron Density', fontsize=16)
    
        # Create an axes instance
    ax = fig.add_subplot(222)    
    sns.violinplot(dens_multilist_max, color="coolwarm_r", lw=2);
    X_labels = ['xfel\nside','Sync 4.1A\nside','Sync 3.6A\nside',
                'xfel\nmain','Sync 4.1A\nmain','Sync 3.6A\nmain']
    ## Custom x-axis labels
    ax.set_xticklabels(X_labels)  
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()   
    ## Custom title,x-label,y-label    
    ax.set_title(str(res_type)+':max',
                 fontsize=20)                
    plt.xlabel('Dataset', fontsize=18)
    plt.ylabel('Electron Density', fontsize=16)
    
        # Create an axes instance
    ax = fig.add_subplot(223)    
    sns.violinplot(dens_multilist_mean, color="coolwarm_r", lw=2);
    X_labels = ['xfel\nside','Sync 4.1A\nside','Sync 3.6A\nside',
                'xfel\nmain','Sync 4.1A\nmain','Sync 3.6A\nmain']
    ## Custom x-axis labels
    ax.set_xticklabels(X_labels)  
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()   
    ## Custom title,x-label,y-label    
    ax.set_title(str(res_type)+':mean',
                 fontsize=20)                
    plt.xlabel('Dataset', fontsize=18)
    plt.ylabel('Electron Density', fontsize=16)
    
        # Create an axes instance
    ax = fig.add_subplot(224)    
    sns.violinplot(dens_multilist_median, color="coolwarm_r", lw=2);
    X_labels = ['xfel\nside','Sync 4.1A\nside','Sync 3.6A\nside',
                'xfel\nmain','Sync 4.1A\nmain','Sync 3.6A\nmain']
    ## Custom x-axis labels
    ax.set_xticklabels(X_labels)  
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()   
    ## Custom title,x-label,y-label    
    ax.set_title(str(res_type)+':median',
                 fontsize=20)                
    plt.xlabel('Dataset', fontsize=18)
    plt.ylabel('Electron Density', fontsize=16)
        
   

    ## Save the figure
    fig.set_size_inches(20,20)
    fig.savefig('./plots/combineddatasetsMAINSIDETOGETHER/' + 
                    'OliPro' + '_'+str(res_type) + '_violin_.png',bbox_inches='tight')
        
        