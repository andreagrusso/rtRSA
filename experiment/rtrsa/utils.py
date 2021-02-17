# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:17:37 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""

import os
import numpy as np
from expyriment_stash.extras.expyriment_io_extras import tbvnetworkinterface
from rtrsa.nfrsa import rtRSA
import glob
import numpy_indexed as npi
import matplotlib.pyplot as plt
import json

#%%


def TBV_value_extractor(TBVip,voi,ctr,outdir,basename):
    
    """
    This function extracts the tvalues from data loaded in TBV and save them
    as binary files. 
    Up to this date there is no possibility to extract the beta values.

    Parameters
    ----------
    TBVip : string
        IP index to access the TBV processed data. When TBV is running on the
        same machine please use 'localhost'
    voi : integer
        A value ranging from 0 to infinite. It indicates the index of the ROI
        used in analysis.
    ctr : integer
        A value ranging from 0 to infinite. It indicates the index of the 
        contrast of interest. The contrast can be defined manually using a .ctr
        file or it can be defined automatically by TBV. Usually the contrast
        '0' corresponds to the map 'first predictors vs. baseline'.
    outdir : string
        The directory for the output files.
    basename : string
        Basename for the outputs.

    Returns
    -------
    None.

    """


    TBV = tbvnetworkinterface.TbvNetworkInterface(TBVip,55555)    
    
                   
    if TBV.get_current_time_point()[0] == TBV.get_expected_nr_of_time_points()[0]: 
        print('Extracting t-values...')
         
    #coordinates of voxels of the roi
        coord_roi_voxels = np.array(TBV.get_all_coords_of_voxels_of_roi(voi)[0])
        
        tvals = np.array([TBV.get_map_value_of_voxel(voi,coord)[0] 
                        for coord in coord_roi_voxels]).reshape(-1,1)

      
        output_tval = np.concatenate((coord_roi_voxels,tvals),axis=1)
            
        print('Saving t-vals data...')
        with open(os.path.join(outdir,basename+'_voi'+str(voi)+'.tvals'), 'w') as outfile:
            np.savetxt(outfile, output_tval)
        
        print('Everything have been estimated! Goodbye!')
    
#%%    
   

def intersect_coords(raw_tvals):
    
    
    """
    This function is needed to intersect the set of coordinates beloging
    to different base stimuli and different version of the same ROI
    

    Parameters
    ----------
    raw_tvals : TYPE
        A list of coordinates and corresponding tvalues. Each item of the
        list is a numpy matrix where the first three columns are the fucntional
        coordinates and the last column the tvalue.

    Returns
    -------
    cc : Numpy matrix
        The set of common coordinates between the givens sets.

    """


    for i in range(len(raw_tvals)-1):
        if i == 0:
            cc = npi.intersection(raw_tvals[i][:,:-1], raw_tvals[i+1][:,:-1])            
        else:            
            cc = npi.intersection(cc,raw_tvals[i+1][:,:-1])
            
    return cc


#%%
def merge_tmaps(name,dist_metric,n_comp,inputdir,outdir,basename):
    
    """
    This function takes the single t-maps associated to the base stimuli and 
    first, it combines them into an RDM, then it estimates both the Representational
    Space and the corresponding inversion matrix. 
    The main ouptut is an rtRSA object that contains all the data that is saved
    as a binary file.

    Parameters
    ----------
    name : string
        Name of the rtRSA object
    
    dist_metric : string
        'pearson'--> estimate the distance as 1-Pearson correlation.
        'euclidean' --> estimate the distance using the formula of the 
                        euclidean distance.
        'mean_diff'--> estimate the distance as the mean activation
                        difference between the brain patterns.
    n_comp : integer
        Number of dimension of the representational space. It is usually 2.
    inputdir : string
        Directory were the single t-maps are stored.
    outdir : string
        Directory where the output data are saved.
    basename : string
        Basename for the output.
        

    Returns
    -------
    None.

    """
    
    #first create an rtRSA object
    rtRSAObj = rtRSA(name,n_comp,dist_metric)
    
    #glob return the full path of the file
    tmaps = glob.glob(os.path.join(inputdir,'*.tvals'))
    tmaps.sort()
    
    #list of the names of the base stimuli
    conditions = []
        
    #tvals contains the common coordinates and the selected tvalues for
    #the corresponding base stimuli
    tvals = dict()
    
    #raw tvals is a list and it can contains t-values matrices of different
    #lenghts
    raw_tvals = []
        
    #stimuli is loaded in alphabetical order
    for filename in tmaps:
        conditions.append(os.path.basename(filename).split('_')[0])
        raw_tvals.append(np.loadtxt(open(filename,'r')))
    
    #load condition names    
    rtRSAObj.load_conditions(conditions)
    
    #intersect coordinates of the maps
    #we know how many base stimuli we have bu we cannot code directly
    #the intersection. Therefore we can implement a recursive technique
    
    #set of commmon cooordinates
    cc = intersect_coords(raw_tvals)
    #sort by 3rd column, 2nd column, 1st column
    sorted_cc = cc[np.lexsort((cc[:,2], cc[:,1],cc[:,0]))]
    
    #storing the common coords (functional space) in the rtRSA object
    rtRSAObj.load_func_coords(sorted_cc)
    
    
    tvals_mat = np.empty((len(cc),len(raw_tvals)))

    #use this common set of coordinates to extract the tvalues
    for i,data in enumerate(raw_tvals):
        
        #get the indices of the common coords
        for j,coord in enumerate(data[:,:-1]):
            index = np.where((sorted_cc == coord).all(axis=1))[0]
            #get the tvalues corresponding to the common coords
            tvals_mat[index,i] = data[j,-1]
            
    
    tvals['ROI'] = tvals_mat 
    
    #storing the t values of the base stimuli in the rtRSA object
    rtRSAObj.load_base_stimuli(tvals_mat)
    
    #create the RS
    rdm = rtRSAObj.createRS(tvals_mat)
    
    #write on disk all the results
    class_outdir = os.path.join(outdir,basename)
    os.mkdir(class_outdir)
    rtRSAObj.saveAs(class_outdir,basename)
    
   
    
    
    plt.figure()
    plt.title('RDM')
    plt.imshow(rdm,cmap='RdBu_r')
    plt.colorbar()
    plt.xticks(range(len(conditions)),conditions)
    plt.yticks(range(len(conditions)),conditions)
    plt.savefig(os.path.join(outdir,basename + '_RDM.png'))
    
    
    plt.figure(facecolor='darkslategray')
    ax = plt.axes()
    plt.title('RSA space',size=20)
    ax.scatter(rtRSAObj.RS_coords[:,0],rtRSAObj.RS_coords[:,1],
               s=200, c='gold',edgecolors = 'black')
    # Setting the background color
    ax.set_facecolor('grey')
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')
    for label, x, y in zip(conditions,rtRSAObj.RS_coords[:,0],rtRSAObj.RS_coords[:,1]):
        plt.annotate(label, xy=(x, y),size=15,
                     textcoords='offset points',xytext=(20, -20), ha='right', va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5', fc='gold', alpha=1),
                     arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0')) 
    plt.savefig(os.path.join(outdir,basename + '_RSA_coords.png'),facecolor='darkslategray', 
                edgecolor='none', dpi=256)

    
    
    
        
    