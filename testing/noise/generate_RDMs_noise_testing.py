# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 14:17:33 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""

import os
import pickle
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from itertools import combinations 

#%%############################################################################
# Functions
###############################################################################

def cmdscale(D,n_comp):
    """                                                                                       
    Classical multidimensional scaling (MDS)                                                  
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    D : (n, n) array                                                                          
        Symmetric distance matrix.

    n : int
        Number of components requested
                                                            
                                                                                               
    Returns                                                                                   
    -------                                                                                   
    Y : (n, p) array                                                                          
        Configuration matrix. Each column represents a dimension. Only the                    
        p dimensions corresponding to positive eigenvalues of B are returned.                 
        Note that each dimension is only determined up to an overall sign,                    
        corresponding to a reflection.                                                        
                                                                                               
    e : (n,) array                                                                            
        Eigenvalues of B.                                                                     
                                                                                               
    """
    # Number of points                                                                        
    n = len(D)
 
    # Centering matrix                                                                        
    H = np.eye(n) - np.ones((n, n))/n
 
    # YY^T                                                                                    
    B = -H.dot(D**2).dot(H)/2
 
    # Diagonalize                                                                             
    evals, evecs = np.linalg.eigh(B)
 
    # Sort by eigenvalue in descending order                                                  
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]
 
    #computing the matrix only for the positive-eigenvalued components
    w, = np.where(evals > 0)
    L  = np.diag(np.sqrt(evals[w]))
    V  = evecs[:,w]
    Y  = V.dot(L)
    
    
    
    #transformation matrix for only the first 2 components    
    coords = Y[:,:n_comp]
    inv_mat = evecs[:,:n_comp].T/(np.sqrt(evals[:n_comp]))[:,None]

    return coords, inv_mat



def first_order_rdm(data, dist_metric):
    """ 
    Return a Representation Dissimilarity Matrix estimated with a specific
    metric.
    Input:
        data = matrix of the values extracted from the ROI. The format accepted
        is NxV (N condition and V voxels)
        dist_metric = a number that specifies the metric to use
                    1 = Pearson correlation
                    2 = Euclidian distance
                    3 = Absolute activation difference
    Output:
        RDM = a NxN matrix where the dissimilarity between each condition is 
        stored
    """

    data = data.T

    if dist_metric == 1:
        # Use correlation distance
        RDM = 1-np.corrcoef(data)

    elif dist_metric == 2:
        # Use Eucledian distance
        RDM = cdist(data,data,'euclidean')

    elif dist_metric == 3:
        # Use absolute activation difference
        means = np.mean(data,axis=1) # Determine mean activation per condition
        m, n = np.meshgrid(means,means) # Create all possible combinations
        RDM = abs(m-n) # Calculate difference for each combination

    return RDM
    




#%%############################################################################
#                                INPUTS                                       #
###############################################################################

#change thi
input_folders = input('Insert the names of the folders, separated by a folder, \
                      where are the data:\n \
                      (e.g. cat,dog,...): \n')
                
input_folders = input_folders.split(',')

#the name of the data are constant across the folders
input_data = os.listdir(input_folders[0])

###############################################################################
#                                OUTPUTS                                      #
###############################################################################

outdir = 'rdms'


#%%############################################################################
# Main loop
###############################################################################


noise_levels = np.arange(0.1,1.2,0.2)



#loop on noise level
plt.figure(figsize=(19.2,10.8))
plt.title('RDMs rank correlation for different noise levels')

for k,noise in enumerate(noise_levels):

    #a dictionary to store the RDM corresponding to each timepoint
    rdms = np.zeros((len(input_folders),len(input_folders),len(input_data)))
    #a dictionary to store the RSA coordinates correspondig to each timepoints
    rsa_coords = np.empty((len(input_folders),len(input_folders)-2,len(input_data)))
    #a dictionary to store the inversion matrix corresponding to each timepoints
    inv_mats = np.empty((len(input_folders)-2,len(input_folders),len(input_data)))  
    
    #loop on each data (related to a specific timepoint)
    for idx_tp,tp in enumerate(input_data):
        

        noisy_tvals = np.empty((100,len(noise_levels),len(input_folders)))

        for i,folder in enumerate(input_folders):
            noisy_tvals[:,:,i] = np.array(pickle.load(open(os.path.join(folder,tp),'rb')))
        

        
        tvals = noisy_tvals[:,k,:]
        
        #selecting one level of noise for all the stimuli
        
        #estimate RDM
        rdms[:,:,idx_tp] = first_order_rdm(tvals,1) 
        
        #transform in 2D space
        rsa_coords[:,:,idx_tp], inv_mats[:,:,idx_tp] = cmdscale(rdms[:,:,idx_tp],2)
        
        
#%% visualization for the current noise level
    rnkcorr_vals = np.empty((len(input_data[:-1])))
    rnkcorr_p_vals = np.empty((len(input_data[:-1])))
    
    ref_rdm = rdms[:,:,-1]
    
    #index of the upper part of the RDM. The ranksum correlation needs only the upper
    #or the lower part the RDM which is symmetric
    idx2select = list(combinations(range(0,len(input_folders)),2))
    
    #vector of the ref rdm
    ref_vec = np.triu(ref_rdm)
    ref_vec = np.array([ref_vec[id] for id in idx2select])
        
    #loop across the key of the dictionary
    for i in range(len(input_data[:-1])):
        
        #upper part selection
        tmp_rdm = np.triu(rdms[:,:,i])
        rdm_vec = np.array([tmp_rdm[id] for id in idx2select])
     
        #estimate rank correlation
        rnkcorr_vals[i] = spearmanr(ref_vec,rdm_vec)[0] 
        rnkcorr_p_vals[i] = spearmanr(ref_vec,rdm_vec)[1] 
    
    
    subplotid = '32'+str(k+1)
    plt.subplot(int(subplotid))
    plt.plot(rnkcorr_vals,linewidth=3)
    plt.xticks(range(len(input_data[:-1])))
    plt.title('noise level '+ str(noise)[:3])
    plt.ylim((-1,1))
    plt.xlabel('Trials')
    plt.ylabel('Rank correlation')
    plt.tight_layout()

        
        

    pickle.dump(rdms,open(os.path.join(outdir,'rdm_noise_level'+str(k)+'.pkl'),'wb'))
    pickle.dump(rsa_coords,open(os.path.join(outdir,'rsa_coords_noise_level'+str(k)+'.pkl'),'wb'))
    pickle.dump(inv_mats,open(os.path.join(outdir,'inv_mats_noise_level'+str(k)+'.pkl'),'wb'))

# #%%############################################################################
# # Reciprocal distances evaluation
# ###############################################################################

plt.savefig(os.path.join('rank_corr_between_RDMs_across_trials.png'),
            figsize=(19.8,10.2),dpi=256)
    












 

    