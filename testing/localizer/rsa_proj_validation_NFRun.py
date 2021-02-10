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
import matplotlib.lines as mlines
from scipy.spatial.distance import euclidean
from scipy.stats import spearmanr
from itertools import combinations 

#%%############################################################################
# Functions
###############################################################################

def estimate_distances(base_stim,new_stim, method):
    
    """ This function estimate the distances between the t-values relative 
    to the base stimuli (estimated externally) and the t-values of the new
    stimulus.
    
    Input:
        base_stim = matrix containing the t-values of the stimuli used to 
        estimate the embedding space (VxN, where V is the number of voxels and
        N is the number of conditions/stimuli)
        new_stim = vector that contains the t-values relative to the current
        stimulus (1xV, where V is the number of voxels)
        method = (1,1-Pearson),(2,Euclidean distance),(3, Mean activity difference)
        
    Output:
        distances = 1xN vector containing the distances (estimated with the 
        proposed method) between the current stimulus and the base stimuli

    """
    
    #the number of stimuli are the columns of the matrix
    NrOfBaseStim = np.shape(base_stim)[1]

    if method == 1:
        # Use the Pearson correlation and transform it in a distance
        distances = [1-np.corrcoef(base_stim[:,idx], new_stim)[0][1] 
        for idx in range(NrOfBaseStim)]

    if method == 2:
        # Use the classical euclidean distance
        distances = [cdist(base_stim[:,idx],new_stim,'euclidean') 
        for idx in range(NrOfBaseStim)]

    if method == 3:
        # Use absolute activation difference
        distances = [np.mean(base_stim[:,idx])-np.mean(new_stim) 
        for idx in range(NrOfBaseStim)]

    
    return distances


def target_positioning (distances, inv_mat, base_rdm):
    
    """ This function uses the distances between the base stimuli
        and the new one, the projection matrix and the RDM calculated on
        the base stimuli (both estimated externally) to project the new
        stimulus in the RSA embedding space defined by the base stimuli
        
        Input:
            distances = 1xN vector that encodes the distances between the t-values
            of the current stimulus and the t-values of the base stimuli
            inv_mat = matrix that allow to project the new stimulus in the 
            embedded space
            base_rdm = the RDM matrix calculated on the t-values of the base 
            stimuli
        Output:
            x,y coordinates of the new stimulus in the embedded space
        
        """
    
    #square the distances of the new stimuli
    sq_dist = np.square(distances)

    #compute mean of the square distances
    rdm_mean = np.sum(np.square(base_rdm),axis=1)/len(base_rdm)
    
    #compute the actual coords in the embedding space 
    coords = -0.5*inv_mat.dot((sq_dist - rdm_mean).T)

    return coords[0],coords[1]   



#%%############################################################################
#                                INPUTS                                       #
###############################################################################


rdms = pickle.load(open(os.path.join('output','rdm_files.pkl'),'rb'))
rsa_coords = pickle.load(open(os.path.join('output','rsa_coords_files.pkl'),'rb'))
inv_mats = pickle.load(open(os.path.join('output','inv_mats_files.pkl'),'rb'))

#the key of the three dictionaries are the identical
data_name = list(rdms.keys())


input_folders = input('Insert the names of the folders, separated by a comma, \
                      where are the data:\n \
                      (e.g. NFDog,...): \n')
                
input_folders = input_folders.split(',')

#the name of the data are constant across the folders
input_data = os.listdir(input_folders[0])


localizer_folders = input('Insert the names of the localizer folders, separated by a comma, \
                      where are the data:\n \
                      (e.g. cat,dog,...): \n')
                
localizer_folders = localizer_folders.split(',')


target_stim = input('Insert the names of the target stimulus:\n \
                      (e.g. cat,dog,...): \n')
                      
idx_target = localizer_folders.index(target_stim)

#%%############################################################################
#                                Reference data                               #
###############################################################################

#data relative to the reference pattern are the last in every dictionary
ref_rdm = rdms[data_name[-1]]
ref_coords = rsa_coords[data_name[-1]]
ref_inv_mat = inv_mats[data_name[-1]]




#selecting the tvalues realtive to the last timepoint
base_stim = np.array([pickle.load(open(os.path.join(folder,data_name[-1]),'rb'))
                  for folder in localizer_folders]).T

#%%############################################################################
# Main loop
###############################################################################

#creating a variable to store the positions at the intermediate timepoints
stimulus_positions = np.empty((len(input_data),ref_coords.shape[1]))

#loop on timepoint excluding the last that is used as reference
for i,tp in enumerate(input_data):
    
    tvals = np.array([pickle.load(open(os.path.join(folder,tp),'rb'))
                      for folder in input_folders])
    values = tvals[0]
    distances = estimate_distances(base_stim, values, 1)
    
    #estimate the position of the new stimulus in the space
    stimulus_positions[i,:] = target_positioning(distances, ref_inv_mat, ref_rdm)        


#%%############################################################################
# Visualization
###############################################################################


plt.ion()
fig = plt.figure()

plt.figure(facecolor='dimgray')
ax = plt.axes()
ax.scatter(ref_coords[:,0],ref_coords[:,1],s=150, c='yellow',edgecolors = 'black')
# Setting the background color
ax.set_facecolor('dimgray')
plt.xticks([])
plt.yticks([])
plt.axis('off')
for label, x, y in zip(localizer_folders, ref_coords[:,0],ref_coords[:,1]):
    plt.annotate(label, xy=(x, y), xytext=(-10, 10),size=15,
                 textcoords='offset points', ha='right', va='bottom',
                 bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=1),
                 arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0')) 


yellow_circle = mlines.Line2D([], [], color='yellow',  linestyle='None', marker='o',markeredgecolor='black', markeredgewidth=0.5,
                          markersize=15, label='Base sitmuli')
red_star = mlines.Line2D([], [], color='red',linestyle='None',  marker='*',markeredgecolor='black', markeredgewidth=0.5,
                          markersize=15, label='Current stimulus')

plt.legend(handles=[yellow_circle, red_star],
    loc='upper center', bbox_to_anchor=(0.5, -0.05),ncol=4)




sent = 'Which stimulus do you want to see: ', input_folders, ': \n'
stimulus_name = input(sent)


stimulus_coords = stimulus_positions

for i,coords in enumerate(stimulus_coords[:,:]):
    
    if i == 0:
        plt.scatter(coords[0],coords[1], 
                    marker = '*',s=200, color = 'red', edgecolors='black')
        ax.set_facecolor('dimgray')
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
    else:
        plt.scatter(stimulus_coords[:i+1,0],stimulus_coords[:i+1,1], 
                    marker = '*',s=200, color = 'darkgray')
        plt.scatter(coords[0],coords[1], 
                    marker = '*',s=200, color = 'red', edgecolors='black')
        for label, x, y in zip(range(i),stimulus_coords[:i+1,0],stimulus_coords[:i+1,1]):
            plt.annotate(str(label+1),xy=(x, y), xytext=(5,-5),textcoords='offset points')        
        #plotting the trajectory
        plt.plot(stimulus_coords[:i+1,0],stimulus_coords[:i+1,1], '-',
                    color = 'green')
        ax.set_facecolor('dimgray')
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
    
    #save figure
    plt.savefig(os.path.join(stimulus_name+ '_trial'+str(i+1)+'.png'),
           facecolor='dimgray', edgecolor='none', dpi=600,
           tight_layout=True,figsize=((19.2,10.8)))


#%%############################################################################
#plot distances for every trial after projection
###############################################################################
    
#storing variable NtrialsXNrStimuli
distances = np.empty(len(stimulus_positions))

for i,coords in enumerate(stimulus_positions):
    
    #distances between pairs of two matrices
    distances[i] = euclidean(coords,ref_coords[idx_target,:])    



plt.figure('Distances from reference point in reference space',figsize=((19.2,10.8)))
plt.plot(np.append(np.roll(distances,1),distances[-1]),linewidth=4,color='b')
#plt.legend(input_folders)
#plt.title('Distances from reference point in reference space')
plt.ylabel('Distances',size=40)
plt.xlabel('Task-blocks',size=40)
plt.xlim(1,len(distances)+1)
xticks = [str(i) for i in np.arange(1,len(distances),1)]
xticks.append('Full \n timeseries')
plt.xticks(np.arange(1,len(distances)+1,1),xticks,size=40)
plt.yticks(size=40)
plt.savefig(os.path.join('output','NFRun_distances_in_ref_space.png'),dpi=600,
            tight_layout=True,figsize=((19.2,10.8)))


