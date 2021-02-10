# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 12:31:22 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""

import os
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist, mahalanobis, euclidean
from scipy.stats import spearmanr#,kendalltau
from itertools import combinations 
import pickle



#%%############################################################################
#                                      functions
###############################################################################



def create_noise_signal(design_mat,NrOfVoxels,NrOfStimuli,chunks,sigma):
    
    pred = design_mat[:,0]
    
    
    #empty variable where we have to store the signal
    noise_stimuli = np.empty((NrOfVoxels,len(pred),NrOfStimuli))

    
    task_duration = 20
    #onset of the task in the prt
    task_ranges = np.arange(20,len(pred)-20,40)


    #separate the voxels in responding and not responding voxels
    #the voxel index will be first shuffled and then divided 
    #in chuncks of equal lenght. The number of chuncks is the
    #number of stimuli

    
    
    #%%############################################################################
    # adding noise across time for each signal (mean=0,sigma=1)
    ###############################################################################

    for i in range(NrOfStimuli):
        
        #pick the index of the voxels that will have a response to the stimulus
        voxels_used = chunks[i]
        #adding the expected signal to the voxel interested
        
        for j in range(NrOfVoxels):
            
            if j in voxels_used:
                noise_stimuli[j,:,i] = pred + np.random.normal(0,sigma,len(pred))

                [noise_stimuli[j,onset:onset+task_duration,i]*np.random.uniform(-1,3) 
                                                                for onset in task_ranges]
                [noise_stimuli[j,onset-20:onset,i]*np.random.uniform(-1,1.5) 
                                                for onset in task_ranges]
                # # # plt.plot(noise_stimuli[j,:,i])
                # plt.plot(pred)
            else:
                noise_stimuli[j,:,i] = np.random.normal(0,sigma,len(pred))
                # plt.plot(noise_stimuli[j,:,i])
                [noise_stimuli[j,onset:onset+task_duration,i]*np.random.uniform(-1,2) 
                                                for onset in task_ranges]
                [noise_stimuli[j,onset-20:onset,i]*np.random.uniform(-1,2) 
                                                                for onset in task_ranges]                
           

    
    return noise_stimuli



def incr_glm(stimulus, design_mat, tp):
    
    """" This function receives as input the signal with noise and the 
         timepoint of interest (timepoint at wich perform the GLM)
        It returns the corresponding t-values per voxel for that timepoint
    """
  
    #storing variable for the t-values
    tstats = np.empty((stimulus.shape[0]))
    all_betas = np.empty((stimulus.shape[0]))    
    all_res = np.empty((stimulus.shape[0]))

        
    #one constrast
    contrast = np.hstack((1,np.zeros(design_mat.shape[1]-1)))     
    #corresponding design matrix    
    X = design_mat[:tp,:]  

            
    ##%%###########################################################################
    # Statistics 
    ###############################################################################
                    
    #picking one voxel at time
    for idx_v,values in enumerate(stimulus):
           
        #glm for the voxel
        glm = sm.OLS(values[:tp],X).fit()
        
        #saving betas 
        all_betas[idx_v] = glm.params[0]
        betas = glm.params
    
        #saving the variance of the residuals for the next t-stats
        all_res[idx_v] = np.var(glm.resid)
        res = np.var(glm.resid)

        
        #%%############################################################################   
        #evaluate contrasts
        ###############################################################################
        
        #result of the matrix multiplication in the denominator of the tstats
        inv_X = np.linalg.inv(np.dot(X.T,X))

        #estimating numerator and denominator
        t_num = np.dot(contrast,betas)
        t_den = np.sqrt((res*(np.linalg.multi_dot((contrast,inv_X,contrast)))))
    
        #storing t-values
        tstats[idx_v] = t_num/t_den

                                   
    return tstats 


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
        
    elif dist_metric == 4:
        
        RDM = np.zeros((data.shape[0],data.shape[0]))
        
        #creating the idenx for the couple interested
        idx2select = list(combinations(range(0,data.shape[0]),2))
        
        for id in idx2select:
            
            #betas relative to two stimuli
            b1 = data[id[0],:]
            b2 = data[id[1],:]
            #pseudo-inverse beacuse the covariance is not invertible
            iv = np.linalg.inv(np.cov(np.vstack((b1,b2)).T))
            
            m_dist = mahalanobis(b1,b2,iv)
            RDM[id[0],id[1]] = m_dist
            RDM[id[1],id[0]] = m_dist
            

        

    return RDM
    


            

def rdm_corr(ref_rdm,rdm):
    
    """ This function estimates the rank-correlation between the given
        RDM and the given reference RDM
    """
    
    idx2select = list(combinations(range(0,rdm.shape[0]),2))
    #upper part selection
    
    #vector of the ref rdm
    ref_vec = np.triu(ref_rdm)
    ref_vec = np.array([ref_vec[id] for id in idx2select])
    
    #vec of the intermediate RDM
    tmp_rdm = np.triu(rdm)
    rdm_vec = np.array([tmp_rdm[id] for id in idx2select])
 
    #estimate rank correlation
    rnkcorr_vals = spearmanr(ref_vec,rdm_vec)[0] 
    rnkcorr_p_vals = spearmanr(ref_vec,rdm_vec)[1]
    
    dist_rdm = ref_vec-rdm_vec
    
    return dist_rdm, rnkcorr_vals, rnkcorr_p_vals 

    # return rnkcorr_vals, rnkcorr_p_vals            
            

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
        distances = [cdist(np.array(base_stim[:,idx]).reshape(1,-1),
                           np.array(new_stim).reshape(1,-1),'euclidean')[0][0] 
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
# MAIN 
##############################################################################

TBV_TargetFolder = 'C:/Users/andre/Documents/TBVData/rtRSA_5Dec2019/TBVFiles/run_0_cat'
ProjName = os.path.basename(TBV_TargetFolder)

#accessing full design matrix
#first accessing predictors
f = open(os.path.join(TBV_TargetFolder,'Prepared_Realtime-DesignMatrix.txt'),'r')
design_mat = np.genfromtxt(f)

design_mat = design_mat[:,:2]
#prolonging desing matrix
# design_mat = np.vstack((design_mat[:-10],design_mat[20:]))

print('Full design matrix loaded!')

         

tps = np.arange(40,430,40)
tps = np.append(tps,430)

NrOfVoxels = 100
NrOfStimuli = 5

NSim = 1000
 
rdm_corr_val = np.empty((NSim,len(tps)))
p_corr_val = np.empty((NSim,len(tps)))


idx2select = list(combinations(range(0,NrOfStimuli),2))


dist_rdms = np.empty((len(idx2select),len(tps),NSim))
proj_distances = np.empty((len(tps),NrOfStimuli,NSim))
# #the responding voxels for each stimuli are fixed in advance
# chunks = []

#variance of the noise of the created stimuli
sigma = 3

all_noise_stimuli = np.empty((NrOfVoxels,tps[-1],NrOfStimuli, NSim))
all_full_glm_tvals = np.empty((NrOfVoxels,NrOfStimuli, NSim))
all_glm_tvals = np.empty((NrOfVoxels,NrOfStimuli,len(tps), NSim))


for i in range(NSim):
    
    #the responding voxels for each stimuli are fixed in advance
    chunks = []
    #in this loop we are trying to recreate different spatial pattern
    #3 stimuli (0,2,4) should be much correlated compared to the other two
    #(1,3) that will show also less activation
    for k in range(NrOfStimuli):
        #stimuli with odd index will have 60 voxels activated
        if np.mod(k,2) == 0:
            # voxels_ind = np.arange(NrOfVoxels-30)
            voxels_ind = np.arange(NrOfVoxels)
            np.random.shuffle(voxels_ind)
            chunks.append(voxels_ind[:np.random.randint(30,60)])
        else:
            #even stimuli will have only 20 voxels activated
            # voxels_ind = np.arange(NrOfVoxels-30,NrOfVoxels)
            voxels_ind = np.arange(NrOfVoxels)
            np.random.shuffle(voxels_ind)
            chunks.append(voxels_ind[:np.random.randint(20,50)])   
        
    print('Simulation #',str(i))

    
    #create noisy stimuli both in variance and amplitude 
    noise_stimuli = create_noise_signal(design_mat,NrOfVoxels,NrOfStimuli,chunks,sigma)
    all_noise_stimuli[:,:,:,i] = noise_stimuli
    
    #full glm
    full_glm_tvals = np.array([incr_glm(noise_stimuli[:,:,i],design_mat,tps[-1]) 
                                  for i in range(NrOfStimuli)]).T
    all_full_glm_tvals[:,:,i] = full_glm_tvals    
    
    #create the reference rdm
    ref_rdm = first_order_rdm(full_glm_tvals,1)
    #estimating reference coordinates and space  
    ref_rsa_coords, ref_inv_mat = cmdscale(ref_rdm,2)
    # plt.plot(full_glm_tvals[:,0])
 
    

    #loop over the timepoints
    for idx_tp,tp in enumerate(tps):
        glm_tvals = np.array([incr_glm(noise_stimuli[:,:tp,j],design_mat,tp) 
                            for j in range(NrOfStimuli)]).T
        
        all_glm_tvals[:,:,idx_tp,i] = glm_tvals
        # plt.plot(glm_tvals[:,0])
        # print(np.corrcoef(glm_tvals[:,0],full_glm_tvals[:,0]))
        #print(tp)
        rdm = first_order_rdm(glm_tvals,1)
        
        for j in range(NrOfStimuli):
            distances = estimate_distances(full_glm_tvals, glm_tvals[:,j], 1) 
        
            #estimate the position of the new stimulus in the space
            stimulus_positions = np.array(target_positioning(distances, ref_inv_mat, ref_rdm))
            proj_distances[idx_tp,j,i] = euclidean(stimulus_positions,ref_rsa_coords[j,:])    
        
        # rdm_corr_val[i,idx_tp] , p_corr_val[i,idx_tp] = rdm_corr(ref_rdm,rdm)
        dist_rdms[:,idx_tp, i],rdm_corr_val[i,idx_tp] , p_corr_val[i,idx_tp] = rdm_corr(ref_rdm,rdm)  

#%%############################################################################
#   saving
        
f=open('simulated_stimuli_sigma_'+str(sigma)+'.pkl','wb')
pickle.dump(all_noise_stimuli,f)
f.close()

f=open('rdm_corr_val_'+str(sigma)+'.pkl','wb')
pickle.dump(rdm_corr_val,f)
f.close()

f=open('p_corr_val'+str(sigma)+'.pkl','wb')
pickle.dump(p_corr_val,f)
f.close()
  
f=open('dist_rdms'+str(sigma)+'.pkl','wb')
pickle.dump(dist_rdms,f)
f.close()
 
f=open('proj_distances'+str(sigma)+'.pkl','wb')
pickle.dump(proj_distances,f)
f.close() 

#%%############################################################################
#visualization
###############################################################################
        
        
plt.figure(figsize=(19.2,10.8))

#averaging over the comparisons
average_dist_rdms = np.average(dist_rdms,axis=0)       

err = np.sqrt(np.std(average_dist_rdms,axis=1)/(average_dist_rdms.shape[1]-1))
up = np.average(average_dist_rdms,axis=1)+err
down = np.average(average_dist_rdms,axis=1)-err      
plt.plot(np.arange(average_dist_rdms.shape[0]),np.average(average_dist_rdms,axis=1))
plt.fill_between(np.arange(average_dist_rdms.shape[0]),up,down,alpha=0.5, edgecolor='black')

plt.ylim = [-1,1]
plt.xticks(list(range(dist_rdms.shape[1])))

#plt.tight_layout()
plt.xlabel('Trials',size=20)
plt.ylabel('Average RDMs difference with standard error',size=20)

plt.title('Average RDMs difference over trials',size=50)

plt.savefig('Rank_corr_RDMs_spatial_pattern_average_stims_overlap_chuncks_sigma'+str(sigma)+'.png',figsize=(19.2,10.8),dpi=300) 

#%%
plt.figure(figsize=(19.2,10.8))

#averaging over the comparisons
#average_dist_rdms = np.average(dist_rdms,axis=0)       

err = np.sqrt(np.std(dist_rdms,axis=2)/(dist_rdms.shape[2]-1))
mean_dist = np.average(dist_rdms,axis=2)
up = mean_dist + err
down = mean_dist - err
for i in range(len(err)):      
    plt.plot(np.arange(dist_rdms.shape[1]),mean_dist[i,:])
    plt.fill_between(np.arange(mean_dist.shape[1]),up[i,:],down[i,:],
                     alpha=0.5, edgecolor='black')

plt.ylim = [-1,1]
plt.xticks(list(range(dist_rdms.shape[1])))

#plt.tight_layout()
plt.xlabel('Trials',size=20)
plt.ylabel('Average RDMs difference with standard error',size=20)

plt.title('Average RDMs difference over trials',size=50)

plt.savefig('Difference_in_RDMs_spatial_pattern_all_stims_overlap_chunks_sigma'+str(sigma)+'.png',figsize=(19.2,10.8),dpi=300) 


#%%

plt.figure(figsize=(19.2,10.8))
       
err = np.sqrt(np.std(rdm_corr_val,axis=0)/(rdm_corr_val.shape[0]-1))      
plt.plot(range(rdm_corr_val.shape[1]),np.average(rdm_corr_val,axis=0))
plt.fill_between(range(rdm_corr_val.shape[1]),
                 np.average(rdm_corr_val,axis=0)-err,
                 np.average(rdm_corr_val,axis=0)+err,
                 alpha=0.5, edgecolor='black')

plt.ylim = [-1,1]
plt.xticks(list(range(rdm_corr_val.shape[1])))

#plt.tight_layout()
plt.xlabel('Trials',size=20)
plt.ylabel('Rank correlation with standard error',size=20)

plt.title('RDMs correlation over trials',size=50)
plt.savefig('RankCorr_RDMs_spatial_pattern_all_stims_overlap_chunks_sigma'+str(sigma)+'.png',figsize=(19.2,10.8),dpi=300) 


#%%

plt.figure(figsize=(19.2,10.8))
       
err = np.sqrt(np.std(proj_distances,axis=2)/(proj_distances.shape[2]-1))
mean_dist = np.average(proj_distances,axis=2)
up = mean_dist + err
down = mean_dist - err
for i in range(err.shape[1]):      
    plt.plot(np.arange(proj_distances.shape[0]),mean_dist[:,i])
    plt.fill_between(np.arange(mean_dist.shape[0]),up[:,i],down[:,i],
                     alpha=0.5, edgecolor='black')

plt.ylim = [-1,1]
plt.xticks(list(range(proj_distances.shape[0])))

#plt.tight_layout()
plt.xlabel('Task-blocks',size=20)
plt.ylabel('Average distance from target and standard error',size=20)
plt.xticks(np.arange())

#plt.title('Average distance from target over trials',size=50)
plt.savefig('AvgDistTarget_RDMs_spatial_pattern_all_stims_overlap_chunks_sigma'+str(sigma)+'.png',figsize=(19.2,10.8),dpi=300) 
