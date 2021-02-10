# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 08:50:42 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""

import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import scipy as sp
from scipy.stats import combine_pvalues

directory = 'C:/Users/andre/_NeuroImaging/scripts/maastricht/pyRSA-BV/testing/noise_testing'

sigma = 3

f=open(os.path.join(directory,'proj_distances'+str(sigma)+'.pkl'),'rb')
proj_distances = pickle.load(f)
f.close() 

f=open(os.path.join(directory,'rdm_corr_val_'+str(sigma)+'.pkl'),'rb')
rdm_corr_val = pickle.load(f)
f.close() 

f=open(os.path.join(directory,'p_corr_val'+str(sigma)+'.pkl'),'rb')
p_corr_val = pickle.load(f)
f.close() 

rdm_mean = np.average(rdm_corr_val,axis=0)
p_mean = combine_pvalues(p_corr_val[:,0],method='fisher')#(rdm_corr_val,axis=0)


#%% fitting exponential to all the simulations output

fitted_data = np.empty((rdm_corr_val.shape[0],rdm_corr_val.shape[1]))
r_sq = np.empty((rdm_corr_val.shape[0],1))

#loop on repetitions
for idx,vec_y in enumerate(rdm_corr_val):
    print(idx)
    
    vec_x = np.arange(len(vec_y))


    slope, intercept, r_sq[idx], p_value, std_err = sp.stats.linregress(vec_x, vec_y)
    fitted_data[idx,:] = slope*np.exp(vec_x) + intercept    

#%% plotting
plt.figure(figsize=(19.2,10.8))
       
err = np.sqrt(np.std(rdm_corr_val,axis=0)/(rdm_corr_val.shape[0]-1)) 
data_avg = np.average(rdm_corr_val,axis=0)
neg_err = data_avg-err
pos_err = data_avg+err
     
plt.plot(np.append(np.roll(data_avg,1),data_avg[-1]))
plt.fill_between(range(pos_err.shape[0]+1),
                  np.append(np.roll(neg_err,1),neg_err[-1]),
                  np.append(np.roll(pos_err,1),pos_err[-1]),
                  alpha=0.5, edgecolor='black')

plt.ylim = [-1,1]
plt.xlim(1,len(data_avg)+1)
xticks = [str(i) for i in np.arange(1,len(data_avg),1)]
xticks.append('Full \n timeseries')
plt.xticks(np.arange(1,len(data_avg)+1,1),xticks,size=40)
#plt.tight_layout()
plt.xlabel('Task-blocks',size=30)
plt.ylabel('Rank correlation with standard error',size=30)
#plt.xticks(np.arange(rdm_corr_val.shape[1]),size=20)
plt.yticks(size=40)

# plt.title('RDMs correlation over trials',size=50)
namesave = 'RankCorr_RDMs_spatial_pattern_all_stims_overlap_chunks_sigma'+str(sigma)+'.png'
plt.savefig(os.path.join(directory,namesave),tight_layout=True,figsize=(19.2,10.8),dpi=300) 

#%%


plt.figure(figsize=(19.2,10.8))
       
err = np.sqrt(np.std(proj_distances,axis=2)/(proj_distances.shape[2]-1))
mean_dist = np.average(proj_distances,axis=2)
up = mean_dist + err
down = mean_dist - err
for i in range(err.shape[1]):
    tmp_dist = mean_dist[:,i]
    tmp_up = up[:,i]
    tmp_down = down[:,i]      
    plt.plot(np.append(np.roll(tmp_dist,1),tmp_dist[-1]))
    plt.fill_between(np.arange(mean_dist.shape[0]+1),
                     np.append(np.roll(tmp_up,1),tmp_up[-1]),
                     np.append(np.roll(tmp_down,1),tmp_down[-1]),
                     alpha=0.5, edgecolor='black')

    plt.ylim = [-1,1]
    plt.ylim = [-1,1]
    plt.xlim(1,len(tmp_dist)+1)
    xticks = [str(i) for i in np.arange(1,len(tmp_dist),1)]
    xticks.append('Full \n timeseries')
    plt.xticks(np.arange(1,len(tmp_dist)+1,1),xticks,size=40)
    plt.xlabel('Task-blocks',size=30)
    plt.ylabel('Average distance from target and standard error',size=30)
    plt.yticks(size=40)
    plt.legend(['Stimulus #1','Stimulus #2','Stimulus #3','Stimulus #4','Stimulus #5'],fontsize=20)
    #plt.title('Average distance from target over trials',size=50)
    namesave = 'AvgDistTarget_RDMs_spatial_pattern_all_stims_overlap_chunks_sigma'+str(sigma)+'.png'
    plt.savefig(os.path.join(directory,namesave),tight_layout=True,figsize=(19.2,10.8),dpi=300) 
