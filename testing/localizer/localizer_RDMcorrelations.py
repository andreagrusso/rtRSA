# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 11:02:41 2021

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""

import pickle
import os
import numpy as np
from scipy.stats import spearmanr
from itertools import combinations
import matplotlib.pyplot as plt

wdir = 'C:/Users/andre/_NeuroImaging/scripts/maastricht/pyRSA-BV/Testing/output'

rdms = pickle.load(open(os.path.join(wdir,'rdm_files.pkl'),'rb'))

rnkcorr_vals = np.empty((len(rdms.keys())-1))
rnkcorr_p_vals = np.empty((len(rdms.keys())-1))
    
ref_rdm = rdms['tvals_tp430.pkl']
    
    #index of the upper part of the RDM. The ranksum correlation needs only the upper
    #or the lower part the RDM which is symmetric
idx2select = list(combinations(range(0,4),2))
    
    #vector of the ref rdm
ref_vec = np.triu(ref_rdm)
ref_vec = np.array([ref_vec[id] for id in idx2select])
      
keys = list(rdms.keys())
#loop across the key of the dictionary
for i in range(len(rdms.keys())-1):
    
    #upper part selection
    tmp_rdm = np.triu(rdms[keys[i]])
    rdm_vec = np.array([tmp_rdm[id] for id in idx2select])
 
    #estimate rank correlation
    rnkcorr_vals[i] = spearmanr(ref_vec,rdm_vec)[0] 
    rnkcorr_p_vals[i] = spearmanr(ref_vec,rdm_vec)[1] 
    

plt.plot(rnkcorr_vals)