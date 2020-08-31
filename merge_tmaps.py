# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:20:30 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""


from rtrsa import utils
from tkinter import filedialog
import os, glob

name = input('Insert name of the newly created rtRSA object:   ')

metric = input('Insert name of the metric to use for the RSA (pearson, euclidean, mean_diff):   ')

if not metric: 
    metric ='pearson'


single_maps = filedialog.askdirectory(title='Select the directory where the single t-values maps are')
merged_maps = filedialog.askdirectory(title='Select the output directory')

max_comp = len(glob.glob(os.path.join(single_maps,'*.tvals')))-1

n_comp = int(input('Insert the number of dimension of the RS (min 2, max '+ str(max_comp) + '):   '))

if not n_comp: 
    metric = 2



out_name = input('Insert name of the output:   ')


utils.merge_tmaps(name,metric,n_comp,single_maps,merged_maps,out_name)

