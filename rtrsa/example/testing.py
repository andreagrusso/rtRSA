# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 15:58:08 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""
import os
import glob
from rtrsa.nfrsa import rtRSA
from rtrsa import utils
import numpy as np
from expyriment_stash.extras.expyriment_io_extras import tbvnetworkinterface


#%%
#working direcrtory
wdir = os.getcwd()

single_maps = os.path.join(wdir,'single_tmaps')
merged_maps =  os.path.join(wdir,'merged_tmaps')
config_file = os.path.join(wdir,'config.json')

#%%
#extract tvalues. This command needs to used for each run
utils.TBV_value_extractor('localhost',0,0,single_maps,'run-NF')

#%%
#load raw tvalues ("raw" because the array where 
#they are stored have different lengths )

raw_tvals_files = glob.glob(os.path.join(single_maps,'*.tvals'))

raw_tvals = [np.loadtxt(open(f,'r')) for f in raw_tvals_files]
    
cc=utils.intersect_coords(raw_tvals)

#%%
utils.merge_tmaps('test','pearson',2,single_maps,merged_maps,'class_test')


#%%
newRSA = rtRSA('ff',4,'euclidean')
newRSA.load(config_file)

#%%

nf_tvals_f = glob.glob(os.path.join(single_maps,'*NF*.tvals'))
nf_tvals = np.loadtxt(open(nf_tvals_f[0],'r'))

#set of coordinates to use for the extraction during the NF run
new_nf_coords = newRSA.match_coords(nf_tvals)[:,:-1]

#test target positioning
TBV = tbvnetworkinterface.TbvNetworkInterface('localhost',55555)    

#get new tvalues
new_stim = [TBV.get_map_value_of_voxel(0,coords)[0] 
                                    for coords in new_nf_coords]
#get coordinates in the RS for the new tvalues
x,y = newRSA.target_positioning(new_stim)

