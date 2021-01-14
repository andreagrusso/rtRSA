# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:27:04 2019

@author: Andrea Gerardo Russo
Biomedical Engineer, PhD candidate in Neuroscience
University of Salerno, Italy
"""

import os
import numpy as np 
from brainvoyagertools import prt
import pickle


#%%#############################PATHS##########################################

wdir = os.getcwd()

outdir = os.path.join(wdir,'prt')
settings = dict()

#%%############################################################################

#first block limits
baseline_duration = 20
task_duration = 20


block_duration = baseline_duration + task_duration

baselines = np.array([1, baseline_duration]).reshape(1,-1)
tasks = np.array([baseline_duration+1,baseline_duration + task_duration]).reshape(1,-1)



n_rip = 10

for i in range(1,n_rip):
    baselines = np.vstack((baselines,baselines[0,:]+i*block_duration))
    tasks = np.vstack((tasks,tasks[0,:]+i*block_duration))

#additional last baseline
baselines = np.vstack((baselines,baselines[-1,:]+block_duration))



feedbacks = np.array([np.arange(tasks[i,0]+7,tasks[i,1],2)
                         for i in range(len(tasks))]).reshape(-1,1)

#%%############################################################################

conditions = ['gatto','sedia','cane','martello']

  

#incremental prt    
for im in conditions:
    im_name = im.split('.')[0]    
    protocol = prt.StimulationProtocol(experiment_name="rtRSA " + im_name, 
                                       time_units="Volumes")
    #baselines
    protocol.add_condition(prt.Condition("Baseline", baselines,colour=[0, 0, 255]))

    #task                                     
    protocol.add_condition(prt.Condition("Task", tasks ,colour=[255, 0, 0]))

        
    protocol.save(os.path.join(outdir, 'continuous', im_name +'_RSA_incremental_task.prt')) 

    
f=open(os.path.join(outdir,'continuous','baselines.pkl'),'wb')
pickle.dump(baselines,f)


f=open(os.path.join(outdir,'continuous','tasks.pkl'),'wb')
pickle.dump(tasks,f)


f=open(os.path.join(outdir,'continuous','feedbacks.pkl'),'wb')
pickle.dump(feedbacks,f)


# #non incremental prt
# for im in conditions:
#     im_name = im.split('.')[0]    
#     protocol = prt.StimulationProtocol(experiment_name="RSA Pilot 1st run " + im_name, 
#                                        time_units="Volumes")
#     #baseline
#     len_baseline = baselines.shape[0]
#     baselines[0,:] = np.array([1,8])
#     baselines = np.insert(baselines,2,np.array([11,20])).reshape(len_baseline+1,2)
#     protocol.add_condition(prt.Condition("Baseline", baselines,colour=[0, 0, 255]))
#     #task                                     
#     for i in range(tasks.shape[0]):
#         tmp_task = np.vstack((np.array([9,10]),tasks[i,:]))
#         protocol.add_condition(prt.Condition(im_name + '_' + str(i), tmp_task, colour=[255, i, 0]))

        
#     protocol.save(os.path.join(outdir,'continuous', im_name +'_RSA_individual_task.prt'))