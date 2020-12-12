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
baseline_duration, task_duration = 20,20


block_duration = baseline_duration + task_duration

onset = np.array([1,1 + baseline_duration])
offset = onset + baseline_duration-1

init_onset = onset
init_offset = offset

n_rip = 10

for i in range(1,n_rip):
    onset = np.append(onset,init_onset+i*block_duration)
    offset = np.append(offset,init_offset+i*block_duration)

#additional last baseline
onset = np.append(onset,offset[-1]+1)
offset = np.append(offset,onset[-1]+19)




baseline_index = np.arange(0,len(onset),2)
baselines = np.vstack((onset[baseline_index],offset[baseline_index])).T
task_index = np.arange(1,len(onset),2)
tasks = np.vstack((onset[task_index],offset[task_index])).T

feedbacks_w_baseline = np.arange(tasks[0,0]+5,baselines[-1,0],5)
feedbacks_wo_baseline = np.array([np.arange(tasks[i,0]+5,tasks[i,1],5)
                         for i in range(len(tasks)-1)]).reshape(-1,1)

#%%############################################################################

conditions = ['cat','chair','dog','hammer']

#non incremental prt
for im in conditions:
    im_name = im.split('.')[0]    
    protocol = prt.StimulationProtocol(experiment_name="RSA Pilot 1st run " + im_name, 
                                       time_units="Volumes")
    #baseline
    protocol.add_condition(prt.Condition("Baseline", baselines,colour=[0, 0, 255]))
    #task                                     
    for i in range(tasks.shape[0]):
        protocol.add_condition(prt.Condition(im_name + '_' + str(i), tasks[i,:].reshape((1,2)),colour=[255, i, 0]))

        
    protocol.save(os.path.join(outdir, im_name +'_RSA_individual_task.prt'))  

#incremental prt    
for im in conditions:
    im_name = im.split('.')[0]    
    protocol = prt.StimulationProtocol(experiment_name="RSA Pilot 1st run " + im_name, 
                                       time_units="Volumes")
    #baselines
    protocol.add_condition(prt.Condition("Baseline", baselines,colour=[0, 0, 255]))

    #task                                     
    protocol.add_condition(prt.Condition("Task", tasks ,colour=[255, 0, 0]))

        
    protocol.save(os.path.join(outdir, im_name +'_RSA_incremental_task.prt')) 

    
f=open(os.path.join(outdir,'new_paradigm','baselines.pkl'),'wb')
pickle.dump(baselines,f)


f=open(os.path.join(outdir,'new_paradigm','tasks.pkl'),'wb')
pickle.dump(tasks,f)

f=open(os.path.join(outdir,'new_paradigm','feedbacks_with_baseline.pkl'),'wb')
pickle.dump(feedbacks_w_baseline,f)
f=open(os.path.join(outdir,'new_paradigm','feedbacks_without_baseline.pkl'),'wb')
pickle.dump(feedbacks_wo_baseline,f)